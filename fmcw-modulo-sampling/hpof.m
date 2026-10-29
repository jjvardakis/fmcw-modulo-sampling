function [rI, rQ, info] = hpof(zI, zQ, lam, g, phi)
%HPOF  Hilbert-Pencil-of-Function unfolding (Alg. 1, Feuillen & Bhandari, SSP 2025).
%   [rI, rQ, info] = hpof(zI, zQ, lam)            paper version (ideal IQ)
%   [rI, rQ, info] = hpof(zI, zQ, lam, g, phi)    IQ-imbalance-aware version (g, phi known from calibration)
%
%   zI, zQ : N x 1 modulo samples of the I and Q channels (N even)
%   lam    : modulo threshold; only used to snap the fold heights to the lattice 2*lam*Z
%   g, phi : Q-branch gain error and phase error,  Q' = (1+g)*(rQ*cos(phi) + rI*sin(phi))
%   rI, rQ : N x 1 unfolded channels (r = z + eps)
%   info   : struct with intermediate quantities for plotting (g-signal, singular values, folds)
if nargin < 4, g = 0; phi = 0; end
zI = zI(:);  zQ = zQ(:);  N = numel(zI);       % N x 1 each
L  = floor(N/3) - 1;                           % pencil parameter (paper: N/3 - 1)
k  = (0:N-1).';                                % N x 1 DFT bin index
s  = sign(k - N*(k > N/2));                    % N x 1 sign(w_k), same as in dht.m

% ---- Step 1 (Eq. 3): remove the unknown HDR signal with Lemma 1 ---------------------------
% ideal IQ:  d = zI - H(zQ) = H(epsQ) - epsI          (the signal r cancels exactly)
% imbalance: d = cos(phi) zI + sin(phi) H(zI) - H(zQ)/(1+g)   (also cancels r exactly)
d = cos(phi)*zI + sin(phi)*dht(zI) - dht(zQ)/(1+g);           % N x 1

% ---- Step 2 (Eq. 4-5): first difference -> spikes at the fold instants; DFT --------------
gd = circshift(d, -1) - d;                     % N x 1, gd(n) = d(n+1) - d(n) (circular)
G  = fft(gd);                                  % N x 1 spectrum = sum of exponentials in k (poles = fold times)
yp = G(1:N/2+1);                               % (N/2+1) x 1  : w >= 0  (y+)
ym = G(N/2+2:N);                               % (N/2-1) x 1  : w <  0  (y-)

% ---- Step 3: stacked Hankel matrix (both halves share the same poles) -------------------
Yp = hankel(yp(1:end-L), yp(end-L:end));       % (N/2+1-L) x (L+1)
Ym = hankel(ym(1:end-L), ym(end-L:end));       % (N/2-1-L) x (L+1)
Y  = [Yp; Ym];                                 % (N-2L) x (L+1), e.g. 88 x 85 for N = 256
[~, S, V] = svd(Y, 'econ');                    % V: (L+1) x (L+1) right singular vectors
sv = diag(S);                                  % (L+1) x 1, sorted descending
[~, K] = max(log10(sv(1:end-1)+realmin) - log10(sv(2:end)+realmin));  % model order = biggest SV gap

% ---- Step 4: matrix pencil -> poles p_m = exp(-1i*2*pi*m/N) -> fold instants m ----------
W = conj(V(:, 1:K));                           % (L+1) x K signal subspace (conj: V' rows -> columns)
p = eig(pinv(W(1:end-1,:))*W(2:end,:));        % K x 1 poles; pinv((L) x K) is K x L, times L x K -> K x K
m = unique(mod(round(-angle(p)*N/(2*pi)), N)); % K' x 1 fold instants (0-based), rounded to the sample grid

% ---- Step 5: fold heights by least squares on the whole spectrum; snap to 2*lam*Z -------
A  = exp(-1i*2*pi*k*m.'/N);                    % N x K': (N x 1)*(1 x K') outer product; column = one fold time
MI = -(cos(phi) + 1i*s*sin(phi)).*A;           % N x K' spectral signature of an I fold (ideal: -A)
MQ = (1i*s/(1+g)).*A;                          % N x K' spectral signature of a Q fold (ideal: +1i*sign(w)*A)
Mx = [MI, MQ];                                 % N x 2K'
x  = [real(Mx); imag(Mx)] \ [real(G); imag(G)];% 2K' x 1 real unknowns [cI; cQ];  (2N x 2K') \ (2N x 1)
x  = 2*lam*round(x/(2*lam));                   % fold heights live on the lattice 2*lam*Z

% ---- Step 6: integrate the steps back into residues (anti-difference S) -----------------
nm = numel(m);
dI = zeros(N,1);  dI(m+1) = x(1:nm);           % N x 1 Delta epsI (spikes), m+1: 0-based -> MATLAB index
dQ = zeros(N,1);  dQ(m+1) = x(nm+1:end);       % N x 1 Delta epsQ
eI = [0; cumsum(dI(1:end-1))];                 % N x 1 residue, up to a constant
eQ = [0; cumsum(dQ(1:end-1))];
eI = eI + 2*lam*round(-mean(zI + eI)/(2*lam)); % constant fixed by zero-mean r (no target at DC)
eQ = eQ + 2*lam*round(-mean(zQ + eQ)/(2*lam));
rI = zI + eI;                                  % N x 1  r = z + eps
rQ = zQ + eQ;                                  % N x 1  (= Q' if an imbalance was given)

info.d = d;  info.g = gd;  info.sv = sv;  info.K = K;
info.folds = find(dI ~= 0 | dQ ~= 0) - 1;      % detected fold instants (0-based)
info.dI = dI;  info.dQ = dQ;
info.candidate_folds = m; % distinct rounded locations before height fitting/rounding
end
