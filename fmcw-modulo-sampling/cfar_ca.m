function hit = cfar_ca(P, Nt, Ng, Pfa)
%CFAR_CA  Cell-averaging CFAR on a power spectrum (circular in range).
%   P   : N x 1 power |X|^2 of the range profile
%   Nt  : training cells per side,  Ng : guard cells per side,  Pfa : design false-alarm probability
%   hit : N x 1 logical, true where a target is declared
%   Threshold = alpha * mean(training cells), alpha = 2Nt*(Pfa^(-1/(2Nt)) - 1)  (exponential noise power)
N     = numel(P);
alpha = 2*Nt*(Pfa^(-1/(2*Nt)) - 1);
lags  = [-(Ng+Nt):-(Ng+1), (Ng+1):(Ng+Nt)];    % 1 x 2Nt offsets of the training cells
hit   = false(N,1);
for n = 1:N
    idx    = mod(n-1 + lags, N) + 1;           % 1 x 2Nt circular indices
    hit(n) = P(n) > alpha*mean(P(idx));
end
end
