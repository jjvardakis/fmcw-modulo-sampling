%% EXPERIMENT 1 - Conventional (clipping) ADC vs modulo ADC + Hilbert-PoF, strong + weak target
%  Noiseless, ideal IQ: isolates what clipping/folding do to the signal and to target estimation.
close all; setup_scene; rng(1);
quickRun = exist('v2_quick','var') && v2_quick;

% ---- 1. analog IQ beat signal --------------------------------------------------------------
a  = amp.*exp(1i*2*pi*rand(1,2));           % 1 x 2 complex amplitudes (random scattering phase)
rC = fmcw_beat(R, a, fc, B, Tc, N);          % N x 1 true HDR signal
rI = real(rC);  rQ = imag(rC);               % N x 1 each

% ---- 2. the two ADCs (same input range +-lam) ----------------------------------------------
cI = clip(rI);         cQ = clip(rQ);        % N x 1 conventional ADC output
zI = modlam(rI, lam);  zQ = modlam(rQ, lam); % N x 1 modulo ADC output

% ---- 3. unfolding --------------------------------------------------------------------------
[hI, hQ, info] = hpof(zI, zQ, lam);          % Hilbert-PoF (paper)
uI = itoh(zI);         uQ = itoh(zQ);        % classic unwrapping baseline (success depends on sample differences)

% ground-truth folds, for the printout only
epsI = round((rI - zI)/(2*lam));  epsQ = round((rQ - zQ)/(2*lam));   % N x 1 integer residues
MC   = union(find(circshift(epsI,-1) ~= epsI), find(circshift(epsQ,-1) ~= epsQ));
fprintf('Peak |r| = %.2f V vs ADC range %.1f V  (x%.1f)\n', max(abs([rI; rQ])), lam, max(abs([rI; rQ]))/lam);
fprintf('Folds: |M_C| = %d (theoretical reference N/3 = %.1f), detected by HPoF = %d\n', numel(MC), N/3, numel(info.folds));
fprintf('NMSE  clipped ADC        : %7.1f dB\n', nmse(cI + 1i*cQ, rC));
fprintf('NMSE  modulo + unwrapping: %7.1f dB\n', nmse(uI + 1i*uQ, rC));
fprintf('NMSE  modulo + HPoF      : %7.1f dB\n\n', nmse(hI + 1i*hQ, rC));

% ---- 4. target estimation from the range profile -------------------------------------------
sig  = {rC, cI + 1i*cQ, zI + 1i*zQ, hI + 1i*hQ};          % 1 x 4 cell, each N x 1
name = {'ideal ADC', 'clipped ADC', 'modulo, not unfolded', 'modulo + HPoF'};
fprintf('%-22s %10s %10s %12s %12s %16s\n', 'receiver', 'R1 [m]', 'R2 [m]', '|a1| [V]', '|a2| [V]', 'worst ghost [dB]');
for i = 1:4
    X  = abs(fft(sig{i}))/N;                 % N x 1 amplitude spectrum (on-grid tone of amp a -> |a|)
    Xp = X(2:N/2);                           % (N/2-1) x 1 positive bins 1..N/2-1
    [~, ord] = sort(Xp, 'descend');
    pk = sort(ord(1:2));                     % 2 strongest bins (1-based in Xp -> bin index = pk)
    ghost = X;  ghost(kbin+1) = 0;           % everything that is not a true target
    fprintf('%-22s %10.1f %10.1f %12.4f %12.4f %16.1f\n', name{i}, pk(1)*dR, pk(2)*dR, ...
            X(pk(1)+1), X(pk(2)+1), 20*log10(max(ghost)/abs(a(2))));
end
fprintf('%-22s %10.1f %10.1f %12.4f %12.4f   (truth; ghost relative to the weak target)\n\n', 'TRUE', R, abs(a));

% ---- 5. robustness check: Lemma 1 is exact only for on-grid (integer-bin) tones ------------
offs = [0 0; 0 0.5; 0.1 0; 0.3 0.3];         % 4 x 2 range offsets in bins (strong, weak)
for i = 1:size(offs,1)
    r2 = fmcw_beat(R + offs(i,:)*dR, a, fc, B, Tc, N);
    [h2I, h2Q] = hpof(modlam(real(r2), lam), modlam(imag(r2), lam), lam);
    fprintf('off-grid by [%.1f %.1f] bin -> Lemma-1 error %.1e, HPoF NMSE = %7.1f dB\n', offs(i,:), ...
            norm(dht(imag(r2)) - real(r2))/norm(real(r2)), nmse(h2I + 1i*h2Q, r2));
end

% ---- 6. figures ----------------------------------------------------------------------------
n = (0:N-1).';
figure('Name', 'Exp1: signals', 'Position', [50 50 1100 750]);
subplot(3,1,1); plot(n, rI, 'k', n, cI, 'r', n, zI, 'b.-'); yline([-lam lam], ':');
legend('true r_I', 'clipped ADC', 'modulo ADC'); title('I channel'); xlabel('sample n'); ylabel('V');
subplot(3,1,2); plot(n, rQ, 'k', n, cQ, 'r', n, zQ, 'b.-'); yline([-lam lam], ':');
legend('true r_Q', 'clipped ADC', 'modulo ADC'); title('Q channel'); xlabel('sample n'); ylabel('V');
subplot(3,1,3); plot(n, rI, 'k', n, hI, 'go', n, uI, 'm--');
legend('true r_I', 'HPoF', 'unwrapping (Itoh)'); title('Recovered I channel'); xlabel('sample n'); ylabel('V');

figure('Name', 'Exp1: range profiles', 'Position', [80 80 900 450]);
plot(rng_axis, rprof(rC), 'k', 'LineWidth', 2); hold on
plot(rng_axis, rprof(cI + 1i*cQ), 'r');
plot(rng_axis, rprof(zI + 1i*zQ), 'Color', [0.6 0.6 0.6]);
plot(rng_axis, rprof(hI + 1i*hQ), 'g--', 'LineWidth', 1.5);
xline(R, ':'); xlim([0 (N/2-1)*dR]); ylim([-120 20]);
legend('ideal ADC', 'clipped ADC', 'modulo, not unfolded', 'modulo + HPoF');
xlabel('range [m]'); ylabel('|X| [dBV]'); title('Range profile: strong target at 3 m, weak (-40 dB) at 30 m');

figure('Name', 'Exp1: singular values', 'Position', [110 110 750 420]);
semilogy(info.sv / info.sv(1), 'o-');
xline(info.K + 0.5, 'r--', sprintf('K = %d', info.K));
grid on;
title('Singular values of the stacked Hankel matrix');
xlabel('Singular-value index');
ylabel('Normalized singular value');
