%% EXPERIMENT 3 - NMSE, probability of detection and false alarms vs SNR (thermal noise)
%  Noise is added BEFORE the ADC (front-end noise):  y = r + w,  w ~ CN(0, sigma^2) per sample.
%  SNR is defined on the STRONG target per sample: SNR = |a_strong|^2/sigma^2 ;  weak target = SNR - 40 dB.
%  Detection: CA-CFAR on |FFT|^2, Pd = P(weak-target bin declared), Pfa_meas = declared bins / empty bins.
%  Two strong-echo levels show how the noise threshold of HPoF depends on the number of folds.
close all; setup_scene; rng(4);
quickRun = exist('v2_quick','var') && v2_quick;
T      = 200;                                % trials per SNR point
snr_dB = 10:5:70;                            % 1 x 13 per-sample SNR of the strong target
levels = [2.5 4]*lam;                        % 1 x 2 strong-echo amplitudes (fewer / more folds)
if quickRun, T=2; snr_dB=[20 40 60]; end
Pfa = 1e-3;  Nt = 16;  Ng = 2;               % CFAR settings
others = setdiff(0:N-1, kbin) + 1;           % 1 x (N-2) indices of empty bins
names  = {'ideal ADC (infinite range)', 'clipping ADC', 'modulo ADC + HPoF'};

NM = zeros(3, numel(snr_dB), T, numel(levels));   % NMSE [dB]       3 receivers x SNR x trials x levels
PD = zeros(3, numel(snr_dB), T, numel(levels));   % weak target detected (0/1)
FA = zeros(3, numel(snr_dB), T, numel(levels));   % fraction of empty bins declared
OK = zeros(numel(snr_dB), T, numel(levels));      % HPoF found every fold exactly (0/1)
for il = 1:numel(levels)
    a0 = levels(il)*(R(1)./R).^2;            % 1 x 2 amplitudes for this level
    for is = 1:numel(snr_dB)
        sigma = a0(1)/10^(snr_dB(is)/20);    % noise std (complex) [V]
        for t = 1:T
            rC = fmcw_beat(R, a0.*exp(1i*2*pi*rand(1,2)), fc, B, Tc, N);     % N x 1 noiseless
            y  = rC + sigma/sqrt(2)*(randn(N,1) + 1i*randn(N,1));          % N x 1 noisy analog input
            [hI, hQ] = hpof(modlam(real(y), lam), modlam(imag(y), lam), lam);
            out = {y, clip(real(y)) + 1i*clip(imag(y)), hI + 1i*hQ};        % 1 x 3 cell of N x 1
            OK(is, t, il) = max(abs(out{3} - y)) < 1e-6;
            for j = 1:3
                hit = cfar_ca(abs(fft(out{j})).^2, Nt, Ng, Pfa);            % N x 1 logical
                NM(j, is, t, il) = nmse(out{j}, rC);                        % error w.r.t. the NOISELESS signal
                PD(j, is, t, il) = hit(kbin(2) + 1);
                FA(j, is, t, il) = mean(hit(others));
            end
        end
        fprintf('level %.1f lam, SNR %2d dB: HPoF exact-unfolding rate %.2f\n', levels(il)/lam, snr_dB(is), mean(OK(is,:,il)));
    end
end

% ---- figures (solid = 2.5*lam scene, dashed = 4*lam scene) ---------------------------------
col = {'k', 'r', 'g'};  sty = {'-', '--'};
figure('Name', 'Exp3: SNR sweep', 'Position', [40 40 1200 800]);
for il = 1:numel(levels)
    for j = 1:3
        subplot(2,2,1); hold on
        plot(snr_dB, median(NM(j,:,:,il), 3), [col{j} sty{il} 'o'], 'DisplayName', sprintf('%s, %.1f\\lambda', names{j}, levels(il)/lam));
        subplot(2,2,2); hold on
        plot(snr_dB, mean(PD(j,:,:,il), 3), [col{j} sty{il} 'o'], 'DisplayName', sprintf('%s, %.1f\\lambda', names{j}, levels(il)/lam));
        subplot(2,2,3); hold on
        semilogy(snr_dB, mean(FA(j,:,:,il), 3) + 1e-6, [col{j} sty{il} 'o'], 'DisplayName', sprintf('%s, %.1f\\lambda', names{j}, levels(il)/lam));
    end
    subplot(2,2,4); hold on
    plot(snr_dB, mean(OK(:,:,il), 2), ['g' sty{il} 'o'], 'DisplayName', sprintf('strong echo %.1f\\lambda', levels(il)/lam));
end
subplot(2,2,1); grid on; xlabel('SNR of strong target [dB]'); ylabel('median NMSE [dB]'); title('Reconstruction error vs noiseless r'); legend('Location', 'southwest');
subplot(2,2,2); grid on; xlabel('SNR of strong target [dB]  (weak target = SNR - 40 dB)'); ylabel('P_D (weak target)'); title(sprintf('Detection, CA-CFAR, P_{fa} = %g', Pfa));
subplot(2,2,3); grid on; set(gca, 'YScale', 'log'); yline(Pfa, ':', 'design P_{fa}');
xlabel('SNR of strong target [dB]'); ylabel('measured false-alarm rate per bin'); title('Ghost targets');
subplot(2,2,4); grid on; xlabel('SNR of strong target [dB]'); ylabel('fraction of trials'); title('HPoF: all folds recovered exactly'); legend('Location', 'southeast');
