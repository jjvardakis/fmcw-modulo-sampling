%% EXPERIMENT 2b - Noise floor set by the ADC (b-bit quantization)
%  Same bit budget b for every ADC:
%    conventional ADC, gain set so the strong echo fits (no clipping): step = 2*FS/2^b, FS = peak of r
%    clipping ADC with range +-lam                                   : step = 2*lam/2^b  (fine, but clips)
%    modulo ADC with range +-lam, then HPoF                          : step = 2*lam/2^b  (fine, no clipping)
%  -> the modulo ADC lowers the quantization floor by 20*log10(FS/lam) dB without clipping.
close all; setup_scene; rng(3);
quickRun = exist('v2_quick','var') && v2_quick;
T     = 30;                                  % trials per bit value
bits  = 3:12;                                % 1 x 10 ADC resolutions
if quickRun, T=2; bits=[6 10 12]; end
sig_w = 1e-4;                                % tiny thermal noise [V] (acts as dither; quantization dominates)
FS    = sum(amp);                            % conventional ADC full scale = largest possible |rI|,|rQ| [V]
quant = @(x, F, b) min(max( (2*F/2^b)*(floor(x/(2*F/2^b)) + 0.5), -F + F/2^b), F - F/2^b);  % mid-rise b-bit, N x 1
others = setdiff(0:N-1, kbin) + 1;           % 1 x (N-2) MATLAB indices of bins without a target
out = cell(1, 3);                            % the 3 receiver outputs, each N x 1
nm  = zeros(3, numel(bits), T);              % 3 receivers x bits x trials : NMSE [dB]
wtf = zeros(3, numel(bits), T);              % weak-target-to-floor ratio in the range profile [dB]
for ib = 1:numel(bits)
    b = bits(ib);
    for t = 1:T
        rC = fmcw_beat(R, amp.*exp(1i*2*pi*rand(1,2)), fc, B, Tc, N);
        y  = rC + sig_w/sqrt(2)*(randn(N,1) + 1i*randn(N,1));          % N x 1 analog input
        out{1} = quant(real(y), FS, b) + 1i*quant(imag(y), FS, b);    % conventional, full range
        out{2} = quant(clip(real(y)), lam, b) + 1i*quant(clip(imag(y)), lam, b);   % clipping ADC
        [hI, hQ] = hpof(quant(modlam(real(y), lam), lam, b), quant(modlam(imag(y), lam), lam, b), lam);
        out{3} = hI + 1i*hQ;                                          % modulo ADC + HPoF
        for j = 1:3
            X  = abs(fft(out{j})).^2;                                 % N x 1 power spectrum
            nm(j, ib, t)  = nmse(out{j}, rC);
            wtf(j, ib, t) = 10*log10( X(kbin(2)+1) / mean(X(others)) );
        end
    end
end
fprintf('expected floor gain of the modulo ADC: 20log10(FS/lam) = %.1f dB (= %.2f bits)\n', 20*log10(FS/lam), log2(FS/lam));
col = {'k', 'r', 'g'};  lab = {'conventional ADC (full range)', 'clipping ADC (+-\lambda)', 'modulo ADC + HPoF'};
figure('Name', 'Exp2b: quantization floor', 'Position', [50 50 1100 420]);
subplot(1,2,1); hold on
for j = 1:3, plot(bits, median(squeeze(nm(j,:,:)), 2), [col{j} 'o-']); end
grid on; xlabel('ADC bits b'); ylabel('median NMSE [dB]'); legend(lab, 'Location', 'southwest'); title('Reconstruction error');
subplot(1,2,2); hold on
for j = 1:3, plot(bits, median(squeeze(wtf(j,:,:)), 2), [col{j} 'o-']); end
grid on; xlabel('ADC bits b'); ylabel('weak target / mean floor [dB]'); title('Weak (-40 dB) target above the floor');

% ---- range profiles at one resolution ------------------------------------------------------
b  = 6;
rC = fmcw_beat(R, amp.*exp(1i*2*pi*rand(1,2)), fc, B, Tc, N);
y  = rC + sig_w/sqrt(2)*(randn(N,1) + 1i*randn(N,1));
[hI, hQ] = hpof(quant(modlam(real(y), lam), lam, b), quant(modlam(imag(y), lam), lam, b), lam);
figure('Name', 'Exp2b: range profiles', 'Position', [80 80 900 420]); hold on
plot(rng_axis, rprof(quant(real(y), FS, b) + 1i*quant(imag(y), FS, b)), 'k');
plot(rng_axis, rprof(quant(clip(real(y)), lam, b) + 1i*quant(clip(imag(y)), lam, b)), 'r');
plot(rng_axis, rprof(hI + 1i*hQ), 'g');
xline(R, ':'); xlim([0 (N/2-1)*dR]); ylim([-90 20]); grid on
legend(lab); xlabel('range [m]'); ylabel('|X| [dBV]'); title(sprintf('%d-bit ADCs', b));
