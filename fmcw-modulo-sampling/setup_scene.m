%% SETUP_SCENE  Shared FMCW scene + helpers for all experiments.
%  This is a SCRIPT (not a function): every variable below lands in the workspace
%  of the script that calls it. Edit the numbers here and every experiment follows.

c   = 3e8;            % speed of light [m/s]
fc  = 77e9;           % chirp start frequency [Hz]
B   = 150e6;          % sweep bandwidth [Hz]      -> range bin dR = c/(2B) = 1 m
N   = 256;            % ADC samples per chirp      -> usable bins 1..N/2-1, R_max = 127 m
Tc  = 25.6e-6;        % chirp duration [s]         -> ADC rate fs = N/Tc = 10 MHz (per channel)
dR  = c/(2*B);        % range resolution [m]

lam = 2.4;            % ADC input range [V]: clipping ADC clips at +-lam, modulo ADC folds into [-lam, lam)

% Two-target scene: a strong near echo and a weak far echo (equal RCS).
R          = [3 30];                   % 1 x 2 ranges [m] (strong, weak); multiples of dR = on the range grid
amp_strong = 2.5*lam;                  % strong-echo amplitude [V] (2.5x above the ADC range)
amp        = amp_strong*(R(1)./R).^2;  % 1 x 2 amplitudes: radar eq. -> amplitude ~ 1/R^2  (weak = -40 dB)
kbin       = round(R/dR);              % 1 x 2 FFT bins of the targets (0-based) = [3 30]

% ---- small helper one-liners (anonymous functions, no classes) ----
clip  = @(x) min(max(x, -lam), lam);                       % conventional ADC (clips at +-lam)
nmse  = @(rh, r) 10*log10( 0.5*norm(real(rh-r))^2/norm(real(r))^2 ...
                         + 0.5*norm(imag(rh-r))^2/norm(imag(r))^2 );   % paper's NMSE [dB]
rprof = @(x) 20*log10(abs(fft(x))/N + 1e-9);               % N x 1 range profile [dBV], fft(N x 1) -> N x 1
% Itoh unwrapping: wrap every sample-to-sample jump back into [-lam, lam), then cumsum (N x 1).
% Valid when |r(n+1) - r(n)| < lam. The constant is fixed like in hpof.m (r has zero mean).
itoh0 = @(z) z(1) + [0; cumsum(diff(z) - 2*lam*round(diff(z)/(2*lam)))];
itoh  = @(z) itoh0(z) + 2*lam*round(-mean(itoh0(z))/(2*lam));
rng_axis = (0:N-1).'*dR;                                   % N x 1 range axis [m] for bin index 0..N-1
