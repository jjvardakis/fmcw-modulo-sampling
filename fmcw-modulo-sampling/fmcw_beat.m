function rC = fmcw_beat(R, a, fc, B, Tc, N)
%FMCW_BEAT  Complex (IQ) dechirped beat signal of one FMCW chirp, TX .* conj(RX).
%   R  : 1 x T target ranges [m]
%   a  : 1 x T complex amplitudes [V] (magnitude = echo strength, angle = scattering phase)
%   rC : N x 1 complex samples, rC = rI + 1i*rQ
%
%   Per target: phase(t) = 2*pi*S*tau*t + 2*pi*fc*tau - pi*S*tau^2,   tau = 2R/c,  S = B/Tc
%   -> a tone at beat frequency S*tau = (R/dR)/Tc, i.e. FFT bin k = R/dR  (always POSITIVE)
c   = 3e8;
S   = B/Tc;                                    % chirp slope [Hz/s]
t   = (0:N-1).'*(Tc/N);                        % N x 1 sample times; (0:N-1).' is a column
tau = 2*R(:).'/c;                              % 1 x T round-trip delays; (:).' forces a row
phs = 2*pi*S*t*tau + 2*pi*fc*tau - pi*S*tau.^2;% N x T: (N x 1)*(1 x T) outer product, + 1 x T row broadcast
rC  = exp(1i*phs)*a(:);                        % (N x T)*(T x 1) -> N x 1 sum over targets
end
