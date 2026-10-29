function y = dht(x)
%DHT  Discrete Hilbert transform of Lemma 1:  H(x) = IDFT( 1i*sign(w) .* DFT(x) ).
%   x : N x 1 (N even).  y : N x 1 (complex in general: the Nyquist bin gets sign(pi) = +1).
%   For a one-sided (positive-frequency) signal rC = rI + 1i*rQ:  H(rQ) = rI  and  H(rI) = -rQ.
x = x(:);                          % N x 1
N = numel(x);
k = (0:N-1).';                     % N x 1 FFT bin index
s = sign(k - N*(k > N/2));         % N x 1 sign(w): 0 at DC, +1 for bins 1..N/2 (incl. Nyquist), -1 for N/2+1..N-1
y = ifft(1i*s.*fft(x));            % N x 1
end
