# Method and metrics

## Signal model

For chirp slope S = B/Tc and round-trip delay tau = 2R/c, the dechirped return has phase `2*pi*S*tau*t + 2*pi*fc*tau - pi*S*tau^2`. Its beat frequency is S*tau. The implementation uses transmitted signal times conjugated received signal, giving positive beat frequencies for stationary targets.

Target amplitudes can be complex; their arguments are scattering phases. Random-phase experiments hold the magnitudes and ranges fixed and draw independent uniform phases. Different sweep points draw new phases using the same reproducible random stream.

## Reconstruction

The code defines H through the DFT multiplier `1i*sign(frequency)`, assigning +1i to the Nyquist bin. For the ideal on-grid input, H(rQ) = rI. This convention is the negative of another common Hilbert-transform convention. The Nyquist multiplier can make H applied to a general real sequence complex.

Write each modulo channel as z = r - epsilon, where epsilon is an integer multiple of 2*lambda. The signal cancels in `zI - H(zQ)`. Circular differencing exposes residue jumps; the subsequent Fourier-domain model is a sum of exponentials encoding their sample locations. The time-domain cancellation difference need not be sparse, because H spreads Q-channel impulses.

For N=256 the pencil parameter is L=84 and the stacked matrix has 88 rows and 85 columns. The implemented largest-gap rule selects orders 1–84. The algorithm rounds poles to sample indices, merges duplicate indices, fits real jump heights, rounds them to multiples of 2*lambda, and integrates. Zero mean resolves the additive constant for the assumed input class.

Under I/Q mismatch, xI=rI and xQ=(1+g)*(rQ*cos(phi)+rI*sin(phi)). The gain in hpof's argument is the increment g; the linear gain gamma is 1+g. Known-mismatch HPoF modifies both cancellation and jump fitting. Its output remains the imbalanced input channels.

## Success and error measures

| Experiment | Reference and measure |
|---|---|
| 1 | Waveform error relative to clean input; true union of I/Q residue-jump locations; target and distortion peaks |
| 2a | Success if the equal-channel NMSE is below -100 dB, relative to the imbalanced input |
| 2b | Median NMSE relative to the clean input; weak-target-bin power divided by mean power over all non-target bins |
| 3 | Exact unfolding if maximum complex-sample error relative to the noisy input is below 1e-6 V; additional NMSE is relative to the clean signal |
| 4 | Reconstruction success if maximum complex-sample error relative to the noiseless input is below 1e-6 V |

The equal-channel NMSE is `10*log10(0.5*norm(errorI)^2/norm(referenceI)^2 + 0.5*norm(errorQ)^2/norm(referenceQ)^2)`. It weights relative I and Q errors equally. It assumes nonzero reference-channel energies.

Reconstruction success does not require the initially selected order to equal the true count. Extra components can disappear when fitted heights round to zero. Conversely, a correct order alone does not guarantee correct locations or heights.

In Experiment 2b, the full-range ADC half-range is sum(amp)=6.06 V; modulo and clipping use 2.4 V. Quantization steps are 2*range/2^b. The uniform-error approximation predicts about 8.05 dB lower quantization-error power for the smaller range, conditional on successful unfolding and quantization-dominated error. The background metric includes distortion and is not a detection probability.

Experiment 3 has no quantizer. Complex Gaussian noise is added before the ADC model, with variance sigma^2 per complex sample. SNR is defined relative to the strong target. The CFAR detection curves are supplementary outputs; the report's thermal figure uses the exact-unfolding fraction.

## Regime map interpretation

The peak normalization is the maximum absolute sample across both channels, not a separate normalization of each channel. The actual number of targets is min(requested targets, available bins).

The colour in each cell blends the two marginal success rates. It is not the measured probability of a joint outcome. The cell labels give the measured rates. The plotted lines separate cells on opposite sides of mean normalized step size 1 and mean fold count N/3; they are not per-trial guarantees.

The sufficient adjacent-sample condition for Itoh recovery and the ideal HPoF fold-count condition describe different restrictions. Neither method is asserted to dominate all scenes. Finite-precision order estimation, conditioning and lattice rounding also affect the empirical results.
