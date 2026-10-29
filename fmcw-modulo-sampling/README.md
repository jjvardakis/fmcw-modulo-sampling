# Modulo sampling via spectrum estimation in FMCW radar

MATLAB code for a Statistical Digital Signal Processing student project comparing conventional ADC acquisition with quadrature modulo acquisition and Hilbert–Pencil-of-Function (HPoF) reconstruction.

The experiments examine clipping, I/Q imbalance, ADC resolution, thermal noise and comparison with Itoh modulo unwrapping. The code uses a single FMCW chirp with stationary targets. The principal experiments use positive, integer-bin beat frequencies so that the periodic discrete Hilbert relation holds.

## Run the project

Tested with MATLAB R2025b. The code uses MATLAB functions and does not require a third-party MATLAB package. Open this folder as the MATLAB current folder, then run:

```matlab
verify_recovery                    % Algebra and noiseless reconstruction checks
run_experiments('quick')            % Small execution checks of all five experiments
run_experiments('full')             % All report experiments at their full trial counts
```

To run and export just one experiment:

```matlab
run_experiments('full', 'exp2a_iq_imbalance')
```

To explore an experiment interactively, run its script directly. Figures are displayed; the runner additionally saves them. If a quick-mode flag was set manually, clear it first:

```matlab
clear v2_quick
exp2a_iq_imbalance
```

The full regime experiment runs 12,000 scenes. Quick runs are execution checks and must not be used as report statistics. Each script resets its stated random seed. Results can depend on numerical linear algebra near singular-value ties; the largest-gap estimator is deliberately retained from the evaluated implementation.

## Report settings

| Experiment | Scene and sweep | Full trials | Seed |
|---|---|---:|---:|
| 1 — clipping versus modulo | 3 m and 30 m; amplitudes 6 V and 0.06 V | One reference scene, plus illustrative off-grid perturbations | 1 |
| 2a — I/Q imbalance | 8 m and 30 m; amplitudes 6 V and 0.426667 V; gain 0:0.25:3 dB and phase 0:5:30 degrees, separately | 50 per point | 2 |
| 2b — ADC resolution | Reference 3 m/30 m scene; 3–12 bits; complex noise standard deviation 1e-4 V | 30 per point | 3 |
| 3 — thermal noise | Reference ranges; strong amplitude 2.5 or 4 times lambda; strong-target SNR 10:5:70 dB | 200 per SNR and amplitude | 4 |
| 4 — regime comparison | Random on-grid tones; requested target counts [1 2 4 8 16]; maximum bins [5 10 20 30 45 60 90 120]; peaks [1.25 1.5 2] times lambda | 100 per cell; 12,000 total | 5 |

Shared settings are N = 256, chirp duration 25.6 microseconds, sampling frequency 10 MHz, carrier frequency 77 GHz, chirp bandwidth 150 MHz and modulo threshold lambda = 2.4 V. Range-bin spacing is 1 m.

`setup_scene.m` defines the 3 m/30 m reference scene. Experiment 2a explicitly replaces the ranges and recomputes amplitudes and bins after loading it. Experiment 4 constructs independent random scenes. No manual scene changes are needed between report experiments.

## Reading order

1. `setup_scene.m` — scene parameters, error metric and Itoh baseline.
2. `fmcw_beat.m` — delayed chirps become beat-frequency tones.
3. `modlam.m` — centered modulo acquisition in each real channel.
4. `dht.m` — periodic discrete Hilbert operator and its sign convention.
5. `hpof.m` — cancellation, circular differences, FFT, stacked Hankel matrix, model-order selection, matrix pencil, jump-height fitting and integration.
6. `exp1_clip_vs_modulo.m` through `exp4_regime_map.m` — experiments.

See [method and metrics](docs/method_and_metrics.md) for notation and interpretation, and [provenance](docs/provenance.md) for the relationship to the supplied implementation.

## Saved outputs

`run_experiments` writes to `results/quick/` or `results/full/`, with a separate folder per experiment:

- `results.mat`: numerical arrays and intermediate structures, retaining their MATLAB dimensions.
- `console.txt`: numerical summaries.
- `.pdf`, `.png`, `.fig`: vector figures, raster previews and editable MATLAB figures.
- `run_summary.csv`: execution status and elapsed time for all selected experiments.

`completed=true` means the script executed successfully, not that HPoF reconstructed every trial. An `error.txt` records a failed script. See [validation](docs/validation.md) for what was executed when packaging this version.

The figures use the experiment layouts. They contain the data behind the report, but their typography or panel arrangement may differ from the final LaTeX layout.

## Scope and limitations

- Model-order selection uses the largest consecutive singular-value ratio. It is a heuristic and does not explicitly select zero order.
- Selected order, distinct rounded locations and final nonzero folds can differ.
- Known-mismatch recovery receives the true gain and phase; it does not estimate calibration.
- Successful unfolding reconstructs the analog input, including its noise and imbalance. It does not denoise or correct I/Q mismatch afterward.
- The ideal fold-count condition is a sufficient theoretical reference, not a guarantee for every finite-precision run of this implementation.
- Independent channel noise and finite-bit quantization disturb the ideal Hilbert relation. Their effect is evaluated empirically.
- The regime map uses Gaussian tone amplitudes scaled to a fixed global channel peak, rather than the inverse-square range law of the reference scene.

## Reference and citation

T. Feuillen and A. Bhandari, “Quadrature Radars Naturally Benefit from Unlimited Sensing,” IEEE Statistical Signal Processing Workshop, 2025, pp. 251–255. DOI: [10.1109/SSP64130.2025.11073462](https://doi.org/10.1109/SSP64130.2025.11073462).

This repository is a student simulation project based on the supplied MATLAB implementation. It is not an official release by the paper's authors. See [CITATION.md](CITATION.md) for how to reference the uploaded repository in the report.
