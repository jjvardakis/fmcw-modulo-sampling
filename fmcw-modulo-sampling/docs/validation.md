# Validation of the upload bundle

Checked on 29 October 2026 using MATLAB R2025b.

## Executed checks

```matlab
verify_recovery
s = run_experiments('quick');
assert(all(s.completed));
```

All five experiment scripts completed, including MAT, PDF, PNG and FIG export. The quick mode deliberately uses fewer trials and sweep points; its success percentages are not the report results.

The baseline checks verified the Hilbert relation, residue cancellation, circular differencing, noiseless recovery, fold support, Itoh reconstruction and recovery with known I/Q mismatch. The noiseless HPoF maximum waveform error was 9.93e-16 V; selected order and true fold count were both 24, with exact support recovery. The Itoh maximum error was 5.02e-15 V.

Machine-readable records: [baseline checks](validation/baseline_checks.csv) and [quick execution summary](validation/quick_run_summary.csv).

## Scope of verification

The full Monte Carlo sweeps were not rerun during packaging. Passing the quick checks confirms execution and basic consistency; it does not independently reproduce every numerical claim in the submitted report. Run `run_experiments('full')` to generate the full results from this bundle.

The five core MATLAB functions match the working project's files byte for byte. Generated quick results are retained locally but excluded from the upload archive to avoid confusion with the full report statistics.
