# Provenance and packaging changes

The project working copy originated from a supplied archive named `hpof_matlab.zip` and was subsequently adapted for the student experiments. The method is described by Feuillen and Bhandari (SSP 2025), cited in the README. The package does not assert that every supplied script was authored by those paper authors.

This upload bundle includes only the five experiments discussed in the report and their dependencies. It excludes the old working copies, temporary diagnostics, chat files and unrelated experiments.

Packaging changes:

- Restored the 3 m/30 m common reference and made Experiment 2a's 8 m/30 m scene explicit.
- Retained the evaluated numerical core functions unchanged; their hashes are recorded in `core_provenance.md`.
- Gave the final three-level regime experiment the consistent filename `exp4_regime_map.m`.
- Removed workspace-clearing commands that interfered with the runner and added reduced quick-run settings.
- Replaced the obsolete result exporter with one runner that exports each experiment's actual figures and numerical arrays.
- Removed the misleading raw-Q-spike overlay from Experiment 1's internal figure; the singular-value plot remains.
- Corrected comments that described sufficient conditions as necessary or guaranteed empirical success.

No new software licence has been assigned to the supplied implementation. The repository's academic citation and provenance notes do not themselves grant a third-party licence. No third-party paper PDF is included.
