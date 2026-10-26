# Reproducing the computational release

Run commands from the repository root. Install the package with development
extras as in the README. The unit tests and guided teaching experiment need
no release asset download:

```sh
python -m pytest -q
python -m lossqcd.walkthrough --output walkthrough_repeat
```

## Recheck the seven core saved studies

The raw archives are retained locally and have not been publicly uploaded.
These commands require matching archives obtained from the maintainer; the
unit tests, plots and readable reports are available without them. To obtain
new evidence independently, use the regeneration commands below.

```sh
python -m lossqcd.evidence --scope core --archive-dir /path/to/archives
python -m lossqcd.reproduction_report --output review_repeat
python -m lossqcd.reproduction_report --check review_repeat
```

`lossqcd.evidence` reads the supplied archive, verifies its SHA-256 and
every member's hash/size, then restores absent files. It refuses conflicting
files and unsafe archive entries. Already identical files are preserved.
Use `--archive-dir` to avoid network access:

```sh
python -m lossqcd.evidence --scope core --archive-dir /path/to/downloads
```

The report recalculates formulas and saved-data summaries using six auditors
plus Figure 2 checks. It does not regenerate all historical random ensembles.
The supplied public report is a dated snapshot; generating a new report records
the current checkout hashes. Historical reports/manifests describe their own
earlier source snapshots and should not be treated as current-checkout checks.

## Regenerate random ensembles

These can take substantially longer and require fresh output folders. The
recorded exact environment is Python 3.14.2 plus `requirements-lock.txt`.
Use this lock only on a compatible Python/platform. A compatible dependency
installation is sufficient for ordinary use but may change exact random arrays
or floating-point results.

```sh
python -m lossqcd.author_comparison --output author_repeat
python -m lossqcd.audit_author_comparison author_repeat
python -m lossqcd.paper_completion --source author_repeat --output paper_repeat
python -m lossqcd.audit_paper_completion paper_repeat --regenerate --regenerate-calibration
python -m lossqcd.verified_figures --output figures_repeat
python -m lossqcd.equation_validation --source figures_repeat --output equations_repeat
python -m lossqcd.audit_equation_validation equations_repeat --regenerate
python -m lossqcd.delay_scaling --output scaling_repeat
python -m lossqcd.audit_delay_scaling scaling_repeat --regenerate
python -m lossqcd.dual_rail_experiment --output dual_rail_repeat
python -m lossqcd.audit_dual_rail dual_rail_repeat --regenerate
```

Figure 2 regeneration is included in `lossqcd.entangled_experiments`; its other
defaults run a larger detection study. For formulas alone, the guided lesson
and combined checker recalculate the saved curves without that Monte Carlo cost.

## Additional studies

```sh
python -m lossqcd.evidence --scope all --archive-dir /path/to/archives
python -m lossqcd.audit_adaptive_validation results_adaptive_validation
python -m lossqcd.audit_null_channel_diagnostics results_null_channel_refined
python -m lossqcd.audit_entangled_null_diagnostics results_entangled_null_diagnostics
python -m lossqcd.audit_robust_thresholds results_robust_thresholds
```

Optional `--regenerate` flags are documented by each module's `--help`. The
conservative-threshold auditor also has `--regenerate-calibration`. The four
checks above are separate from the core paper-report total. Saved historical
failed runs are retained and described as such, not silently relabeled passes.

## Public-source audit and provenance

The small CC BY 4.0 Zenodo tables/script are included and verified against their
published checksums. No MATLAB execution is needed. Manuscript/figure source
archives are not redistributed. To rerun `scripts/audit_arxiv_sources.py`, obtain
v1/v2 archives directly from the source URLs in `source_audit/source_inventory.json`
and place them at `data/arxiv_source_audit/2503.12276v1.tar` and `...v2.tar`.
Check the saved archive hashes before rerunning. The audit is a bounded source
availability observation, not proof that no original simulator exists.

Working JSON/notebook paths are relative to the repository root. The publication
migration record preserves old/new file hashes, while untouched source snapshots
remain in the private local project. Personal local paths are not distributed
as provenance copies. Scientific array/table bytes remain unchanged. The release
assets hold raw arrays; `release/evidence.json` identifies them individually.

ARL intervals are approximate Monte Carlo evidence. Detection delay means
E[T−ν | T≥ν], with early alarms separate; it is not Lorden worst-case delay.
See [FINDINGS.md](FINDINGS.md) before interpreting published-figure agreement.
