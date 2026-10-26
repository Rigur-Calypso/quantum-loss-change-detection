# Publication validation

Validated locally on 26 October 2026 using Python 3.14.2 and the recorded dependencies.

- **421 unit tests pass** in the publication copy.
- **1,869 core saved-evidence checks pass**; **701 file hashes verified**.
- **984 additional checks pass** across adaptive ARL, scalar/correlated initial-channel mismatch, and conservative thresholds.
- Both release archives and all **739 raw-array members** were restored/verified offline with the public loader. Existing identical files were retained.
- Original local project reads were blocked during tests/core/extension audits, except the installed dependency environment. The publication source was imported explicitly.

This validates a relocated source/data copy in the existing numerical environment.
It is not a fresh optical experiment or full historical Monte Carlo regeneration.
GitHub Actions additionally runs the unit tests from a fresh checkout; see the repository Actions page for its live status.

[Machine-readable record](validation.json) · [Unit-test log](validation_tests.txt) · [Core log](validation_core.txt) · [Extension summary](extension_audits/summary.json)

The raw archives remain local; no public raw-data release is available yet.
The clean committed checkout also passed all 421 tests without raw arrays; see [clean-checkout record](clean_checkout.json). GitHub Actions passed on Python 3.11 and 3.14 for the initial published commit.

The current [consolidated report](../paper_reproduction_public/REPORT.md) recalculates saved arrays. Earlier report folders record earlier snapshots and may have stale current-checkout hashes or historical limitations. Always use the current README and FINDINGS.md for the release's scientific claims.

**Successful checks do not establish full published-figure agreement.** See [remaining gaps](../FINDINGS.md).
