from pathlib import Path
import shutil
import numpy as np
import pytest
from lossqcd.authors import author_data, main


def test_author_metadata_and_analytical_ratios_are_recovered_exactly():
    data = author_data("data/authors")
    assert data["figure3"]["max_absolute_kl_error"] < 1e-12
    assert data["figure7"]["cp"] == 1005 and data["figure7"]["saved_n"] == 5000
    assert len(data["figure7"]["energies"]) == 9
    assert data["supplement_raw_latency_negative_count"] == 2
    assert data["figure3"]["energies"][-1] == pytest.approx(np.sinh(0.8) ** 2)


def test_corrupt_source_data_refused(tmp_path):
    source = Path("data/authors")
    for path in source.iterdir():
        shutil.copyfile(path, tmp_path / path.name)
    with (tmp_path / "prl_fig7_data").open("a") as f:
        f.write("\n")
    with pytest.raises(ValueError, match="checksum failed"):
        author_data(tmp_path)


@pytest.mark.parametrize(
    "args",
    [
        ["--delay-trials", "1"],
        ["--batch-size", "1"],
        ["--calibration-cycles", "0"],
        ["--validation-cycles", "123"],
    ],
)
def test_invalid_budgets_do_not_create_output(args, tmp_path):
    output = tmp_path / "results"
    with pytest.raises(SystemExit):
        main(["--output", str(output)] + args)
    assert not output.exists()


def test_small_end_to_end_comparison_and_saved_data_audit(tmp_path):
    import json
    from lossqcd.audit_authors import audit

    output = tmp_path / "comparison"
    main(
        [
            "--output",
            str(output),
            "--calibration-cycles",
            "1000",
            "--validation-cycles",
            "1000",
            "--batch-size",
            "1000",
            "--delay-trials",
            "100",
        ]
    )
    report = json.loads((output / "report.json").read_text())
    assert len(report["models"]) == 61
    assert len(report["figure7_rows"]) == 45 and len(report["figure3_rows"]) == 8
    assert (output / "figure7_overlay.png").exists() and (output / "figure3_overlay.png").exists()
    assert audit(output)["data_audit_passed"]


def test_archived_author_path_resolves_to_renamed_data(tmp_path):
    source = Path("data/authors")
    short = tmp_path / "authors"
    shutil.copytree(source, short)
    assert author_data(tmp_path / "zenodo_17460928") == author_data(short)


def test_archived_author_file_maps_to_short_path(tmp_path):
    from lossqcd.authors import author_path

    short = tmp_path / "authors"
    short.mkdir()
    path = short / "manifest.json"
    path.write_text("verified source")
    assert author_path(tmp_path / "zenodo_17460928/manifest.json") == path


def test_existing_legacy_author_file_is_not_substituted(tmp_path):
    from lossqcd.authors import author_path

    old = tmp_path / "zenodo_17460928"
    old.mkdir()
    path = old / "manifest.json"
    path.write_text("unverified source")
    assert author_path(path) == path
