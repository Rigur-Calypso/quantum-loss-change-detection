import numpy as np
import pytest
from lossqcd.octave_data import parse_octave, load_octave


def test_author_nested_cells_and_range():
    data = load_octave("data/zenodo_17460928/prl_fig7_data")
    np.testing.assert_allclose(data["Na_range"], np.arange(0.1, 0.501, 0.05), atol=1e-15)
    np.testing.assert_array_equal(data["blocksizes"], [[1, 2, 4, 8]])
    assert len(data["mean_sqentlat"]) == 4
    assert data["mean_sqentlat"][3][0, 0] == pytest.approx(19.9756)
    assert data["mean_clat"][0, 0] == pytest.approx(64.624)
    assert data["styles"] == "+*xs"


def test_range_does_not_force_last_value_to_upper_limit():
    x = parse_octave(
        "# name: r\n# type: double_range\n# base, limit, increment\n0.1 0.881373587019543 0.1\n"
    )["r"]
    np.testing.assert_allclose(x, np.arange(1, 9) / 10)
    assert len(x) == 8 and x[-1] < 0.881373587019543


def test_matrix_shape_and_negative_range_increment():
    x = parse_octave("# name: x\n# type: matrix\n# rows: 2\n# columns: 2\n1 2\n3 4\n")
    np.testing.assert_array_equal(x["x"], [[1, 2], [3, 4]])
    x = parse_octave("# name: r\n# type: double_range\n# base, limit, increment\n0.5 0.1 -0.1\n")[
        "r"
    ]
    np.testing.assert_allclose(x, [0.5, 0.4, 0.3, 0.2, 0.1])


@pytest.mark.parametrize(
    "text",
    [
        "# name: x\n# type: object\n",
        '# name: x\n# type: scalar\nsystem("bad")\n',
        "# name: x\n# type: scalar\nnan\n",
        "# name: x\n# type: matrix\n# rows: 2\n# columns: 2\n1 2\n",
        "# name: x\n# type: matrix\n# rows: -1\n# columns: 2\n",
        "# name: x\n# type: scalar\n1 2\n",
        "# name: x\n# type: scalar\n1\n# name: x\n# type: scalar\n2\n",
        "# name: r\n# type: double_range\n# base, limit, increment\n0 1 0\n",
        "# name: c\n# type: cell\n# rows: -1\n# columns: -1\n",
        "# name: s\n# type: sq_string\n# elements: 1\n# length: 4\nx\n",
    ],
)
def test_invalid_or_executable_content_refused(text):
    with pytest.raises(ValueError):
        parse_octave(text)
