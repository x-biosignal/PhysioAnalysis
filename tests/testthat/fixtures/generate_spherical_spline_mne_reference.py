"""Generate the MNE-Python spherical-spline reference as JSON."""

import json
import sys

import mne
import numpy as np
from mne.channels.interpolation import _make_interpolation_matrix


def lift(x, y, center, radius):
    u = (x - center[0]) / radius
    v = (y - center[1]) / radius
    return np.column_stack(
        [u, v, np.sqrt(np.maximum(0.0, 1.0 - u * u - v * v))]
    )


def main(output_path):
    angles = np.arange(8) * (2 * np.pi / 8)
    x = np.r_[0.8 * np.cos(angles), 0.0]
    y = np.r_[0.8 * np.sin(angles), 0.0]
    values = np.array([1.0, 0.5, -0.2, -0.8, -1.0, -0.4, 0.3, 0.9, 0.1])
    query_x = np.array(
        [
            -0.70,
            -0.45,
            -0.20,
            0.05,
            0.30,
            0.55,
            0.72,
            -0.60,
            -0.35,
            -0.10,
            0.15,
            0.40,
            0.62,
            0.00,
            0.25,
        ]
    )
    query_y = np.array(
        [
            -0.35,
            0.10,
            0.55,
            -0.65,
            -0.25,
            0.20,
            0.45,
            0.55,
            -0.55,
            0.25,
            0.60,
            -0.45,
            -0.10,
            0.00,
            0.35,
        ]
    )
    center = np.array([0.0, 0.0])
    radius = 1.0
    matrix = _make_interpolation_matrix(
        lift(x, y, center, radius),
        lift(query_x, query_y, center, radius),
        alpha=None,
    )

    payload = {
        "x": x.tolist(),
        "y": y.tolist(),
        "values": values.tolist(),
        "query_x": query_x.tolist(),
        "query_y": query_y.tolist(),
        "center": center.tolist(),
        "radius": radius,
        "stiffness": 4,
        "n_terms": 50,
        "regularization": 0,
        "expected": (matrix @ values).tolist(),
        "provenance": {
            "reference": (
                "MNE-Python "
                "mne.channels.interpolation._make_interpolation_matrix"
            ),
            "mne_version": mne.__version__,
            "mne_revision": "v1.10.1",
            "source_url": (
                "https://github.com/mne-tools/mne-python/blob/"
                "v1.10.1/mne/channels/interpolation.py"
            ),
            "python_version": sys.version.split()[0],
            "numpy_version": np.__version__,
            "alpha": None,
            "stiffness": 4,
            "n_legendre_terms": 50,
            "generation_date": "2026-07-28",
            "command": (
                "_MNE_FAKE_HOME_DIR=/tmp/mne-home "
                "PYTHONPATH=/tmp/mne-wscb10 python3 "
                "tests/testthat/fixtures/"
                "generate_spherical_spline_mne_reference.py "
                "/tmp/wscb10_mne_fixture.json; "
                "Rscript -e 'x <- jsonlite::read_json("
                "\"/tmp/wscb10_mne_fixture.json\", simplifyVector=TRUE); "
                "saveRDS(x, \"tests/testthat/fixtures/"
                "spherical_spline_mne_reference.rds\", version=3)'"
            ),
        },
    }
    with open(output_path, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} OUTPUT_JSON")
    main(sys.argv[1])
