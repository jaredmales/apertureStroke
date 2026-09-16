#!/usr/bin/env python3
"""Report the maximum adjacent-pixel stroke in each residual-phase FITS cube frame.

This reproduces aperture_stroke::maxAbsPixelDiff: only nonzero pupil pixels
are considered, and only the down and right neighbors are tested.  Coordinates
are zero-based C++ image coordinates: row, column (equivalently FITS y, x).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def maximum_pixel_difference(image: np.ndarray, pupil: np.ndarray) -> tuple[float, int, int, int, int]:
    """Return the same maximum and position selected by maxAbsPixelDiff."""
    maximum = 0.0
    maximum_position: tuple[int, int, int, int] | None = None
    rows, columns = image.shape

    # Preserve apertureStroke's column-major scan order and strict comparison,
    # including its deterministic tie breaking.
    for column in range(columns):
        for row in range(rows):
            if pupil[row, column] == 0:
                continue

            if row + 1 < rows and pupil[row + 1, column] != 0:
                difference = abs(float(image[row + 1, column] - image[row, column]))
                if difference > maximum:
                    maximum = difference
                    maximum_position = (row, column, row + 1, column)

            if column + 1 < columns and pupil[row, column + 1] != 0:
                difference = abs(float(image[row, column + 1] - image[row, column]))
                if difference > maximum:
                    maximum = difference
                    maximum_position = (row, column, row, column + 1)

    if maximum_position is None:
        raise ValueError("no nonzero adjacent pupil-pixel pair was found")

    return maximum, *maximum_position


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cube", type=Path, help="Residual-phase FITS cube to inspect.")
    parser.add_argument("pupil", type=Path, help="Binary pupil FITS image used by apertureStroke.")
    args = parser.parse_args()

    try:
        cube = np.asarray(fits.getdata(args.cube))
        pupil = np.asarray(fits.getdata(args.pupil))
    except OSError as error:
        print(error, file=sys.stderr)
        return 2
    if cube.ndim == 2:
        cube = cube[np.newaxis, :, :]
    if cube.ndim != 3:
        print(f"expected a 2-D image or 3-D cube, got {cube.ndim} dimensions", file=sys.stderr)
        return 2
    if pupil.ndim != 2 or pupil.shape != cube.shape[1:]:
        print(f"pupil shape {pupil.shape} does not match cube frame shape {cube.shape[1:]}", file=sys.stderr)
        return 2

    print("# frame row col neighbor_row neighbor_col stroke_rad")
    cube_maximum: tuple[float, int, int, int, int, int] | None = None
    for frame, image in enumerate(cube):
        try:
            result = maximum_pixel_difference(image, pupil)
        except ValueError as error:
            print(f"{args.cube}: {error}", file=sys.stderr)
            return 2

        stroke, row, column, neighbor_row, neighbor_column = result
        print(f"{frame} {row} {column} {neighbor_row} {neighbor_column} {stroke:.9g}")
        candidate = (stroke, frame, row, column, neighbor_row, neighbor_column)
        if cube_maximum is None or candidate[0] > cube_maximum[0]:
            cube_maximum = candidate

    assert cube_maximum is not None
    stroke, frame, row, column, neighbor_row, neighbor_column = cube_maximum
    print(f"# cube_max frame={frame} row={row} col={column} "
          f"neighbor_row={neighbor_row} neighbor_col={neighbor_column} stroke_rad={stroke:.9g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
