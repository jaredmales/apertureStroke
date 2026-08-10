# apertureStroke

`apertureStroke` generates atmospheric phase screens, removes configured
mode cutoffs, and measures residual P2V, nearest-neighbor pixel difference,
mode amplitudes, and modified-Fourier amplitudes. All reported amplitudes are
microns of optical surface (`wavefront / 2`).

The pupil is either a constructed circular aperture, optionally with a central
obscuration, or a binary mask read from a FITS file. A primary basis can be
generated as true Zernikes or read from a FITS cube. An optional FITS prefix
basis is included at every primary-mode cutoff; this is the generic mechanism
for experiments such as GMT segment piston/tip/tilt plus N Zernikes.

Build and inspect the configuration options:

```bash
make
./apertureStroke --help
```

Run the supplied circular and GMT configurations:

```bash
./apertureStroke -c circle_25m.conf
./apertureStroke -c gmt_25m.conf
./apertureStroke -c magaox_6p5m.conf
```

Command-line options override values from the configuration file. For a quick
setup check that writes the pupil and optional basis without constructing fit
matrices or generating turbulence, set `--simulation.trials=0`.

Long options follow the same `section.key` convention as configuration-file
targets. For example:

```bash
--atmosphere.subharmonicLevel=4 --output.directory=output/circle_25m_sh4
```

Turbulence is generated with a single unit-strength layer; the LCO profile is
used only to supply default seeing and outer scale values. The piston and
tip/tilt PSD-removal transfer functions are controlled independently. To
generate the unfiltered Kolmogorov PSD used by the Noll prediction, set both
to false:

```bash
--atmosphere.psdSubtractPiston=false \
--atmosphere.psdSubtractTipTilt=false
```

Output filenames identify these settings with the tags `1layer`, `psdP`, and
`psdTT`.

To test independent modal projection instead of a simultaneous least-squares
fit, use `basis.fit=projection`. Every coefficient is measured against the
same piston-removed input screen, and the fixed coefficients are then
subtracted cumulatively for the requested cutoffs:

```bash
./apertureStroke -c circle_25m.conf \
    --atmosphere.subharmonicLevel=4 \
    --basis.fit=projection \
    --output.directory=output/circle_25m_sh4_projection
```

The exact configured `output.directory` is created recursively. Pupil, phase,
run-parameter, histogram, and statistics products are written beneath it.

To prepend a saved GMT segment piston/tip/tilt cube to every Zernike cutoff,
add the following to a GMT configuration:

```ini
[basis]
type=zernike
modes=1000
prefixFile=/path/to/gmt_segment_modes_182x182.fits
prefixName=segment
```

Here, cutoff 0 fits only the segment modes and cutoff 100 fits the segment
modes plus 100 Zernikes. Prefix and primary amplitude statistics are written
to separate tables.

For a 50-pixel MagAO-X grid with a 48-pixel, 6.5 m pupil, use:

```ini
[pupil]
arraySize=50
diameterPixels=48
diameterMeters=6.5
centralObscuration=0.27
```

The turbulence grid diameter is calculated as
`diameterMeters * arraySize / diameterPixels`, preserving the specified
pixels per physical pupil even when the pupil is padded in its array.

## Circular Zernike theory

`zernikeTheory` calculates the Noll Kolmogorov prediction for an unobstructed
circular aperture. It constructs and normalizes the same discrete pupil and
Zernike basis as `apertureStroke`, then converts modal P2V to microns surface.

The simulation configuration can be used directly:

```bash
make zernikeTheory
./zernikeTheory -c circle_25m.conf
```

This writes `output/circle_25m/zernike_theory_circle_25m.dat`. Columns 2 and 3
are the predicted mean absolute modal P2V and RMS about that mean, matching the
`mean` and `rms` columns in the measured `zernike_p2v` tables. Seeing is
converted to r0 at its reference wavelength and then scaled to the configured
phase wavelength before applying the Noll variance.

The utility intentionally rejects FITS pupils, central obscurations, nonzero
outer scale, and non-Zernike bases because those are outside the circular Noll
prediction.
