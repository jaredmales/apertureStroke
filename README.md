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
matrices or generating turbulence, set `--trials=0`.

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
