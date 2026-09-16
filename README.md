# apertureStroke

`apertureStroke` generates atmospheric phase screens, removes configured
mode cutoffs, and measures residual P2V, nearest-neighbor pixel difference,
mode amplitudes, and optional modified-Fourier amplitudes. All reported amplitudes are
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

Trials are independent and run in parallel. By default,
`simulation.threads=0` uses OpenMP's worker limit; set an explicit value to
bound the per-worker turbulence memory, for example:

```bash
--simulation.threads=8
```

Each worker owns a separate turbulence generator, while the basis setup is
shared.

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

`output.writeBasis=true` writes the normalized input modes, while
`output.writeOrthogonalBasis=true` additionally writes
`basis_orthogonalized_<label>.fits`. The latter is a pupil-masked,
pupil-orthonormal version of the full combined basis, made with stabilized
Gram--Schmidt in basis order (prefix modes first). It is a diagnostic only and
does not change the basis used by any fit. For a cutoff, use the first
`prefix-mode count + cutoff` planes.

Primary modes can be smoothed before they are combined with any prefix basis:

```ini
[basis]
smoothingFwhm=5
smoothingStart=100
```

`smoothingStart` is a zero-based primary-mode index, so this leaves the first
100 primary modes unchanged. Each later mode is nearest-neighbor extrapolated
past the binary pupil by the Gaussian support, filtered, and cropped back to
the pupil before the usual RMS normalization. Prefix modes (for example GMT
segment PTT) are not smoothed. Set `smoothingFwhm=0` to disable preprocessing.

For a sampling-defined alternative, use a circular hard Fourier cutoff in
cycles per pixel instead of Gaussian smoothing:

```ini
[basis]
smoothingFwhm=0
lowPassCutoff=0.2
smoothingStart=100
```

The cutoff is radial (so `0.2` cycles/pixel corresponds to a five-pixel
period) and retains spatial frequencies at or below the cutoff. Before the
FFT, each mode is nearest-neighbor extrapolated by one pupil-array width in
every direction; it is then filtered, cropped, pupil-masked, and RMS
normalized. This moves the periodic FFT boundary away from the aperture and
prevents a zero-valued pupil edge from creating spurious high-frequency power.
`smoothingFwhm` and `lowPassCutoff` are mutually exclusive.

For a hybrid primary basis, select an inclusive Noll-Zernike range followed by
mxlib's rectangular modified-Fourier sequence:

```ini
[basis]
type=hybrid
hybridZernikeMin=2
hybridZernikeMax=10
hybridFourierN=30
```

The Fourier portion uses `mx::sigproc::makeFourierModeFreqs_Rect` followed by
`makeModifiedFourierMode`, so its plane count is determined by the rectangular
frequency cutoff, not by `basis.modes`. With `hybridFourierN=30`, the hybrid
primary basis has nine Zernikes (Noll 2--10) and 960 Fourier modes, for 969
planes. The actual count is written as `basis.primaryModes` in the run
metadata; `basis.cutoffs=default` includes the final generated count.

`run_gmagaox_stroke_exploration_hybrid.sh` runs the corresponding full
comparison under `output/gmagaox_stroke_exploration_hybrid` by default. Its
hybrid cutoffs retain the Zernike-only control, then increase the rectangular
Fourier linear DOF by 2 from `8` through `32`. This gives primary-mode cutoffs
`0,9,89,129,177,233,297,369,449,537,633,737,849,969,1097` for Zernikes Noll
2--10 followed by the Fourier sequence. The 6.5 m reference cases remain
Zernike bases and use `0,5,10,20,40,60,80,100` cutoffs.

For the actual `basis.fit=pinv` decomposition, set
`output.writePinvSvdBasis=true`. This writes one cube per cutoff as
`basis_pinv_svd_<label>_<cutoff>modes.fits`. Its planes are the pupil-domain
left singular modes used by the pseudo-inverse, in descending singular-value
order. Components rejected by the configured pseudo-inverse are zero planes.
Unlike the Gram--Schmidt diagnostic, these cubes depend on each cutoff.

To inspect a set of residual realizations across every requested modal cutoff,
set the number of frames to retain per cutoff:

```bash
./apertureStroke -c exploration_gmt_25m_segmentPTT.conf \
    --output.residualCubeFrames=100 \
    --output.directory=output/gmagaox_stroke_exploration/gmt_25m_segmentPTT/L0inf/seeing0p79
```

This writes one cube per cutoff as
`phase_residual_cube_<label>_<cutoff>modes_<run-tag>.fits`. Each cube contains
the requested number of trials (or all trials when fewer were run), with one
masked residual-phase image per plane. FITS values are phase in radians,
matching the existing `phase_input` and `phase_residual` products. A value of
zero, the default, disables these cubes.

To find where each residual-cube frame reaches its maximum inter-actuator
stroke, use the pupil mask from the simulation:

```bash
python3 max_interactuator_stroke.py \
    phase_residual_cube_gmt_25m_100modes_<run-tag>.fits \
    /home/jrmales/Source/mxWork/GMT/pupil/gmt-pupil-182x182.fits
```

The output gives zero-based `row col` coordinates for both pixels in the
maximum-difference pair, matching `maxAbsPixelDiff` in `apertureStroke`.

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
are the predicted mean absolute modal P2V and true RMS, matching the `mean`
and `rms` columns in the measured `zernike_p2v` tables. Column 4 is the
standard deviation about the mean. Seeing is
converted to r0 at its reference wavelength and then scaled to the configured
phase wavelength before applying the Noll variance.

The utility intentionally rejects FITS pupils, central obscurations, nonzero
outer scale, and non-Zernike bases because those are outside the circular Noll
prediction.
