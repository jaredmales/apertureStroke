# GMT–MagAO-X Stroke Exploration

## Selected configuration

This report records the basis-conditioning study for the 25.4 m GMT pupil and
the configuration selected for the full `gmagaox_stroke_exploration` run.
The residual diagnostic is the maximum absolute difference between horizontally
or vertically adjacent valid pupil pixels, reported in reflective-surface
microns.

```ini
[basis]
fit=pinv
pinvAlpha=0.001
pinvMaxCondition=1e6
smoothingFwhm=12.5
smoothingStart=100
lowPassCutoff=0
```

Gaussian smoothing applies only to primary Zernikes with zero-based index 100
and above; prefix modes are not smoothed. For the GMT segment PTT + Zernike
case, the 21 segment PTT modes remain at the front of each fit. Before
convolution, each primary mode is nearest-neighbor extrapolated beyond the
binary pupil, filtered, cropped, pupil-masked, and RMS-normalized. This avoids
creating high spatial frequency content by zero-padding the aperture edge.

The workstation sweep tags all 25.4 m output cases and files as `fwhm12p5`,
so they cannot be confused with earlier unfiltered or FWHM 5, 10, and 15
results.

## Diagnostic development

- `output.residualCubeFrames` saves a chosen number of pupil-masked residual
  phase frames at each cutoff. FITS values are phase radians, not surface
  microns.
- `max_interactuator_stroke.py` reproduces the production `maxAbsPixelDiff`
  search and reports the maximum pixel pair for every cube plane.
- Combined-basis, Gram--Schmidt, and pseudo-inverse SVD diagnostic cubes can
  be written without changing the selected pseudo-inverse fit.

Raw-basis, Gram--Schmidt, and residual-cube inspection located the problematic
high-frequency power at exterior segment edges. Gram--Schmidt cleans the low
orders, but high-order modes still accumulate edge-scale power.

## Filter sweep

All cases use 10,000 phase screens, 0.79 arcsec seeing at 500 nm, infinite
outer scale, subharmonic level 8, 0.8 um wavelength, the GMT segment-PTT
prefix, and `basis.fit=pinv`. At 1000 primary modes:

| Gaussian FWHM (pixels) | Mean P2P (um) | Maximum P2P (um) |
| ---: | ---: | ---: |
| 5.0 | 0.366330 | 0.815757 |
| 7.5 | 0.345813 | 0.726554 |
| 10.0 | 0.344800 | 0.661828 |
| **12.5** | **0.336349** | **0.624999** |
| 15.0 | 0.354558 | 0.700740 |

The optimum shifts with cutoff: FWHM 15 is best near 400--600 modes, whereas
FWHM 12.5 is best from 700--1000 modes. FWHM 15 therefore begins to suppress
the available high-order degrees of freedom before the 1000-mode cutoff.
FWHM 12.5 is the best sampled 1000-mode compromise.

A hard circular Fourier cutoff of 0.2 cycles/pixel, with a full pupil-width
nearest-neighbor extension before the FFT, gave 0.362568 um mean and 0.934866
um maximum P2P at 1000 modes. Its sharp transition was inferior to the
Gaussian rolloff, consistent with edge ringing.

## Hybrid Zernike--Fourier basis

The next basis replaces the high-order Zernikes with mxlib's rectangular
modified-Fourier sequence. It retains the segment PTT prefix and uses primary
Zernikes Noll 2--10 (tip/tilt through Noll 10), followed by the full sequence
from `mx::sigproc::makeFourierModeFreqs_Rect(30)`. This produces nine Zernikes
and 960 Fourier modes: 969 primary modes, or 990 modes including the 21
segment-PTT prefix modes. The terminal cutoff is therefore 969, not 1000.

No Gaussian or hard Fourier filter is applied to this hybrid basis. The
combined diagnostic basis is written to
`output/gmagaox_stroke_exploration/gmt_25m_segmentPTT_hybrid_z2to10_fourierRectN30/`.

| 10,000-screen comparison | Mean P2P (um) | Standard deviation (um) | Maximum P2P (um) |
| --- | ---: | ---: | ---: |
| FWHM 12.5 Zernike, 1000 primary modes | 0.336349 | 0.038375 | 0.624999 |
| **Hybrid, 969 primary modes** | **0.291286** | **0.024938** | **0.469965** |
| MagAO-X 6.5 m, 60 modes | 0.264785 | 0.029632 | 0.443668 |

Relative to the FWHM 12.5 Zernike result, the hybrid basis reduces mean P2P
by 13.4% and the observed maximum by 24.8%. It is within 10% in mean and 6%
in maximum of the matched 6.5 m/60-mode result. This strongly supports the
remaining tail being driven by the high-order Zernike edge behavior rather
than only by the larger GMT sampled area.

## Matched 6.5 m comparison

The approximate characteristic modal scales are

\[
25.4 / \sqrt{1000} = 0.803\ \mathrm{m}, \qquad
6.5 / \sqrt{60} = 0.839\ \mathrm{m}.
\]

GMT 1000 modes therefore corresponds roughly to 60--66 modes on the 6.5 m
pupil. The 6.5 m/60-mode reference has 0.264785 um mean and 0.443668 um
maximum P2P. The selected GMT/1000-mode case has 0.336349 um mean and
0.624999 um maximum P2P.

At fixed actuator pitch, opportunities for an extreme scale approximately
with pupil area:

\[
r = (25.4/6.5)^2 = 15.27.
\]

The `p2pdiff` files are histograms of the per-screen spatial maximum. Thus,
under matched local statistics and independent effective samples, the 6.5 m
histogram predicts the GMT maximum CDF as

\[
F_{\mathrm{GMT}}(t) = F_{6.5}(t)^r.
\]

For GMT/1000 modes with FWHM 12.5:

| Quantile of per-screen maximum P2P | Area-scaled 6.5 m prediction (um) | Observed GMT (um) |
| ---: | ---: | ---: |
| 50% | 0.3108 | 0.3158 |
| 90% | 0.3547 | 0.3752 |
| 99% | 0.4066 | 0.4562 |
| 99.9% | 0.4334 | 0.5337 |

Larger sampled GMT area explains most of the central increase. The remaining
high-tail excess is consistent with the exterior-edge basis behavior. It is
retained as a conservative state: a refined basis may reduce that tail, but
hardware manageability must also be preserved. The 10,000-screen 6.5 m sample
has no events above about 0.44 um, so beyond that observed tail a larger sample
or fitted-tail model is required.

## Workstation execution

From `apertureStroke`:

```bash
make
./run_gmagaox_stroke_exploration.sh
```

The optional positional arguments are trial count, worker count, and output
root:

```bash
./run_gmagaox_stroke_exploration.sh 10000 0 output/gmagaox_stroke_exploration
```

Set `FORCE=1` to rerun finished case directories. Set `SEGMENT_MODES` if the
GMT segment-PTT FITS cube is not at the configured default location. The
existing 6.5 m reference output is reused by default; set
`RUN_REFERENCE_CASES=1` only when it must be regenerated.

The 25.4 m cases are `circle_25m_fwhm12p5`, `gmt_25m_fwhm12p5`, and
`gmt_25m_segmentPTT_fwhm12p5`; the 6.5 m reference cases are unchanged.
