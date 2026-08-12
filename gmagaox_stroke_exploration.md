This command: `./apertureStroke -c circle_25m.conf     --atmosphere.subharmonicLevel=8     --atmosphere.psdSubtractPiston=false     --atmosphere.psdSubtractTipTilt=false  --output.directory=output/circle_25m_sh8_os0_svd_rawpsd --atmosphere.subharmonicPrecalc=false --simulation.oversize=0` produced results just as good as with oversize=3, and ran extremely fast on a good workstation with many cores.  This sets the baseline for how we can do the following experiment.

Goal: establish require stroke for a AO system on the 25.4 m GMT (GMagAO-X), with comparison to an AO system on a 6.5 m telescope (MagAO-X) with the same site.  

We want to first establish the theoretical baseline with circular unobstructed apertures and compare to Noll theory for infinite outscale, and show the reduction in required stroke for smaller outerscale.

Then we will move to the real apertures.  
- For MagAO-X this is central obscuration of 0.29 on a 48 pixel diameter, we won't worry about spiders.  Test 0, 20, 40, 60, 80, and 100 modes.
- For GMagAO-X use the 182x182 FITS file pupil. Test 0,20,100-1000 modes in steps of 100.  Do with and without segment modes.

All of this is against seeing of 0.5,0.64,0.79 at 500 nm.

Outerscales of 25, 50, 100, 1000, meters and inf.

We don't need fourier modes.  

Use the SVD fit.  
