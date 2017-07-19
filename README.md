# n253lines
NGC 253 LineID project

To grab a copy of this project from git:

    git clone https://github.com/astroumd/n253lines

Currently this is in a private GIT repo, so you will need git 1.7.10
or above to download this with the right credentials. 
There is still a question if we can commit in 1.7.1 even though
the repo is now public.



Example commands to test scripts:

    ./line_matching.py testCubeSpectrum.tab ngc253_lines.list

    ./cubespectrum2.py ngc6503.cube.fits
    ./cubespectrum2.py ngc6503.cube.fits  244 182
    ./cubespectrum2.py ngc6503.cube.fits  161 128
    ./cubespectrum2.py ngc6503.cube.fits  77 85

Some N253 examples

     # will pick the reference pixel (167,167)
     ./cubespectrum2.py ngc253_fullcube_compact_spw1_clean.ce.fits 

     # pick another more interesting (?) point
     ./cubespectrum2.py ngc253_fullcube_compact_spw1_clean.ce.fits 138 183

     # now switch to plotting in km/s, using an arbitrary 113.0 as the restfreq
     ./cubespectrum2.py ngc253_fullcube_compact_spw1_clean.ce.fits 138 183 113.0

     # now switch to plotting in km/s, using 113.0 as the restfreq and from -500 to 500 km/s
     ./cubespectrum2.py ngc253_fullcube_compact_spw1_clean.ce.fits 138 183 113.0 -500 500

     # of course the km/s scale is wrong, the galaxy doesn't have a velocity of -300 as it
     shows in this plot.   So, to get a better idea, use the line matching now to get the velocity

     # use the last output from cubespectrum2.py in line_matching
     ./line_matching.py Frequency_Flux.tab ngc253_lines.list

     # if you match the two strong lines on the right, it would result in roughly 197 km/s
     # if you want to match the two weaker on the left, 450 km/s may seem better
     #    but now the two strong lines don't match at all.
     # Note however, that the C17O line *will* match nicely with the strong CN lines
     #
     # the complication will occur if spectra are not "simple gaussian", wihch indicates
     # line of sight and/or beam smearing issues for which there is no good solution.

     # repeat this process for spw's 0, 2 and 3 and see if you can agree that 197 is a good velocity
     # for this position.
