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


