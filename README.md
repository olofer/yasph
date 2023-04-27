# yasph
Plain `C` smoothed particle hydrodynamics in 2D using basic `OpenMP`. For educational purposes. Employs either hash-based nearest-neighbour indexing or quad tree indexing (default is hash-based).

## Compile & test
Standard method:
```
mkdir build
cd build
cmake ..
make
ctest
```

## Demo & visualize
If compiled and tested and all seems well, do this (in the `build` folder):
```
./yasph state.yasph params.yasph dt=1.0e-6 steps=5e4 threads=4 frame-file=barriers.yasph frame-steps=50 epsbn=0.50 barrier-ball:0=0,0,3 barrier-ball:1=-1,-1,-1
```

This should take about a minute. Then open `browser/frame-player.html` in a browser. Select the newly created file `barriers.yasph` and watch the particles evolve (the browser program shows intructions at the bottom of its page on how you can interact with the playback).

### Demo not depending on post-test files
In the `build` folder run the script `../scripts/demo.sh 4000` as a self-contained illustration of how to run and use the program. Run the browser program on the output file `octagon.yasph` to playback simulation results.

### Diagnostics trace files
The browser program `browser/trace-viewer.html` can be used to view time traces of column data in the output files specified with the `trace-file` (and `trace-steps`) options. These files are basic CSV files. The above example script generates `octagon-trace.yasph`.

### Gravitational separation demo
The demonstration program `../scripts/demo-gravity.sh 1000 5e5` sets up a mixture of two low-pressure noble gases that separate in a gravitational field within a circular container. Use the visualization tools under `browser` to playback results. The simulation takes a while.

## Useful references

* J. J. Monaghan _Smoothed particle hydrodynamics_
https://doi.org/10.1146/annurev.aa.30.090192.002551

* D. J. Price _Smoothed particle hydrodynamics and magnetohydrodynamics_
https://doi.org/10.1016/j.jcp.2010.12.011
