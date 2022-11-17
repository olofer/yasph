# yasph
Plain `C` smooth particle hydrodynamics in 2D using basic `OpenMP`.

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
In the `build` folder run the script `../scripts/demo.sh` as a self-contained illustration of how to run and use the program. Run the browser program on the output file `octagon.yashp` to playback simulation results.

### Diagnostics trace files
The browser program `browser/trace-viewer.html` can be used to view time traces of column data in the output files specified with the `trace-file` (and `trace-steps`) options. These files are basic CSV files. 

## Useful references
TBD.
