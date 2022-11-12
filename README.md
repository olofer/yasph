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

This should take about a minute. Then open `browser/frame-player.html` in a browser. Select the newly created file `barriers.yasph` and watch the particles evolve.
