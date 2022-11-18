#!/bin/bash

# usage: ./demo.sh numparticles

# Generate initial set of particles and a parameters file
python3 ../scripts/yasph-init-box.py --n $1 --particle-file demo-parts.yasph --parameter-file demo-params.yasph --barrier-n 8 --barrier-radius 1.50 --xmin -1 --xmax 1 --ymin -0.5 --ymax 0.5 
./yasph demo-parts.yasph demo-params.yasph steps=0 param-file=demo-params-explicit.yasph frame-steps=50 frame-file=octagon.yasph viscosity=1 epsbn=0.0 epsbt=0.0 trace-steps=25 trace-file=octagon-trace.yasph verbosity=0

# Pass in a comprehensive set of parameters (and additionally override a few options on command line; typical usage)
./yasph demo-parts.yasph demo-params-explicit.yasph steps=2e4 dt=1.0e-6 threads=4 verbosity=1 > octagon-log.yasph 
