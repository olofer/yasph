#!/bin/bash

# Low-pressure Helium and Xenon
# usage: ./demo-gravity.sh numpergas

python3 ../scripts/yasph-init-ball.py --num-per-gas $1 --particle-file gdemo-parts.yasph --parameter-file gdemo-params.yasph --barrier-radius 1.00 --density 0.1786 5.851 --pressure 10 10
./yasph gdemo-parts.yasph gdemo-params.yasph steps=0 param-file=gdemo-params-explicit.yasph frame-steps=50 frame-file=gravity.yasph viscosity=1 epsbn=0.0 epsbt=0.0 trace-steps=50 trace-file=gravity-trace.yasph verbosity=0

# Pass in a comprehensive set of parameters (and additionally override a few options on command line; typical usage)
./yasph gdemo-parts.yasph gdemo-params-explicit.yasph steps=5e4 dt=1.0e-5 threads=4 verbosity=1 gravity=0.0,-9.82 > gravity-log.yasph 
