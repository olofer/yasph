#!/bin/bash

# Confined low-pressure ideal gas in a gravitational field
# usage: ./demo-gravity.sh numpergas

# example: ../scripts/demo-gravity.sh 500 2e5

# 0.1786 5.851

python3 ../scripts/yasph-init-ball.py --num-per-gas $1 --particle-file gdemo-parts.yasph --parameter-file gdemo-params.yasph --lattice --barrier-radius 1.00 --density 0.90 1.78 --pressure 20 20
./yasph gdemo-parts.yasph gdemo-params.yasph steps=0 param-file=gdemo-params-explicit.yasph frame-steps=100 frame-file=gravity.yasph viscosity=1 epsbn=0.0 epsbt=0.0 trace-steps=100 trace-file=gravity-trace.yasph verbosity=0

# Pass in a comprehensive set of parameters (and additionally override a few options on command line; typical usage)
./yasph gdemo-parts.yasph gdemo-params-explicit.yasph steps=$2 dt=2.0e-5 threads=4 verbosity=1 gravity=0.0,-9.82 > gravity-log.yasph 
