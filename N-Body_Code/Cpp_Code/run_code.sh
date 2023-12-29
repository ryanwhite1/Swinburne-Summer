#!/bin/bash

g++ -o nbody c_implementation.cpp -std=c++11
./nbody

python plot_sim.py