#!/bin/bash

python3 generate_SG_disk.py 8 0.5 0.01

make
./rebound <<'EOF'
8
0.5
0.01
EOF
python3 plot_sim.py