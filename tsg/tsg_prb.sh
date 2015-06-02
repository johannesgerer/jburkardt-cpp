#!/bin/bash
#
g++ -O3 -ffast-math -fopenmp -mtune=native -otsg_prb -I/$HOME/include tsg_prb.cpp \
  -L/$HOME/libcpp/$ARCH -ltsg -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tsg_prb.o"
  exit
fi
#
./tsg_prb > tsg_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tsg_prb."
  exit
fi
rm tsg_prb
#
echo "Program output written to tsg_prb_output.txt"
