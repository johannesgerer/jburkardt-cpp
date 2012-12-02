#!/bin/bash
#
cp brownian_motion_simulation.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include brownian_motion_simulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brownian_motion_simulation.cpp"
  exit
fi
rm compiler.txt
#
mv brownian_motion_simulation.o ~/libcpp/$ARCH/brownian_motion_simulation.o
#
echo "Library installed as ~/libcpp/$ARCH/brownian_motion_simulation.o"
