#!/bin/bash
#
cp fem2d_bvp_linear.hpp /$HOME/include
#
g++ -c -I /$HOME/include fem2d_bvp_linear.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_linear.cpp"
  exit
fi
#
mv fem2d_bvp_linear.o ~/libcpp/$ARCH/fem2d_bvp_linear.o
#
echo "Library installed as ~/libcpp/$ARCH/fem2d_bvp_linear.o"
