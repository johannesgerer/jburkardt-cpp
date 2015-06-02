#!/bin/bash
#
cp fem2d_bvp_serene.hpp /$HOME/include
#
g++ -c -I /$HOME/include fem2d_bvp_serene.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_serene.cpp"
  exit
fi
#
mv fem2d_bvp_serene.o ~/libcpp/$ARCH/fem2d_bvp_serene.o
#
echo "Library installed as ~/libcpp/$ARCH/fem2d_bvp_serene.o"
