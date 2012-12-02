#!/bin/bash
#
g++ -c -g fem2d_navier_stokes.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_navier_stokes.cpp"
  exit
fi
rm compiler.txt
#
mv fem2d_navier_stokes.o ~/libcpp/$ARCH
#
echo "The fem2d_navier_stokes partial program has been created."
