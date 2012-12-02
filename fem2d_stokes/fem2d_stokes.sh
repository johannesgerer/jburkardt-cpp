#!/bin/bash
#
g++ -c -g fem2d_stokes.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_stokes.cpp"
  exit
fi
rm compiler.txt
#
mv fem2d_stokes.o ~/libcpp/$ARCH
#
echo "Object code installed as ~/libcpp/$ARCH/fem2d_stokes.o"
