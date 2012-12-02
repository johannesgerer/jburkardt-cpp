#!/bin/bash
#
cp fem3d_pack.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include fem3d_pack.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem3d_pack.cpp"
  exit
fi
rm compiler.txt
#
mv fem3d_pack.o ~/libcpp/$ARCH/fem3d_pack.o
#
echo "Library installed as ~/libcpp/$ARCH/fem3d_pack.o"
