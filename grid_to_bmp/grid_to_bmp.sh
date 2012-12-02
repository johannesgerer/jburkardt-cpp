#!/bin/bash
#
g++ -c grid_to_bmp.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grid_to_bmp.cpp"
  exit
fi
rm compiler.txt
#
g++ grid_to_bmp.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grid_to_bmp.o."
  exit
fi
#
rm grid_to_bmp.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/grid_to_bmp
#
echo "Executable installed as ~/bincpp/$ARCH/grid_to_bmp"
