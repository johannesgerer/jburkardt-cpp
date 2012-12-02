#!/bin/bash
#
g++ -c -g -I$HOME/include grid_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grid_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ grid_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grid_dataset.o."
  exit
fi
#
rm grid_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/grid_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/grid_dataset"
