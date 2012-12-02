#!/bin/bash
#
g++ -c -g -I$HOME/include hex_grid_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hex_grid_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ hex_grid_dataset.o $HOME/libcpp/$ARCH/hex_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hex_grid_dataset.o."
  exit
fi
#
rm hex_grid_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/hex_grid_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/hex_grid_dataset"
