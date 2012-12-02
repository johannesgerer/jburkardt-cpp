#!/bin/bash
#
g++ -c -g -I$HOME/include sparse_grid_open_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_open_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ sparse_grid_open_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_open_dataset.o."
  exit
fi
#
rm sparse_grid_open_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/sparse_grid_open_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/sparse_grid_open_dataset"
