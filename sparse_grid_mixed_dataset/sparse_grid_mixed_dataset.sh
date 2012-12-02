#!/bin/bash
#
g++ -c -I/$HOME/include sparse_grid_mixed_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_mixed_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ sparse_grid_mixed_dataset.o /$HOME/libcpp/$ARCH/sparse_grid_mixed.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_mixed_dataset.o."
  exit
fi
#
rm sparse_grid_mixed_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/sparse_grid_mixed_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/sparse_grid_mixed_dataset"
