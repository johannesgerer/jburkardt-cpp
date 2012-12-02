#!/bin/bash
#
g++ -c -g -I/$HOME/include sparse_grid_mixed_size_table.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_mixed_size_table.cpp"
  exit
fi
rm compiler.txt
#
g++ sparse_grid_mixed_size_table.o /$HOME/libcpp/$ARCH/sparse_grid_mixed.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_mixed_size_table.o."
  exit
fi
#
rm sparse_grid_mixed_size_table.o
#
mv a.out sparse_grid_mixed_size_table
./sparse_grid_mixed_size_table > sparse_grid_mixed_size_table_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_grid_mixed_size_table."
  exit
fi
rm sparse_grid_mixed_size_table
#
echo "Program output written to sparse_grid_mixed_size_table_output.txt"
