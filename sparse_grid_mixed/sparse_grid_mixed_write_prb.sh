#!/bin/bash
#
g++ -c -g -I/$HOME/include sparse_grid_mixed_write_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_mixed_write_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sparse_grid_mixed_write_prb.o /$HOME/libcpp/$ARCH/sparse_grid_mixed.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_mixed_write_prb.o."
  exit
fi
#
rm sparse_grid_mixed_write_prb.o
#
mv a.out sparse_grid_mixed_write_prb
./sparse_grid_mixed_write_prb > sparse_grid_mixed_write_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_grid_mixed_write_prb."
  exit
fi
rm sparse_grid_mixed_write_prb
#
echo "Program output written to sparse_grid_mixed_write_prb_output.txt"
#
#  Move sparse grid files to dataset directory.
#
mv *_a.txt ../../datasets/sparse_grid_mixed
mv *_b.txt ../../datasets/sparse_grid_mixed
mv *_r.txt ../../datasets/sparse_grid_mixed
mv *_w.txt ../../datasets/sparse_grid_mixed
mv *_x.txt ../../datasets/sparse_grid_mixed
#
echo "Program output files moved to ../../datasets/sparse_grid_mixed"
