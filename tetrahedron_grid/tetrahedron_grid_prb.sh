#!/bin/bash
#
g++ -c -g tetrahedron_grid_prb.C >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_grid_prb.C."
  exit
fi
rm compiler.txt
#
g++ tetrahedron_grid_prb.o /$HOME/libcpp/$ARCH/tetrahedron_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_grid_prb.o."
  exit
fi
#
rm tetrahedron_grid_prb.o
#
mv a.out tetrahedron_grid_prb
./tetrahedron_grid_prb > tetrahedron_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tetrahedron_grid_prb."
  exit
fi
rm tetrahedron_grid_prb
#
echo "Program output written to tetrahedron_grid_prb_output.txt"
