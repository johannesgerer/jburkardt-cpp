#!/bin/bash
#
g++ -c -g -I/$HOME/include tri_surface_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tri_surface_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ tri_surface_io_prb.o /$HOME/libcpp/$ARCH/tri_surface_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tri_surface_io_prb.o."
  exit
fi
#
rm tri_surface_io_prb.o
#
mv a.out tri_surface_io_prb
./tri_surface_io_prb > tri_surface_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tri_surface_io_prb."
  exit
fi
rm tri_surface_io_prb
#
echo "Program output written to tri_surface_io_prb_output.txt"
