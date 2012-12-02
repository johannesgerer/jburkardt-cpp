#!/bin/bash
#
g++ -c -g -I/$HOME/include nco_triangle_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nco_triangle_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ nco_triangle_prb.o /$HOME/libcpp/$ARCH/nco_triangle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nco_triangle_prb.o."
  exit
fi
#
rm nco_triangle_prb.o
#
mv a.out nco_triangle_prb
./nco_triangle_prb > nco_triangle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nco_triangle_prb."
  exit
fi
rm nco_triangle_prb
#
echo "Program output written to nco_triangle_prb_output.txt"
