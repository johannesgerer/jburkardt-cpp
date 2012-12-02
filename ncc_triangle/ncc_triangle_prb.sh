#!/bin/bash
#
g++ -c -g -I/$HOME/include ncc_triangle_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ncc_triangle_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ncc_triangle_prb.o /$HOME/libcpp/$ARCH/ncc_triangle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ncc_triangle_prb.o."
  exit
fi
#
rm ncc_triangle_prb.o
#
mv a.out ncc_triangle_prb
./ncc_triangle_prb > ncc_triangle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ncc_triangle_prb."
  exit
fi
rm ncc_triangle_prb
#
echo "Program output written to ncc_triangle_prb_output.txt"
