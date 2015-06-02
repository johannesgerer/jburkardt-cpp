#!/bin/bash
#
g++ -c -I/$HOME/include polygon_triangulate_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_triangulate_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ polygon_triangulate_prb.o /$HOME/libcpp/$ARCH/polygon_triangulate.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygon_triangulate_prb.o."
  exit
fi
#
rm polygon_triangulate_prb.o
#
mv a.out polygon_triangulate_prb
./polygon_triangulate_prb > polygon_triangulate_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polygon_triangulate_prb."
  exit
fi
rm polygon_triangulate_prb
#
echo "Program output written to polygon_triangulate_prb_output.txt"
