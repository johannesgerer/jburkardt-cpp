#!/bin/bash
#
g++ -c -g -I/$HOME/include ncc_tetrahedron_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ncc_tetrahedron_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ncc_tetrahedron_prb.o /$HOME/libcpp/$ARCH/ncc_tetrahedron.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ncc_tetrahedron_prb.o."
  exit
fi
#
rm ncc_tetrahedron_prb.o
#
mv a.out ncc_tetrahedron_prb
./ncc_tetrahedron_prb > ncc_tetrahedron_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ncc_tetrahedron_prb."
  exit
fi
rm ncc_tetrahedron_prb
#
echo "Program output written to ncc_tetrahedron_prb_output.txt"
