#!/bin/bash
#
g++ -c -g -I/$HOME/include nco_tetrahedron_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nco_tetrahedron_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ nco_tetrahedron_prb.o /$HOME/libcpp/$ARCH/nco_tetrahedron.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nco_tetrahedron_prb.o."
  exit
fi
#
rm nco_tetrahedron_prb.o
#
mv a.out nco_tetrahedron_prb
./nco_tetrahedron_prb > nco_tetrahedron_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nco_tetrahedron_prb."
  exit
fi
rm nco_tetrahedron_prb
#
echo "Program output written to nco_tetrahedron_prb_output.txt"
