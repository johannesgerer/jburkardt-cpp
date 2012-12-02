#!/bin/bash
#
g++ -c -g -I/$HOME/include tetrahedron_monte_carlo_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_monte_carlo_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ tetrahedron_monte_carlo_prb.o /$HOME/libcpp/$ARCH/tetrahedron_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_monte_carlo_prb.o."
  exit
fi
#
rm tetrahedron_monte_carlo_prb.o
#
mv a.out tetrahedron_monte_carlo_prb
./tetrahedron_monte_carlo_prb > tetrahedron_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tetrahedron_monte_carlo_prb."
  exit
fi
rm tetrahedron_monte_carlo_prb
#
echo "Program output written to tetrahedron_monte_carlo_prb_output.txt"
