#!/bin/bash
#
g++ -c poisson.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson.cpp"
  exit
fi
rm compiler.txt
#
g++ poisson.o
if [ $? -ne 0 ]; then
  echo "Errors loading poisson.o"
  exit
fi
rm poisson.o
#
mv a.out poisson
./poisson > poisson_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poisson"
  exit
fi
rm poisson
#
echo "Program output written to poisson_output.txt"
