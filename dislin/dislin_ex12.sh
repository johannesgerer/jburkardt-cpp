#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex12.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex12.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex12.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex12.o."
  exit
fi
#
rm dislin_ex12.o
#
mv a.out dislin_ex12
./dislin_ex12 > dislin_ex12_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex12."
  exit
fi
rm dislin_ex12
#
echo "Program output written to dislin_ex12_output.txt"
