#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex05.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex05.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex05.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex05.o."
  exit
fi
#
rm dislin_ex05.o
#
mv a.out dislin_ex05
./dislin_ex05 > dislin_ex05_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex05."
  exit
fi
rm dislin_ex05
#
echo "Program output written to dislin_ex05_output.txt"
