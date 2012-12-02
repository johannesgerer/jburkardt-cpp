#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex13.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex13.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex13.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex13.o."
  exit
fi
#
rm dislin_ex13.o
#
mv a.out dislin_ex13
./dislin_ex13 > dislin_ex13_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex13."
  exit
fi
rm dislin_ex13
#
echo "Program output written to dislin_ex13_output.txt"
