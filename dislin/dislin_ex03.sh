#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex03.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex03.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex03.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex03.o."
  exit
fi
#
rm dislin_ex03.o
#
mv a.out dislin_ex03
./dislin_ex03 > dislin_ex03_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex03."
  exit
fi
rm dislin_ex03
#
echo "Program output written to dislin_ex03_output.txt"
