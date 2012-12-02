#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex04.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex04.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex04.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex04.o."
  exit
fi
#
rm dislin_ex04.o
#
mv a.out dislin_ex04
./dislin_ex04 > dislin_ex04_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex04."
  exit
fi
rm dislin_ex04
#
echo "Program output written to dislin_ex04_output.txt"
