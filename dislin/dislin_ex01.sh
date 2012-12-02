#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex01.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex01.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex01.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex01.o."
  exit
fi
#
rm dislin_ex01.o
#
mv a.out dislin_ex01
./dislin_ex01 > dislin_ex01_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex01."
  exit
fi
rm dislin_ex01
#
echo "Program output written to dislin_ex01_output.txt"
