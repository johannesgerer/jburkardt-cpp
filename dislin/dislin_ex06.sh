#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex06.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex06.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex06.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex06.o."
  exit
fi
#
rm dislin_ex06.o
#
mv a.out dislin_ex06
./dislin_ex06 > dislin_ex06_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex06."
  exit
fi
rm dislin_ex06
#
echo "Program output written to dislin_ex06_output.txt"
