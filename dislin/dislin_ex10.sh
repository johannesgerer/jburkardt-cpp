#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex10.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex10.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex10.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex10.o."
  exit
fi
#
rm dislin_ex10.o
#
mv a.out dislin_ex10
./dislin_ex10 > dislin_ex10_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex10."
  exit
fi
rm dislin_ex10
#
echo "Program output written to dislin_ex10_output.txt"
