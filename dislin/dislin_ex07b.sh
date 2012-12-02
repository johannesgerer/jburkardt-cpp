#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex07b.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex07b.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex07b.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex07b.o."
  exit
fi
#
rm dislin_ex07b.o
#
mv a.out dislin_ex07b
./dislin_ex07b > dislin_ex07b_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex07b."
  exit
fi
rm dislin_ex07b
#
echo "Program output written to dislin_ex07b_output.txt"
