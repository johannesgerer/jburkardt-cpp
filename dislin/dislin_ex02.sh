#!/bin/bash
#
g++ -c -I/usr/local/dislin dislin_ex02.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex02.cpp"
  exit
fi
rm compiler.txt
#
g++ dislin_ex02.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex02.o."
  exit
fi
#
rm dislin_ex02.o
#
mv a.out dislin_ex02
./dislin_ex02 > dislin_ex02_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex02."
  exit
fi
rm dislin_ex02
#
echo "Program output written to dislin_ex02_output.txt"
