#!/bin/bash
#
g++ -c -g minesweeper.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling minesweeper.cpp"
  exit
fi
rm compiler.txt
#
g++ minesweeper.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading minesweeper.o"
  exit
fi
rm minesweeper.o
#
mv a.out ~/bincpp/$ARCH/minesweeper
#
echo "Executable installed as ~/bincpp/$ARCH/minesweeper"
