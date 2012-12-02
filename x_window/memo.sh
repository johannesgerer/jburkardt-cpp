#!/bin/bash
#
g++ -c memo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling memo.cpp"
  exit
fi
rm compiler.txt
#
g++ memo.o -L/usr/X11R6/lib -lXt -lXaw
if [ $? -ne 0 ]; then
  echo "Errors linking and loading memo.o"
  exit
fi
#
rm memo.o
mv a.out ~/bincpp/$ARCH/memo
#
echo "Executable installed as ~/bincpp/$ARCH/memo"
