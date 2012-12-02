#!/bin/bash
#
g++ -c goodbye.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling goodbye.cpp"
  exit
fi
rm compiler.txt
#
g++ goodbye.o -L/usr/X11R6/lib -lXt -lXaw
if [ $? -ne 0 ]; then
  echo "Errors loading goodbye.o"
  exit
fi
#
rm goodbye.o
mv a.out ~/bin/$ARCH/goodbye
#
echo "Executable installed as ~/bin/$ARCH/goodbye"
