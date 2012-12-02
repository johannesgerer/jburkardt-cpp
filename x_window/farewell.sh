#!/bin/bash
#
g++ -c farewell.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling farewell.cpp"
  exit
fi
rm compiler.txt
#
g++ farewell.o -L/usr/X11R6/lib -lXt -lXaw
if [ $? -ne 0 ]; then
  echo "Errors linking and loading farewell.o"
  exit
fi
#
rm farewell.o
mv a.out ~/bincpp/$ARCH/farewell
#
echo "Executable installed as ~/bincpp/$ARCH/farewell"
