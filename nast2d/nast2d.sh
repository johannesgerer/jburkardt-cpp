#!/bin/bash
#
g++ -c -g main.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling main.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -g boundary.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling boundary.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -g extras.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling extras.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -g init.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling init.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -g surface.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling surface.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -g uvp.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uvp.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -g visual.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling visual.cpp"
  exit
fi
rm compiler.txt
#
g++ main.o boundary.o extras.o init.o surface.o uvp.o visual.o -lm
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the nast2d object files."
  exit
fi
rm *.o
mv a.out ~/bincpp/$ARCH/nast2d
#
echo "Executable installed as ~/bincpp/$ARCH/nast2d"
