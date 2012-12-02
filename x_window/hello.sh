#!/bin/bash
#
g++ -c hello.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hello.cpp"
  exit
fi
rm compiler.txt
#
g++ hello.o -L/usr/X11R6/lib -lXm -lX11 -lXt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hello.o"
  exit
fi
#
rm hello.o
mv a.out ~/bincpp/$ARCH/hello
#
echo "Executable installed as ~/bincpp/$ARCH/hello"
