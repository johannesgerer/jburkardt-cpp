#!/bin/bash
#
g++ -c -g mandelbrot.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot.cpp"
  exit
fi
rm compiler.txt
#
g++ mandelbrot.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mandelbrot.o"
  exit
fi
rm mandelbrot.o
#
mv a.out ~/bincpp/$ARCH/mandelbrot
#
echo "Executable installed as ~/bincpp/$ARCH/mandelbrot"
