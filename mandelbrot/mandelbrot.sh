#!/bin/bash
#
g++ -c mandelbrot.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot.cpp"
  exit
fi
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
