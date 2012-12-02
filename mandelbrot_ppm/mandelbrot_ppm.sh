#!/bin/bash
#
g++ -c -g mandelbrot_ppm.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot_ppm.cpp"
  exit
fi
rm compiler.txt
#
g++ mandelbrot_ppm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mandelbrot_ppm.o"
  exit
fi
rm mandelbrot_ppm.o
#
mv a.out ~/bincpp/$ARCH/mandelbrot_ppm
#
echo "Executable installed as ~/bincpp/$ARCH/mandelbrot_ppm"
