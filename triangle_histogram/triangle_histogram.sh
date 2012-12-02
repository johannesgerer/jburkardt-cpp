#!/bin/bash
#
g++ -c -g triangle_histogram.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_histogram.cpp"
  exit
fi
rm compiler.txt
#
g++ triangle_histogram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_histogram.o."
  exit
fi
#
rm triangle_histogram.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangle_histogram
#
echo "Executable installed as ~/bincpp/$ARCH/triangle_histogram"
