#!/bin/bash
#
g++ -c -g triangulation_histogram.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_histogram.cpp"
  exit
fi
rm compiler.txt
#
g++ triangulation_histogram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_histogram.o."
  exit
fi
#
rm triangulation_histogram.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangulation_histogram
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_histogram"
