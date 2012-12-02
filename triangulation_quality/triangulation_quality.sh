#!/bin/bash
#
g++ -c -g triangulation_quality.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_quality.cpp"
  exit
fi
rm compiler.txt
#
g++ triangulation_quality.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_quality.o."
  exit
fi
#
rm triangulation_quality.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangulation_quality
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_quality"
