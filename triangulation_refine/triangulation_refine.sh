#!/bin/bash
#
g++ -c -g triangulation_refine.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_refine.cpp"
  exit
fi
rm compiler.txt
#
g++ triangulation_refine.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_refine.o."
  exit
fi
#
rm triangulation_refine.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangulation_refine
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_refine"
