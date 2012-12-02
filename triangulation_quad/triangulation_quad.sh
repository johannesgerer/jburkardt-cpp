#!/bin/bash
#
g++ -c -g triangulation_quad.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_quad.cpp"
  exit
fi
rm compiler.txt
#
g++ triangulation_quad.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_quad.o."
  exit
fi
#
rm triangulation_quad.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangulation_quad
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_quad"
