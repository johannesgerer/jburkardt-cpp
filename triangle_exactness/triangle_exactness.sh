#!/bin/bash
#
g++ -c -g -I$HOME/include triangle_exactness.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_exactness.cpp"
  exit
fi
rm compiler.txt
#
g++ triangle_exactness.o -lm 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_exactness.o."
  exit
fi
#
rm triangle_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangle_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/triangle_exactness"
