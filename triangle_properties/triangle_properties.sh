#!/bin/bash
#
g++ -c -g triangle_properties.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_properties.cpp"
  exit
fi
rm compiler.txt
#
g++ triangle_properties.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_properties.o"
  exit
fi
rm triangle_properties.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangle_properties
#
echo "Executable installed as ~/bincpp/$ARCH/triangle_properties"
