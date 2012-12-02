#!/bin/bash
#
g++ -c -g shallow_water_1d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shallow_water_1d.cpp"
  exit
fi
rm compiler.txt
#
g++ shallow_water_1d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shallow_water_1d.o"
  exit
fi
rm shallow_water_1d.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/shallow_water_1d
#
echo "Program installed as ~/bincpp/$ARCH/shallow_water_1d"
