#!/bin/bash
#
cp cube_integrals.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cube_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv cube_integrals.o ~/libcpp/$ARCH/cube_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/cube_integrals.o"
