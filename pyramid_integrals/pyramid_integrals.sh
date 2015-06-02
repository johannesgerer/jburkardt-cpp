#!/bin/bash
#
cp pyramid_integrals.hpp /$HOME/include
#
g++ -c -I /$HOME/include pyramid_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv pyramid_integrals.o ~/libcpp/$ARCH/pyramid_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/pyramid_integrals.o"
