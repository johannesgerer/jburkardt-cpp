#!/bin/bash
#
cp exactness.hpp /$HOME/include
#
g++ -c -I/$HOME/include exactness.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling exactness.cpp"
  exit
fi
rm compiler.txt
#
mv exactness.o ~/libcpp/$ARCH/exactness.o
#
echo "Library installed as ~/libcpp/$ARCH/exactness.o"
