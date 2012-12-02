#!/bin/bash
#
cp toms462.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include toms462.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms462.cpp"
  exit
fi
rm compiler.txt
#
mv toms462.o ~/libcpp/$ARCH/toms462.o
#
echo "Library installed as ~/libcpp/$ARCH/toms462.o"
