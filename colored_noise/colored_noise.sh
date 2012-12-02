#!/bin/bash
#
cp colored_noise.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include colored_noise.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling colored_noise.cpp"
  exit
fi
rm compiler.txt
#
mv colored_noise.o ~/libcpp/$ARCH/colored_noise.o
#
echo "Library installed as ~/libcpp/$ARCH/colored_noise.o"
