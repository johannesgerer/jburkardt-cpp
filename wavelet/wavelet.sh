#!/bin/bash
#
cp wavelet.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include wavelet.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wavelet.cpp"
  exit
fi
rm compiler.txt
#
mv wavelet.o ~/libcpp/$ARCH/wavelet.o
#
echo "Library installed as ~/libcpp/$ARCH/wavelet.o"
