#!/bin/bash
#
cp chebyshev.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include chebyshev.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev.cpp"
  exit
fi
rm compiler.txt
#
mv chebyshev.o ~/libcpp/$ARCH/chebyshev.o
#
echo "Library installed as ~/libcpp/$ARCH/chebyshev.o"
