#!/bin/bash
#
cp chebyshev_series.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include chebyshev_series.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_series.cpp"
  exit
fi
rm compiler.txt
#
mv chebyshev_series.o ~/libcpp/$ARCH/chebyshev_series.o
#
echo "Library installed as ~/libcpp/$ARCH/chebyshev_series.o"
