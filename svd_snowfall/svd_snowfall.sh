#!/bin/bash
#
cp svd_snowfall.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include svd_snowfall.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_snowfall.cpp"
  exit
fi
rm compiler.txt
#
mv svd_snowfall.o ~/libcpp/$ARCH/svd_snowfall.o
#
echo "Library installed as ~/libcpp/$ARCH/svd_snowfall.o"
