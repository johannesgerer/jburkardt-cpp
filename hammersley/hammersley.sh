#!/bin/bash
#
cp hammersley.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include hammersley.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hammersley.cpp"
  exit
fi
rm compiler.txt
#
mv hammersley.o ~/libcpp/$ARCH/hammersley.o
#
echo "Library installed as ~/libcpp/$ARCH/hammersley.o"
