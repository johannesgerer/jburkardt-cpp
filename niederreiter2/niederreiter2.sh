#!/bin/bash
#
cp niederreiter2.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include niederreiter2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter2.cpp"
  exit
fi
rm compiler.txt
#
mv niederreiter2.o ~/libcpp/$ARCH/niederreiter2.o
#
echo "Library installed as ~/libcpp/$ARCH/niederreiter2.o"
