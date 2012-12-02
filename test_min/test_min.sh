#!/bin/bash
#
cp test_min.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_min.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_min.cpp"
  exit
fi
rm compiler.txt
#
mv test_min.o ~/libcpp/$ARCH/test_min.o
#
echo "Library installed as ~/libcpp/$ARCH/test_min.o"
