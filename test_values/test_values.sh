#!/bin/bash
#
cp test_values.hpp /$HOME/include
#
g++ -c -I /$HOME/include test_values.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling test_values.cpp"
  exit
fi
#
mv test_values.o ~/libcpp/$ARCH/test_values.o
#
echo "Library installed as ~/libcpp/$ARCH/test_values.o"
