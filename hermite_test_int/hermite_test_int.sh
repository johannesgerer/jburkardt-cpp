#!/bin/bash
#
cp hermite_test_int.hpp /$HOME/include
#
g++ -c -I /$HOME/include hermite_test_int.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_test_int.cpp"
  exit
fi
#
mv hermite_test_int.o ~/libcpp/$ARCH/hermite_test_int.o
#
echo "Library installed as ~/libcpp/$ARCH/hermite_test_int.o"
