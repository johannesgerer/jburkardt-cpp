#!/bin/bash
#
cp fem_io.hpp /$HOME/include
#
g++ -c -I /$HOME/include fem_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_io.cpp"
  exit
fi
#
mv fem_io.o ~/libcpp/$ARCH/fem_io.o
#
echo "Library installed as ~/libcpp/$ARCH/fem_io.o"
