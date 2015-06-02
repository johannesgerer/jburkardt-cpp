#!/bin/bash
#
cp hermite_cubic.hpp /$HOME/include
#
g++ -c -I /$HOME/include hermite_cubic.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_cubic.cpp"
  exit
fi
#
mv hermite_cubic.o ~/libcpp/$ARCH/hermite_cubic.o
#
echo "Library installed as ~/libcpp/$ARCH/hermite_cubic.o"
