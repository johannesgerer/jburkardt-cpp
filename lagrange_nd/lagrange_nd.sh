#!/bin/bash
#
cp lagrange_nd.hpp /$HOME/include
#
g++ -c -I /$HOME/include lagrange_nd.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_nd.cpp"
  exit
fi
#
mv lagrange_nd.o ~/libcpp/$ARCH/lagrange_nd.o
#
echo "Library installed as ~/libcpp/$ARCH/lagrange_nd.o"
