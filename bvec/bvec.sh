#!/bin/bash
#
cp bvec.hpp /$HOME/include
#
g++ -c -I /$HOME/include bvec.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling bvec.cpp"
  exit
fi
#
mv bvec.o ~/libcpp/$ARCH/bvec.o
#
echo "Library installed as ~/libcpp/$ARCH/bvec.o"
