#!/bin/bash
#
cp uniform.hpp /$HOME/include
#
g++ -c -I /$HOME/include uniform.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform.cpp"
  exit
fi
#
mv uniform.o ~/libcpp/$ARCH/uniform.o
#
echo "Library installed as ~/libcpp/$ARCH/uniform.o"
