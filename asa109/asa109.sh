#!/bin/bash
#
cp asa109.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa109.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa109.cpp"
  exit
fi
#
mv asa109.o ~/libcpp/$ARCH/asa109.o
#
echo "Library installed as ~/libcpp/$ARCH/asa109.o"
