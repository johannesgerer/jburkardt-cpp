#!/bin/bash
#
cp asa053.hpp /$HOME/include
#
g++ -c -I/$HOME/include asa053.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa053.cpp"
  exit
fi
#
mv asa053.o ~/libcpp/$ARCH/asa053.o
#
echo "Library installed as ~/libcpp/$ARCH/asa053.o"
