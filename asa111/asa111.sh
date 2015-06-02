#!/bin/bash
#
cp asa111.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa111.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa111.cpp"
  exit
fi
#
mv asa111.o ~/libcpp/$ARCH/asa111.o
#
echo "Library installed as ~/libcpp/$ARCH/asa111.o"
