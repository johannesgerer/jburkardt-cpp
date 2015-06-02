#!/bin/bash
#
cp asa006.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa006.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa006.cpp"
  exit
fi
#
mv asa006.o ~/libcpp/$ARCH/asa006.o
#
echo "Library installed as ~/libcpp/$ARCH/asa006.o"
