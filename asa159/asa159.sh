#!/bin/bash
#
cp asa159.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa159.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa159.cpp"
  exit
fi
#
mv asa159.o ~/libcpp/$ARCH/asa159.o
#
echo "Library installed as ~/libcpp/$ARCH/asa159.o"
