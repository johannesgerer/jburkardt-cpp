#!/bin/bash
#
cp asa245.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa245.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa245.cpp"
  exit
fi
#
mv asa245.o ~/libcpp/$ARCH/asa245.o
#
echo "Library installed as ~/libcpp/$ARCH/asa245.o"
