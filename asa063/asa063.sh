#!/bin/bash
#
cp asa063.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa063.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa063.cpp"
  exit
fi
#
mv asa063.o ~/libcpp/$ARCH/asa063.o
#
echo "Library installed as ~/libcpp/$ARCH/asa063.o"
