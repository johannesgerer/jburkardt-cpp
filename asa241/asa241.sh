#!/bin/bash
#
cp asa241.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa241.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa241.cpp"
  exit
fi
#
mv asa241.o ~/libcpp/$ARCH/asa241.o
#
echo "Library installed as ~/libcpp/$ARCH/asa241.o"
