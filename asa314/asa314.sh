#!/bin/bash
#
cp asa314.hpp /$HOME/include
#
g++ -c -I/$HOME/include asa314.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa314.cpp"
  exit
fi
#
mv asa314.o ~/libcpp/$ARCH/asa314.o
#
echo "Library installed as ~/libcpp/$ARCH/asa314.o"
