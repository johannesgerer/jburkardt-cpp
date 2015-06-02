#!/bin/bash
#
cp asa152.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa152.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa152.cpp"
  exit
fi
#
mv asa152.o ~/libcpp/$ARCH/asa152.o
#
echo "Library installed as ~/libcpp/$ARCH/asa152.o"
