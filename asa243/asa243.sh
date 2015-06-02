#!/bin/bash
#
cp asa243.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa243.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa243.cpp"
  exit
fi
#
mv asa243.o ~/libcpp/$ARCH/asa243.o
#
echo "Library installed as ~/libcpp/$ARCH/asa243.o"
