#!/bin/bash
#
cp asa310.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa310.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa310.cpp"
  exit
fi
#
mv asa310.o ~/libcpp/$ARCH/asa310.o
#
echo "Library installed as ~/libcpp/$ARCH/asa310.o"
