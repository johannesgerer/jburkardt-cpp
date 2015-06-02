#!/bin/bash
#
cp asa136.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa136.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa136.cpp"
  exit
fi
#
mv asa136.o ~/libcpp/$ARCH/asa136.o
#
echo "Library installed as ~/libcpp/$ARCH/asa136.o"
