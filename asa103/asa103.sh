#!/bin/bash
#
cp asa103.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa103.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa103.cpp"
  exit
fi
#
mv asa103.o ~/libcpp/$ARCH/asa103.o
#
echo "Library installed as ~/libcpp/$ARCH/asa103.o"
