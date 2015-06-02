#!/bin/bash
#
cp asa047.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa047.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa047.cpp"
  exit
fi
#
mv asa047.o ~/libcpp/$ARCH/asa047.o
#
echo "Library installed as ~/libcpp/$ARCH/asa047.o"
