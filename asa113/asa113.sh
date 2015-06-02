#!/bin/bash
#
cp asa113.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa113.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa113.cpp"
  exit
fi
#
mv asa113.o ~/libcpp/$ARCH/asa113.o
#
echo "Library installed as ~/libcpp/$ARCH/asa113.o"
