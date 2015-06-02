#!/bin/bash
#
cp linplus.hpp /$HOME/include
#
g++ -c -I /$HOME/include linplus.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling linplus.cpp"
  exit
fi
#
mv linplus.o ~/libcpp/$ARCH/linplus.o
#
echo "Library installed as ~/libcpp/$ARCH/linplus.o"
