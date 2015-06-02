#!/bin/bash
#
cp haar.hpp /$HOME/include
#
g++ -c -I /$HOME/include haar.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling haar.cpp"
  exit
fi
#
mv haar.o ~/libcpp/$ARCH/haar.o
#
echo "Library installed as ~/libcpp/$ARCH/haar.o"
