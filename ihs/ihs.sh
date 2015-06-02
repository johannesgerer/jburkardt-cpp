#!/bin/bash
#
cp ihs.hpp /$HOME/include
#
g++ -c -I /$HOME/include ihs.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ihs.cpp"
  exit
fi
#
mv ihs.o ~/libcpp/$ARCH/ihs.o
#
echo "Library installed as ~/libcpp/$ARCH/ihs.o"
