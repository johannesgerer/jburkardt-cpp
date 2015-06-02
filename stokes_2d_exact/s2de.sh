#!/bin/bash
#
cp s2de.hpp /$HOME/include
#
g++ -c -I/$HOME/include s2de.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling s2de.cpp"
  exit
fi
#
mv s2de.o ~/libcpp/$ARCH/s2de.o
#
echo "Library installed as ~/libcpp/$ARCH/s2de.o"
