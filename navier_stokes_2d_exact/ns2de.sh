#!/bin/bash
#
cp ns2de.hpp /$HOME/include
#
g++ -c -I/$HOME/include ns2de.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ns2de.cpp"
  exit
fi
#
mv ns2de.o ~/libcpp/$ARCH/ns2de.o
#
echo "Library installed as ~/libcpp/$ARCH/ns2de.o"
