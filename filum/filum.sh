#!/bin/bash
#
cp filum.hpp /$HOME/include
#
g++ -c -I /$HOME/include filum.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling filum.cpp"
  exit
fi
#
mv filum.o ~/libcpp/$ARCH/filum.o
#
echo "Library installed as ~/libcpp/$ARCH/filum.o"
