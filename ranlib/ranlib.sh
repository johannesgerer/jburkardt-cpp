#!/bin/bash
#
cp ranlib.hpp /$HOME/include
#
g++ -c -I/$HOME/include ranlib.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ranlib.cpp"
  exit
fi
#
mv ranlib.o ~/libcpp/$ARCH/ranlib.o
#
echo "Library installed as ~/libcpp/$ARCH/ranlib.o"
