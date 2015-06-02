#!/bin/bash
#
cp cg.hpp /$HOME/include
#
g++ -c -I/$HOME/include cg.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cg.cpp"
  exit
fi
#
mv cg.o ~/libcpp/$ARCH/cg.o
#
echo "Library installed as ~/libcpp/$ARCH/cg.o"
