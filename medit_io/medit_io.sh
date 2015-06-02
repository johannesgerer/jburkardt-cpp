#!/bin/bash
#
cp medit_io.hpp /$HOME/include
#
g++ -c -I /$HOME/include medit_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling medit_io.cpp"
  exit
fi
#
mv medit_io.o ~/libcpp/$ARCH/medit_io.o
#
echo "Library installed as ~/libcpp/$ARCH/medit_io.o"
