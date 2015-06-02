#!/bin/bash
#
cp truncated_normal.hpp /$HOME/include
#
g++ -c -I /$HOME/include truncated_normal.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling truncated_normal.cpp"
  exit
fi
#
mv truncated_normal.o ~/libcpp/$ARCH/truncated_normal.o
#
echo "Library installed as ~/libcpp/$ARCH/truncated_normal.o"
