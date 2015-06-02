#!/bin/bash
#
cp edge.hpp /$HOME/include
#
g++ -c -I/$HOME/include edge.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling edge.cpp"
  exit
fi
#
mv edge.o ~/libcpp/$ARCH/edge.o
#
echo "Library installed as ~/libcpp/$ARCH/edge.o"
