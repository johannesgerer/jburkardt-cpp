#!/bin/bash
#
cp dream.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include dream.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dream.cpp"
  exit
fi
rm compiler.txt
#
mv dream.o ~/libcpp/$ARCH/dream.o
#
echo "Library installed as ~/libcpp/$ARCH/dream.o"
