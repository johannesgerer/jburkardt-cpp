#!/bin/bash
#
cp wtime.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include wtime.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime.cpp"
  exit
fi
rm compiler.txt
#
mv wtime.o ~/libcpp/$ARCH/wtime.o
#
echo "Library installed as ~/libcpp/$ARCH/wtime.o"
