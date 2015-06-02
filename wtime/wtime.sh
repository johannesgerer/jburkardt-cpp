#!/bin/bash
#
cp wtime.hpp /$HOME/include
#
g++ -c -I /$HOME/include wtime.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime.cpp"
  exit
fi
#
mv wtime.o ~/libcpp/$ARCH/wtime.o
#
echo "Library installed as ~/libcpp/$ARCH/wtime.o"
