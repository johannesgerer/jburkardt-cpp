#!/bin/bash
#
cp st_io.hpp /$HOME/include
#
g++ -c -I/$HOME/include st_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling st_io.cpp"
  exit
fi
#
mv st_io.o ~/libcpp/$ARCH/st_io.o
#
echo "Library installed as ~/libcpp/$ARCH/st_io.o"
