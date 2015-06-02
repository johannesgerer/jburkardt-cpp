#!/bin/bash
#
cp ffmsh_io.hpp /$HOME/include
#
g++ -c -I/$HOME/include ffmsh_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ffmsh_io.cpp"
  exit
fi
#
mv ffmsh_io.o ~/libcpp/$ARCH/ffmsh_io.o
#
echo "Library installed as ~/libcpp/$ARCH/ffmsh_io.o"
