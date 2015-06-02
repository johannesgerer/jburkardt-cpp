#!/bin/bash
#
cp spiral_data.hpp /$HOME/include
#
g++ -c -I/$HOME/include spiral_data.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling spiral_data.cpp"
  exit
fi
#
mv spiral_data.o ~/libcpp/$ARCH/spiral_data.o
#
echo "Library installed as ~/libcpp/$ARCH/spiral_data.o"
