#!/bin/bash
#
cp square_exactness.hpp /$HOME/include
#
g++ -c -I/$HOME/include square_exactness.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_exactness.cpp"
  exit
fi
#
mv square_exactness.o ~/libcpp/$ARCH/square_exactness.o
#
echo "Library installed as ~/libcpp/$ARCH/square_exactness.o"
