#!/bin/bash
#
cp wathen.hpp /$HOME/include
#
g++ -c -I/$HOME/include wathen.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wathen.cpp"
  exit
fi
#
mv wathen.o ~/libcpp/$ARCH/wathen.o
#
echo "Library installed as ~/libcpp/$ARCH/wathen.o"
