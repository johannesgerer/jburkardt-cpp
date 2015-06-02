#!/bin/bash
#
cp i4lib.hpp /$HOME/include
#
g++ -c -I /$HOME/include i4lib.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling i4lib.cpp"
  exit
fi
#
mv i4lib.o ~/libcpp/$ARCH/i4lib.o
#
echo "Library installed as ~/libcpp/$ARCH/i4lib.o"
