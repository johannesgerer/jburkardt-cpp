#!/bin/bash
#
cp hammersley.hpp /$HOME/include
#
g++ -c -I /$HOME/include hammersley.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hammersley.cpp"
  exit
fi
#
mv hammersley.o ~/libcpp/$ARCH/hammersley.o
#
echo "Library installed as ~/libcpp/$ARCH/hammersley.o"
