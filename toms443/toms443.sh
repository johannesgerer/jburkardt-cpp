#!/bin/bash
#
cp toms443.hpp /$HOME/include
#
g++ -c -I/$HOME/include toms443.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling toms443.cpp"
  exit
fi
#
mv toms443.o ~/libcpp/$ARCH/toms443.o
#
echo "Library installed as ~/libcpp/$ARCH/toms443.o"
