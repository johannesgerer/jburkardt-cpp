#!/bin/bash
#
cp hpp.hpp /$HOME/include
#
g++ -c -I /$HOME/include hpp.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hpp.cpp"
  exit
fi
#
mv hpp.o ~/libcpp/$ARCH/hpp.o
#
echo "Library installed as ~/libcpp/$ARCH/hpp.o"
