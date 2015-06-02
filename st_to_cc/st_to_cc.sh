#!/bin/bash
#
cp st_to_cc.hpp /$HOME/include
#
g++ -c -I/$HOME/include st_to_cc.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling st_to_cc.cpp"
  exit
fi
#
mv st_to_cc.o ~/libcpp/$ARCH/st_to_cc.o
#
echo "Library installed as ~/libcpp/$ARCH/st_to_cc.o"
