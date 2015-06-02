#!/bin/bash
#
cp zero_rc.hpp /$HOME/include
#
g++ -c -I/$HOME/include zero_rc.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling zero_rc.cpp"
  exit
fi
#
mv zero_rc.o ~/libcpp/$ARCH/zero_rc.o
#
echo "Library installed as ~/libcpp/$ARCH/zero_rc.o"
