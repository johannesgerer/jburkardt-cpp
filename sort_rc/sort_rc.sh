#!/bin/bash
#
cp sort_rc.hpp /$HOME/include
#
g++ -c -I/$HOME/include sort_rc.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sort_rc.cpp"
  exit
fi
#
mv sort_rc.o ~/libcpp/$ARCH/sort_rc.o
#
echo "Library installed as ~/libcpp/$ARCH/sort_rc.o"
