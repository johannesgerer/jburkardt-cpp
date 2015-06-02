#!/bin/bash
#
cp weekday.hpp /$HOME/include
#
g++ -c -I /$HOME/include weekday.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling weekday.cpp"
  exit
fi
#
mv weekday.o ~/libcpp/$ARCH/weekday.o
#
echo "Library installed as ~/libcpp/$ARCH/weekday.o"
