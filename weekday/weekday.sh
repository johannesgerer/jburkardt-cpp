#!/bin/bash
#
cp weekday.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include weekday.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling weekday.cpp"
  exit
fi
rm compiler.txt
#
mv weekday.o ~/libcpp/$ARCH/weekday.o
#
echo "Library installed as ~/libcpp/$ARCH/weekday.o"
