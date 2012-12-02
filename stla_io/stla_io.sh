#!/bin/bash
#
cp stla_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include stla_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stla_io.cpp"
  exit
fi
rm compiler.txt
#
mv stla_io.o ~/libcpp/$ARCH/stla_io.o
#
echo "Library installed as ~/libcpp/$ARCH/stla_io.o"
