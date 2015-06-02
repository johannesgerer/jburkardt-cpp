#!/bin/bash
#
cp circle_segment.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include circle_segment.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_segment.cpp"
  exit
fi
rm compiler.txt
#
mv circle_segment.o ~/libcpp/$ARCH/circle_segment.o
#
echo "Library installed as ~/libcpp/$ARCH/circle_segment.o"
