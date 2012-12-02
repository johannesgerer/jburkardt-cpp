#!/bin/bash
#
cp spline.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include spline.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spline.cpp"
  exit
fi
rm compiler.txt
#
mv spline.o ~/libcpp/$ARCH/spline.o
#
echo "Library installed as ~/libcpp/$ARCH/spline.o"
