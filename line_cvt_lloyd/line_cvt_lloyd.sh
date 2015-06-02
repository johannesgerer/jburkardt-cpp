#!/bin/bash
#
cp line_cvt_lloyd.hpp /$HOME/include
#
g++ -c -I/$HOME/include line_cvt_lloyd.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_cvt_lloyd.cpp"
  exit
fi
#
mv line_cvt_lloyd.o ~/libcpp/$ARCH/line_cvt_lloyd.o
#
echo "Library installed as ~/libcpp/$ARCH/line_cvt_lloyd.o"
