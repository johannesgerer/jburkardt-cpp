#!/bin/bash
#
cp hb_io.hpp /$HOME/include
#
g++ -c -I /$HOME/include hb_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hb_io.cpp"
  exit
fi
#
mv hb_io.o ~/libcpp/$ARCH/hb_io.o
#
echo "Library installed as ~/libcpp/$ARCH/hb_io.o"
