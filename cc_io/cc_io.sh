#!/bin/bash
#
cp cc_io.hpp /$HOME/include
#
g++ -c -I /$HOME/include cc_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cc_io.cpp"
  exit
fi
#
mv cc_io.o ~/libcpp/$ARCH/cc_io.o
#
echo "Library installed as ~/libcpp/$ARCH/cc_io.o"
