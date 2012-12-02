#!/bin/bash
#
cp ppmb_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include ppmb_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppmb_io.cpp"
  exit
fi
rm compiler.txt
#
mv ppmb_io.o ~/libcpp/$ARCH/ppmb_io.o
#
echo "Library installed as ~/libcpp/$ARCH/ppmb_io.o"
