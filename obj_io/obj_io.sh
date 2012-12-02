#!/bin/bash
#
cp obj_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include obj_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling obj_io.cpp"
  exit
fi
rm compiler.txt
#
mv obj_io.o ~/libcpp/$ARCH/obj_io.o
#
echo "Library installed as ~/libcpp/$ARCH/obj_io.o"
