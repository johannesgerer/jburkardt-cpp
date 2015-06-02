#!/bin/bash
#
cp asa032.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa032.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa032.cpp"
  exit
fi
#
mv asa032.o ~/libcpp/$ARCH/asa032.o
#
echo "Library installed as ~/libcpp/$ARCH/asa032.o"
