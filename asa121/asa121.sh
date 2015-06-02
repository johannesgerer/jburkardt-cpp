#!/bin/bash
#
cp asa121.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa121.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa121.cpp"
  exit
fi
#
mv asa121.o ~/libcpp/$ARCH/asa121.o
#
echo "Library installed as ~/libcpp/$ARCH/asa121.o"
