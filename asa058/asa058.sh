#!/bin/bash
#
cp asa058.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa058.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa058.cpp"
  exit
fi
#
mv asa058.o ~/libcpp/$ARCH/asa058.o
#
echo "Library installed as ~/libcpp/$ARCH/asa058.o"
