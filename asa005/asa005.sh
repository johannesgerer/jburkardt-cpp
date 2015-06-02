#!/bin/bash
#
cp asa005.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa005.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa005.cpp"
  exit
fi
#
mv asa005.o ~/libcpp/$ARCH/asa005.o
#
echo "Library installed as ~/libcpp/$ARCH/asa005.o"
