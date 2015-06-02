#!/bin/bash
#
cp asa172.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa172.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa172.cpp"
  exit
fi
#
mv asa172.o ~/libcpp/$ARCH/asa172.o
#
echo "Library installed as ~/libcpp/$ARCH/asa172.o"
