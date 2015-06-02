#!/bin/bash
#
cp latin_random.hpp /$HOME/include
#
g++ -c -I /$HOME/include latin_random.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_random.cpp"
  exit
fi
#
mv latin_random.o ~/libcpp/$ARCH/latin_random.o
#
echo "Library installed as ~/libcpp/$ARCH/latin_random.o"
