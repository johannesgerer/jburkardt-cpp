#!/bin/bash
#
cp bisection_rc.hpp /$HOME/include
#
g++ -c -I/$HOME/include bisection_rc.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling bisection_rc.cpp"
  exit
fi
#
mv bisection_rc.o ~/libcpp/$ARCH/bisection_rc.o
#
echo "Library installed as ~/libcpp/$ARCH/bisection_rc.o"
