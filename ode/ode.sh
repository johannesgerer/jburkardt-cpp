#!/bin/bash
#
cp ode.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include ode.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ode.cpp."
  exit
fi
rm compiler.txt
#
mv ode.o ~/libcpp/$ARCH/ode.o
#
echo "Library installed as ~/libcpp/$ARCH/ode.o"
