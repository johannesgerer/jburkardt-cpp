#!/bin/bash
#
cp pwl_interp_2d_scattered.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include pwl_interp_2d_scattered.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d_scattered.cpp"
  exit
fi
rm compiler.txt
#
mv pwl_interp_2d_scattered.o ~/libcpp/$ARCH/pwl_interp_2d_scattered.o
#
echo "Library installed as ~/libcpp/$ARCH/pwl_interp_2d_scattered.o"
