#!/bin/bash
#
g++ -c pcl_read.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pcl_read.cpp"
  exit
fi
rm compiler.txt
#
g++ pcl_read.o
if [ $? -ne 0 ]; then
  echo "Errors loading pcl_read.o."
  exit
fi
#
rm pcl_read.o
mv a.out ~/bincpp/$ARCH/pcl_read
#
echo "Executable installed as ~/bincpp/$ARCH/pcl_read"
