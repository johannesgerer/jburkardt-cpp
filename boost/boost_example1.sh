#!/bin/bash
#
g++ -c -I /opt/local/include boost_example1.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling boost_example1.cpp"
  exit
fi
rm compiler.txt
#
g++ boost_example1.o -L/usr/local/lib -lm
if [ $? -ne 0 ]; then
  echo "Errors while loading boost_example1.o"
  exit
fi
rm boost_example1.o
#
mv a.out ~/bincpp/$ARCH/boost_example1
#
echo "Executable installed as ~/bincpp/$ARCH/boost_example1"
