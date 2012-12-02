#!/bin/bash
#
g++ -c -I /opt/local/include boost_example2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling boost_example2.cpp"
  exit
fi
rm compiler.txt
#
g++ boost_example2.o -L/opt/local/lib -lboost_regex-mt
if [ $? -ne 0 ]; then
  echo "Errors while loading boost_example2.o"
  exit
fi
rm boost_example2.o
#
mv a.out ~/bincpp/$ARCH/boost_example2
#
echo "Executable installed as ~/bincpp/$ARCH/boost_example2"
