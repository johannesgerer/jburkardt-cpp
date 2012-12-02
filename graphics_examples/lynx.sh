#!/bin/bash
#
g++ -c -I/usr/local/dislin lynx.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lynx.cpp"
  exit
fi
rm compiler.txt
#
g++ lynx.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lynx.o."
  exit
fi
#
rm lynx.o
#
mv a.out lynx
./lynx > lynx_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lynx."
  exit
fi
rm lynx
#
echo "Program output written to lynx_output.txt"
