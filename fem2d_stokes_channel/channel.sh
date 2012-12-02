#!/bin/bash
#
g++ -c -g channel.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling channel.cpp"
  exit
fi
rm compiler.txt
#
g++ ~/lib/$ARCH/free_fem_stokes.o channel.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading free_fem_stokes.o + channel.o"
  exit
fi
rm channel.o
#
chmod ugo+x a.out
mv a.out channel
./channel nodes6.txt triangles6.txt > channel_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running channel."
  exit
fi
rm channel
#
echo "The channel program has been executed."
