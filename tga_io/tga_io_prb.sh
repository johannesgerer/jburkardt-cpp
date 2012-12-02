#!/bin/bash
#
g++ -c -g -I/$HOME/include tga_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tga_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ tga_io_prb.o /$HOME/libcpp/$ARCH/tga_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tga_io_prb.o."
  exit
fi
#
rm tga_io_prb.o
#
mv a.out tga_io_prb
./tga_io_prb > tga_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tga_io_prb."
  exit
fi
rm tga_io_prb
#
echo "Program output written to tga_io_prb_output.txt"
