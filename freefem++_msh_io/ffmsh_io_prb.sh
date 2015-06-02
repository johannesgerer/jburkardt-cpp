#!/bin/bash
#
g++ -c -I/$HOME/include ffmsh_io_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ffmsh_io_prb.c"
  exit
fi
#
g++ -o ffmsh_io_prb ffmsh_io_prb.o /$HOME/libcpp/$ARCH/ffmsh_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ffmsh_io_prb.o."
  exit
fi
#
rm ffmsh_io_prb.o
#
./ffmsh_io_prb > ffmsh_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ffmsh_io_prb."
  exit
fi
rm ffmsh_io_prb
#
echo "Program output written to ffmsh_io_prb_output.txt"
