#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include d_sample_hb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling d_sample_hb.cpp"
  exit
fi
#
#  Link and load
#
#g++ d_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran d_sample_hb.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ d_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading d_sample_hb.o"
  exit
fi
rm d_sample_hb.o
mv a.out d_sample_hb
#
#  Run
#
./d_sample_hb < sample_rua.txt sample_hb> d_sample_hb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running d_sample_hb"
  exit
fi
rm d_sample_hb
#
#  Terminate.
#
echo "Program output written to d_sample_hb_output.txt"
