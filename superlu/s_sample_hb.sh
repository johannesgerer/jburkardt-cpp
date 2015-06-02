#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include s_sample_hb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling s_sample_hb.cpp"
  exit
fi
#
#  Link and load
#
#g++ s_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran s_sample_hb.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ s_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading s_sample_hb.o"
  exit
fi
rm s_sample_hb.o
mv a.out s_sample_hb
#
#  Run
#
./s_sample_hb < sample_rua.txt > s_sample_hb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running s_sample_hb"
  exit
fi
rm s_sample_hb
#
#  Terminate.
#
echo "Program output written to s_sample_hb_output.txt"
