#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include c_sample_hb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling c_sample_hb.cpp"
  exit
fi
#
#  Link and load
#
#g++ c_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran c_sample_hb.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ c_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c_sample_hb.o"
  exit
fi
rm c_sample_hb.o
mv a.out c_sample_hb
#
#  Run
#
./c_sample_hb < sample_cua.txt > c_sample_hb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c_sample_hb"
  exit
fi
rm c_sample_hb
#
#  Terminate.
#
echo "Program output written to c_sample_hb_output.txt"
