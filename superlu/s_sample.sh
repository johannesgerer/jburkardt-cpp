#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include s_sample.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling s_sample.cpp"
  exit
fi
#
#  Link and load
#
#g++ s_sample.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran s_sample.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ s_sample.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading s_sample.o"
  exit
fi
rm s_sample.o
mv a.out s_sample
#
#  Run
#
./s_sample > s_sample_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running s_sample"
  exit
fi
rm s_sample
#
#  Terminate.
#
echo "Program output written to s_sample_output.txt"
