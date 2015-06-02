#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include z_sample.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling z_sample.cpp"
  exit
fi
#
#  Link and load
#
#g++ z_sample.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran z_sample.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ z_sample.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading z_sample.o"
  exit
fi
rm z_sample.o
mv a.out z_sample
#
#  Run
#
./z_sample > z_sample_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running z_sample"
  exit
fi
rm z_sample
#
#  Terminate.
#
echo "Program output written to z_sample_output.txt"
