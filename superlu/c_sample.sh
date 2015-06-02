#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include c_sample.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling c_sample.cpp"
  exit
fi
#
#  Link and load
#
#g++ c_sample.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran c_sample.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ c_sample.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c_sample.o"
  exit
fi
rm c_sample.o
mv a.out c_sample
#
#  Run
#
./c_sample > c_sample_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c_sample"
  exit
fi
rm c_sample
#
#  Terminate.
#
echo "Program output written to c_sample_output.txt"
