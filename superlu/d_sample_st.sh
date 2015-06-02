#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include d_sample_st.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling d_sample_st.cpp"
  exit
fi
#
#  Link and load
#
#g++ d_sample_st.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran d_sample_st.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ d_sample_st.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading d_sample_st.o"
  exit
fi
rm d_sample_st.o
mv a.out d_sample_st
#
#  Run
#
./d_sample_st < sample_rst.txt > d_sample_st_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running d_sample_st"
  exit
fi
rm d_sample_st
#
#  Terminate.
#
echo "Program output written to d_sample_st_output.txt"
