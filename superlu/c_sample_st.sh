#!/bin/bash
#
#  Compile
#
g++ -c -I/$HOME/include c_sample_st.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling c_sample_st.cpp"
  exit
fi
#
#  Link and load
#
#g++ c_sample_st.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran c_sample_st.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas -lstdc++
#g++ c_sample_st.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c_sample_st.o"
  exit
fi
rm c_sample_st.o
mv a.out c_sample_st
#
#  Run
#
./c_sample_st < sample_cst.txt > c_sample_st_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c_sample_st"
  exit
fi
rm c_sample_st
#
#  Terminate.
#
echo "Program output written to c_sample_st_output.txt"
