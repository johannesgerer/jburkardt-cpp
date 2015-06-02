#!/bin/bash
#
g++ -c -I/$HOME/include cc_to_st_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cc_to_st_prb.cpp"
  exit
fi
#
g++ cc_to_st_prb.o /$HOME/libcpp/$ARCH/cc_to_st.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cc_to_st_prb.o"
  exit
fi
#
rm cc_to_st_prb.o
#
mv a.out cc_to_st_prb
./cc_to_st_prb > cc_to_st_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cc_to_st_prb."
  exit
fi
rm cc_to_st_prb
#
echo "Program output written to cc_to_st_prb_output.txt"
