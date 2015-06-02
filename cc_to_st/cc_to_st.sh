#!/bin/bash
#
cp cc_to_st.hpp /$HOME/include
#
g++ -c -I /$HOME/include cc_to_st.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cc_to_st.cpp"
  exit
fi
#
mv cc_to_st.o ~/libcpp/$ARCH/cc_to_st.o
#
echo "Library installed as ~/libcpp/$ARCH/cc_to_st.o"
