#!/bin/bash
#
g++ -c -I/$HOME/include hypersphere_properties_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hypersphere_properties_prb.cpp"
  exit
fi
#
g++ hypersphere_properties_prb.o /$HOME/libcpp/$ARCH/hypersphere_properties.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hypersphere_properties_prb.o."
  exit
fi
#
rm hypersphere_properties_prb.o
#
mv a.out hypersphere_properties_prb
./hypersphere_properties_prb > hypersphere_properties_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hypersphere_properties_prb."
  exit
fi
rm hypersphere_properties_prb
#
echo "Program output written to hypersphere_properties_prb_output.txt"
