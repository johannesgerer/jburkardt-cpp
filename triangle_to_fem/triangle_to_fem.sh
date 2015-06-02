#!/bin/bash
#
g++ -c triangle_to_fem.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_to_fem.cpp"
  exit
fi
#
g++ triangle_to_fem.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_to_fem.o."
  exit
fi
rm triangle_to_fem.o
#
mv a.out ~/bincpp/$ARCH/triangle_to_fem
#
echo "Executable installed as ~/bincpp/$ARCH/triangle_to_fem"
