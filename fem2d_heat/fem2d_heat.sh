#!/bin/bash
#
g++ -c -g fem2d_heat.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_heat.cpp"
  exit
fi
rm compiler.txt
#
mv fem2d_heat.o ~/libcpp/$ARCH
#
echo "Object code installed as $HOME/libcpp/$ARCH/fem2d_heat.o"

