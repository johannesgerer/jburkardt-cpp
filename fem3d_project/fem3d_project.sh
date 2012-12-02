#!/bin/bash
#
g++ -c -g fem3d_project.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem3d_project.cpp"
  exit
fi
rm compiler.txt
#
g++ fem3d_project.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem3d_project.o."
  exit
fi
#
rm fem3d_project.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem3d_project
#
echo "Executable installed as ~/bincpp/$ARCH/fem3d_project"
