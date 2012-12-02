#!/bin/bash
#
g++ -c -g fem2d_project.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_project.cpp"
  exit
fi
rm compiler.txt
#
g++ fem2d_project.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_project.o."
  exit
fi
#
rm fem2d_project.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem2d_project
#
echo "Executable installed as ~/bincpp/$ARCH/fem2d_project"
