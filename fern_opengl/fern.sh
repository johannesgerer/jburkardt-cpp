#!/bin/bash
#
g++ -c fern.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fern.cpp"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ fern.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
# g++ fern.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fern.o"
  exit
fi
#
rm fern.o
mv a.out ~/bincpp/$ARCH/fern
#
echo "Executable installed as ~/bincpp/$ARCH/fern"
