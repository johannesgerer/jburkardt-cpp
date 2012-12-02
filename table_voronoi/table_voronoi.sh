#!/bin/bash
#
g++ -c table_voronoi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_voronoi.cpp"
  exit
fi
rm compiler.txt
#
g++ table_voronoi.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_voronoi.o."
  exit
fi
#
rm table_voronoi.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/table_voronoi
#
echo "Executable installed as ~/bincpp/$ARCH/table_voronoi"
