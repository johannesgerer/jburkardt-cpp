#!/bin/bash
#
g++ -c delaunay_tree_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling delaunay_tree_2d.cpp"
  exit
fi
rm compiler.txt
#
g++ delaunay_tree_2d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading delaunay_tree_2d.o."
  exit
fi
#
rm delaunay_tree_2d.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/delaunay_tree_2d
#
echo "A new version of delaunay_tree_2d has been created."
