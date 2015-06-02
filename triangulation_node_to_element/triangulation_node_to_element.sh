#!/bin/bash
#
g++ -c -g triangulation_node_to_element.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_node_to_element.cpp"
  exit
fi
rm compiler.txt
#
g++ triangulation_node_to_element.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_node_to_element.o."
  exit
fi
#
rm triangulation_node_to_element.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangulation_node_to_element
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_node_to_element"
