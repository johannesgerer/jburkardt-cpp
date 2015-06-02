#!/bin/bash
#
g++ -c triangle_to_xml.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_to_xml.cpp"
  exit
fi
#
g++ triangle_to_xml.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_to_xml.o."
  exit
fi
rm triangle_to_xml.o
#
mv a.out ~/bincpp/$ARCH/triangle_to_xml
#
echo "Executable installed as ~/bincpp/$ARCH/triangle_to_xml"
