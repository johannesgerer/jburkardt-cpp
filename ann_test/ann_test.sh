#!/bin/bash
#
g++ -c -g -I/$HOME/include ann_test.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ann_test.cpp."
  exit
fi
rm compiler.txt
#
g++ -c -g -I/$HOME/include rand.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rand.cpp"
  exit
fi
rm compiler.txt
#
g++ ann_test.o rand.o -L/$HOME/libcpp/$ARCH -lann -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ann_test + rand + libann."
  exit
fi
#
rm ann_test.o
rm rand.o
#
mv a.out ~/bincpp/$ARCH/ann_test
#
echo "Executable installled as ~/bincpp/$ARCH/ann_test"
