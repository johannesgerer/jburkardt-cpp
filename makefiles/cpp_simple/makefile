#
#  On Chili and Kachina, this is where the C++ compiler is:
#
CPP = /usr/bin/g++

cpp_simple_output.txt: cpp_simple
	./cpp_simple > cpp_simple_output.txt
cpp_simple: cpp_simple.o midpoint.o f.o
	$(CPP) -o cpp_simple cpp_simple.o midpoint.o f.o -lm
cpp_simple.o: cpp_simple.cpp cpp_simple.hpp
	$(CPP) -c cpp_simple.cpp
midpoint.o: midpoint.cpp cpp_simple.hpp
	$(CPP) -c midpoint.cpp
f.o: f.cpp cpp_simple.hpp
	$(CPP) -c f.cpp
clean:
	/bin/rm cpp_simple.out cpp_simple cpp_simple.o midpoint.o f.o



