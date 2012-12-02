CXX=g++
#CXXFLAGS= -O3
CXXFLAGS= -O3 -m64
NETCDF=/opt/netcdf/3.6.3/64

DBGFLAGS= -g

PLATFORM=_MACOS
PLATFORM=_LINUX

ifeq ($(PLATFORM),_LINUX)
	CXXLIBS = -I$(NETCDF)/include -L$(NETCDF)/lib -lGL -lglut -lnetcdf_c++ -lnetcdf -lGLU -lstdc++
endif

ifeq ($(PLATFORM),_MACOS)
	CXXLIBS = -L$(NETCDF)/lib -lnetcdf -lnetcdf_c++ -I$(NETCDF)/include -lstdc++ -framework OpenGL -framework GLUT
endif

all:
	$(CXX) mpas_draw.cpp vec_utils.cpp netcdf_utils.cpp $(CXXLIBS) $(CXXFLAGS) -D$(PLATFORM) -o mpas_draw

debug:
	$(CXX) mpas_draw.cpp vec_utils.cpp netcdf_utils.cpp $(CXXLIBS) $(DBGFLAGS) -D$(PLATFORM) -o mpas_draw

clean:
	rm mpas_draw

