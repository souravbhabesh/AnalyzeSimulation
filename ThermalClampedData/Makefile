SRC = $(wildcard *.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = variables.h
CXX = g++
DEBUG = -g
CXXFLAGS = -Wall -c $(DEBUG) -std=c++11
LFLAGS = $(DEBUG) -O2 -Wall 

GSDCOM : $(OBJS)
	$(CXX) -o GSDCOM $(OBJS) $(LFLAGS)

gsd.o : $(DEPS) gsd.h  gsd_fn.cpp
	$(CXX) $(CXXFLAGS) gsd_fn.cpp

analyze.o : $(DEPS) analyze.h  analyze.cpp 
	$(CXX) $(CXXFLAGS) analyze.cpp

gsd_read.o : $(DEPS) gsd.h gsd_read.h  gsd_read.cpp
	$(CXX) $(CXXFLAGS) gsd_read.cpp

main.o : $(DEPS) gsd.h gsd_read.h analyze.h main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

clean:
	\rm *.o *~ GSDCOM

