CXX=g++
CXXFLAGS= -O2

LIBS=

all: main

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: main.o PhaseSolver.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o main *~
