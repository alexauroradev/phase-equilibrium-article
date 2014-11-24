CXX=g++
CXXFLAGS= -O0 -g

LIBS=

all: main

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: Phases.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o main *~ res.csv
