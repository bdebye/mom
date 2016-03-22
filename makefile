
CC = g++
OPTION = -std=c++11 -fopenmp -O3
OBJECTS = antenna.o clock.o dataio.o feed.o \
			gmsh.o graphic.o main.o mathlib.o \
			rwg_edge.o type.o

INCLUDE = ./headers
VPATH = src

%.o: %.cpp
	$(CC) -c -o $@ $< $(OPTION) -I$(INCLUDE)

antenna_simu : $(OBJECTS)
	$(CC) -o $@ $^ $(OPTION)

.PHONY : clean 
clean : 
	rm antenna_simu $(OBJECTS)

