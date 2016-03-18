EXTLIB = External/lib
EXTINC = External/include

CFLAGS += -I'External/include'
CFLAGS += -std=c++11
CFLAGS += -Wall
LDFLAGS += -L'External/lib'
LDLIBS += -lUtils -lSimplifier -lstdc++

all: bin/main

#.PHONY: clean

bin/main: obj/main.o $(EXTLIB)/libUtils.a $(EXTLIB)/libSimplifier.a
	g++ $(CFLAGS) obj/main.o -o bin/main $(LDFLAGS) $(LDLIBS)

obj/main.o:
	g++ $(CFLAGS) -c main.cpp -o obj/main.o

$(EXTLIB)/libUtils.a: obj/VectorMath.o obj/HydraExport.o
	ar rcs $(EXTLIB)/libUtils.a obj/VectorMath.o obj/HydraExport.o

obj/VectorMath.o:
	g++ $(CFLAGS) -c Utils/VectorMath.cpp -o obj/VectorMath.o

obj/HydraExport.o:
	g++ $(CFLAGS) -c Utils/HydraExport.cpp -o obj/HydraExport.o

$(EXTLIB)/libSimplifier.a: obj/Simplifier.o
	ar rcs $(EXTLIB)/libSimplifier.a obj/Simplifier.o

obj/Simplifier.o:
	g++ $(CFLAGS) -c Simplifier/main.cpp -o obj/Simplifier.o

clean:
	rm -r obj/*
	rm -r External/lib/*