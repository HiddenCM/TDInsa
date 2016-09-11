CC=/usr/bin/g++
#CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread  -fmax-errors=1
LDFLAGS=-pthread


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=readCorrection

all: $(EXEC)

readCorrection:   main.o td2_etudiants.o utils_etudiants.o
	$(CC) -o $@ $^ $(LDFLAGS)


utils_etudiants.o: utils_etudiants.cpp 
	$(CC) -o $@ -c $< $(CFLAGS)

td2_etudiants.o: td2_etudiants.cpp utils_etudiants.h
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: main.cpp td2_etudiants.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
