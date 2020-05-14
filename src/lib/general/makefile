CC=g++ -std=gnu++11
CFLAGS=-I.

OBJECTS = Args.o excepthnd.o HMMalgo.o matvec.o mblexcept.o misc.o TransMatSym2D.o \
          manipulators.o util_genetics.o dir.o  util.o linalg.o eigen.o Ordgeno.o InhVector.o 

all: general.lib

general.lib: $(OBJECTS)$  
	ar rcs general.lib $(OBJECTS)$

Args.o: Args.cpp
	$(CC) -c Args.cpp

excepthnd.o: excepthnd.cpp
	$(CC) -c excepthnd.cpp

HMMalgo.o: HMMalgo.cpp 
	$(CC) -c HMMalgo.cpp

matvec.o: matvec.cpp 
	$(CC) -c matvec.cpp

mblexcept.o: mblexcept.cpp 
	$(CC) -c mblexcept.cpp

misc.o: misc.cpp
	$(CC) -c misc.cpp

TransMatSym2D.o: TransMatSym2D.cpp
	$(CC) -c TransMatSym2D.cpp

manipulators.o: manipulators.cpp	
	$(CC) -c manipulators.cpp

util_genetics.o: util_genetics.cpp
	$(CC) -c util_genetics.cpp 

dir.o: dir.cpp
	$(CC) -c dir.cpp 

util.o: util.cpp
	$(CC) -c util.cpp

linalg.o: linalg.cpp
	$(CC) -c linalg.cpp

eigen.o: eigen.cpp
	$(CC) -c eigen.cpp

OrdGeno.o: OrdGeno.cpp
	$(CC) -c OrdGeno.cpp

InhVector.o: InhVector.cpp
	$(CC) -c InhVector.cpp


clean: 
	del $(OBJECTS)$

