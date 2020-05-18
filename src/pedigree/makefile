CC=g++ -std=gnu++11
CFLAGS=-I.
INC=-I. -I../general
LIB=-L../general

OBJECTS = mainR.o analysis_fam.o Loc.o Markerscore.o \
           crosses.o popt.o read_ldf.o read_map.o

all: pedigree

pedigree: $(OBJECTS)$
	$(CC) $(INC) $(LIB) -o pedigree $(OBJECTS)$  -lgeneral

Loc.o: Loc.cpp
	$(CC) $(INC) -c Loc.cpp

analysis_fam.o: analysis_fam.cpp
	$(CC) $(INC) -c analysis_fam.cpp

crosses.o: crosses.cpp
	$(CC) $(INC) -c crosses.cpp

mainR.o: mainR.cpp
	$(CC) $(INC) -c mainR.cpp

markerscore.o: markerscore.cpp
	$(CC) $(INC) -c markerscore.cpp

popt.o: popt.cpp
	$(CC) $(INC) -c popt.cpp

read_ldf.o: read_ldf.cpp
	$(CC) $(INC) -c read_ldf.cpp

read_map.o: read_map.cpp
	$(CC) $(INC) -c read_map.cpp

clean:
	del $(OBJECTS)$


