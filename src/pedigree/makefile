CC=g++ -std=gnu++11
CFLAGS=-I.
INC=-I. -I../lib/general
LIB=-L../lib/general

OBJECTS = mainR.o analysis_fam.o Loc.o Markerscore.o analysis_sub.o \
           crosses.o IndProp.o pedigree.o popt.o read_ldf.o read_loc.o read_map.o

#read_ped.o analysis_ped.o output.o output_bin.o output_flx.o output_tab.o

all: pedigree

pedigree: $(OBJECTS)$
	$(CC) $(INC) $(LIB) -o pedigree $(OBJECTS)$  -lgeneral

Loc.o: Loc.cpp
	$(CC) $(INC) -c Loc.cpp

analysis_fam.o: analysis_fam.cpp
	$(CC) $(INC) -c analysis_fam.cpp

#analysis_ped.o: analysis_ped.cpp
#	$(CC) $(INC) -c analysis_ped.cpp

analysis_sub.o: analysis_sub.cpp
	$(CC) $(INC) -c analysis_sub.cpp

crosses.o: crosses.cpp
	$(CC) $(INC) -c crosses.cpp

IndProp.o: IndProp.cpp
	$(CC) $(INC) -c IndProp.cpp

mainR.o: mainR.cpp
	$(CC) $(INC) -c mainR.cpp

markerscore.o: markerscore.cpp
	$(CC) $(INC) -c markerscore.cpp

#output.o: output.cpp
	#$(CC) $(INC) -c output.cpp

#output_bin.o: output_bin.cpp
#	$(CC) $(INC) -c output_bin.cpp

#output_flx.o: output_flx.cpp
#	$(CC) $(INC) -c output_flx.cpp

#output_tab.o: output_tab.cpp
#	$(CC) $(INC) -c output_tab.cpp

pedigree.o: pedigree.cpp
	$(CC) $(INC) -c pedigree.cpp

popt.o: popt.cpp
	$(CC) $(INC) -c popt.cpp

read_ldf.o: read_ldf.cpp
	$(CC) $(INC) -c read_ldf.cpp

read_loc.o: read_loc.cpp
	$(CC) $(INC) -c read_loc.cpp

read_map.o: read_map.cpp
	$(CC) $(INC) -c read_map.cpp

#read_ped.o: read_ped.cpp
	#$(CC) $(INC) -c read_ped.cpp

clean:
	del $(OBJECTS)$


