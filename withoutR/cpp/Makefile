

# which compiler to use 
#CC = gcc -O3 -Wall -Wunused 

#if for relase
#CC = g++ -O3
#if for debug
CC = g++ -O3 -g -Wall -Wunused -DDEBUG_R


#############################################################
################## TOP LEVEL TARGETS       ##################
#############################################################

all: test_iso

clean:
	rm -f *.o \
	rm -f test_iso 

test_iso: test_iso.o riso.o
	$(CC) test_iso.o riso.o -o test_iso
	
#############################################################
################## C++ OBJECTS	           ##################
#############################################################

test_iso.o: test_iso.cpp riso.h
	$(CC) -c test_iso.cpp
	
	
	
#############################################################
################## C++ OBJECTS in R src    ##################
################## symbolic links!!        ##################
################## ln -s /path/to/file /path/to/symlink #####
#############################################################


riso.o: riso.cpp riso.h
	$(CC) -c riso.cpp -o riso.o
	

