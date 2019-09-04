CC = mpic++ 
DEBUG = -g
LAPACK_DIR = /home/syuan/local/lapack
P4EST_INC = /home/syuan/p4est/include
P4EST_LIB = /home/syuan/p4est/lib
INCS = -I $(P4EST_INC)
LIBS = -L $(P4EST_LIB)
CFLAGS = -Wall -std=c++11 -c  $(DEBUG) $(INCS) $(LIBS) 
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS)

OBJS = lp_main.o particle_data.o initializer.o octree_manager.o
all: $(OBJS) lp

lp_main.o: lp_main.cpp
	$(CC) $(CFLAGS) lp_main.cpp
particle_data.o: particle_data.cpp particle_data.h
	$(CC) $(CFLAGS) particle_data.cpp
initializer.o: initializer.cpp initializer.h
	$(CC) $(CFLAGS) initializer.cpp
octree_manager.o: octree_manager.cpp octree_manager.h
	$(CC) $(CFLAGS) octree_manager.cpp

lp: $(OBJS)  
	$(CC) $(LFLAGS) $(OBJS) -o lp -lsc -lp4est 

clean:
	rm *.o
	rm lp
