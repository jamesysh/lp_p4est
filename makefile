CC = mpic++ 
DEBUG = -g
LAPACK_DIR = /home/syuan/local/lapack
P4EST_INC = /home/syuan/p4est/include
P4EST_LIB = /home/syuan/p4est/lib
MAIN_DIR:=${CURDIR}
BOUNDARY_DIR=$(MAIN_DIR)/boundary/
STATE_DIR=$(MAIN_DIR)/state/
GEOMETRY_DIR=$(MAIN_DIR)/geometry/

INCS = -I $(P4EST_INC) -I $(BOUNDARY_DIR) -I $(STATE_DIR) -I $(GEOMETRY_DIR) -I $(MAIN_DIR)
LIBS = -L $(P4EST_LIB)
CFLAGS = -Wall -std=c++11 -c  $(DEBUG) $(INCS) $(LIBS) 
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS)
vpath %.h $(GEOMETRY_DIR) $(STATE_DIR) $(BOUNDARY_DIR)

MAIN_OBJS = lp_main.o particle_data.o initializer.o octree_manager.o registrar.o lp_solver.o eos.o particle_viewer.o
GEOMETRY_OBJS = geometry.o geometry_pellet.o
STATE_OBJS = state.o state_pellet.o 

B_OBJS := $(foreach OBJ,$(BOUNDARY_OBJS),$(addprefix $(BOUNDARY_DIR),$(OBJ)))
S_OBJS := $(foreach OBJ,$(STATE_OBJS),$(addprefix $(STATE_DIR),$(OBJ)))
G_OBJS := $(foreach OBJ,$(GEOMETRY_OBJS),$(addprefix $(GEOMETRY_DIR),$(OBJ)))

OBJS = $(B_OBJS) $(S_OBJS) $(G_OBJS) $(MAIN_OBJS)

all: $(OBJS) lp

G_OBJS:
	cd $(GEOMETRY_DIR)&&make;
S_OBJS:
	cd $(STATE_DIR)&&make;



lp_main.o: lp_main.cpp initializer.h particle_data.h
	$(CC) $(CFLAGS) lp_main.cpp
particle_data.o: particle_data.cpp particle_data.h geometry.h initializer.h
	$(CC) $(CFLAGS) particle_data.cpp
initializer.o: initializer.cpp initializer.h
	$(CC) $(CFLAGS) initializer.cpp
octree_manager.o: octree_manager.cpp octree_manager.h particle_data.h
	$(CC) $(CFLAGS) octree_manager.cpp
registrar.o: registrar.h registrar.cpp
	$(CC) $(CFLAGS) registrar.cpp
lp_solver.o: lp_solver.h lp_solver.cpp initializer.h particle_data.h
	$(CC) $(CFLAGS) lp_solver.cpp
eos.o: eos.h eos.cpp
	$(CC) $(CFLAGS) eos.cpp

particle_viewer.o: particle_data.h particle_viewer.cpp particle_viewer.h
	$(CC) $(CFLAGS) particle_viewer.cpp
lp: $(OBJS)  
	$(CC) $(LFLAGS) $(OBJS) -o lp -lsc -lp4est 


clean:
	
	cd $(GEOMETRY_DIR)&&make clean; 
	cd $(STATE_DIR)&&make clean;
	rm *.o
	rm lp
