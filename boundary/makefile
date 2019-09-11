CC     = mpic++
DEBUG  = -g
#OMP =
LAPACK_DIR=/gpfs/home/shyyuan/local/lapack-3.8.0
MAIN_DIR:=${CURDIR}/../
BOUNDARY_DIR=$(MAIN_DIR)/boundary/
STATE_DIR=$(MAIN_DIR)/state/
GEOMETRY_DIR=$(MAIN_DIR)/geometry/
P4EST_INC = /home/syuan/p4est/include
INCS   = -I $(BOUNDARY_DIR) -I $(STATE_DIR) -I $(GEOMETRY_DIR) -I $(MAIN_DIR) -I $(P4EST_INC)
LIBS   = -L $(LAPACK_DIR) 
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP) $(INCS)
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS) $(OMP)
OBJS = boundary.o boundary_pellet.o 
all: $(OBJS)

boundary.o: boundary.h boundary.cpp
	${CC} ${CFLAGS} boundary.cpp
boundary_pellet.o: boundary.h boundary_pellet.h boundary_pellet.cpp
	${CC} ${CFLAGS} boundary_pellet.cpp




clean:
	rm *.o


