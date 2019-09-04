#include <iostream>

#include "sc.h"
#include "mpi.h"
#include "particle_data.h"

using namespace std;

Global_Data:: Global_Data(Initializer* init){

    initlevel = init->initlevel;
    maxlevel = init->maxlevel;
    elem_particles = init->elem_particles;
    
    


}


Global_Data:: ~Global_Data(){}
