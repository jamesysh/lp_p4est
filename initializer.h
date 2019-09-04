#ifndef __INITIALIZER_H__
#define __INITIALIZER_H__


class Initializer {

    public:
    
        Initializer();
        ~Initializer();




        int initlevel ; //init level of octree
        int maxlevel ;
        int elem_particles; //max number of particles per octant
        double cfl_coefficient;
        double endt;
        double localspacing;

};

#endif // __INITIALIZER_H__
