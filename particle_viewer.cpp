#include <iostream>
#include "particle_viewer.h"
#include <string.h>
#include <cassert>

using namespace std;
ParticleViewer::ParticleViewer(Global_Data *g, const std::string &filename,int numd ){
    gdata = g;
    outputfilename = filename;
    numdigit = 4;
    writestep = 0;
}

string ParticleViewer::rightFlush(size_t numDigits) {
	
	assert(pow(10,numDigits) > writestep);

	string result;

	if(writestep == 0) numDigits--;
	for(size_t i=writestep; i!=0; i /= 10) numDigits--;

	for( ; numDigits>0; numDigits--) result.push_back('0');
	
	result += to_string(writestep); 
	
	return result;

}

void ParticleViewer:: writeResult(double t){
    
   static bool FIRST = true;
   size_t lpnum = gdata->particle_data->elem_count;
   size_t li;
   int mpirank = gdata->mpirank;
   string filename_t = outputfilename + rightFlush(numdigit)+"_";
   string filename = filename_t + to_string(mpirank) + ".vtk";
    pdata_t *pad;
    
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Unable to open file: %s\n",filename.c_str()); 
		return;
	}
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",0.0);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",(long int)lpnum);
    
    pad = (pdata_t *)gdata->particle_data->array;
    //pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){
    fprintf(outfile,"%.16g %.16g %.16g\n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
    pad++ ;
    }
  
    fprintf(outfile,"POINT_DATA %ld\n",(long int)lpnum);

	fprintf(outfile,"VECTORS Velocity double\n");
    
    pad = (pdata_t *)gdata->particle_data->array;
    //pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){
    fprintf(outfile,"%.16g %.16g %.16g\n",pad->v[0],pad->v[1],pad->v[2]);
    pad++ ;
    }

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    pad = (pdata_t *)gdata->particle_data->array;
    for(li = 0; li<lpnum; li++){
    fprintf(outfile,"%.16g\n",pad->pressure);
    pad++ ;
    }
    
    fclose(outfile);
    if(mpirank == 0){
    FILE *visitfile;
    string fname = "output_data.visit";
    visitfile = fopen(fname.c_str(),"a");
   
	if(visitfile==nullptr) {
		printf("Unable to open file: %s\n",fname.c_str()); 
		return;
	}
    if(FIRST){
        FIRST = false;
        fprintf(visitfile,"!NBLOCKS %d\n",gdata->mpisize);
    
        }
        for(int i=0;i<gdata->mpisize;i++){
            string name  = filename_t +to_string(i) + ".vtk\n";
            fprintf(visitfile,name.c_str());
            
            }
    
    fclose(visitfile);
    }   
    writestep ++;

}
