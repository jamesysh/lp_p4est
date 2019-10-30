#include <iostream>
#include "particle_viewer.h"
#include <string.h>
#include <cassert>
#include "boundary_pellet.h"
using namespace std;
ParticleViewer::ParticleViewer(Global_Data *g, const std::string &filename,int numd ){
    gdata = g;
    outputfilename = filename;
    numdigit = numd;
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


void ParticleViewer:: writeGhost(int step){
    
   writestep = step; 
   static bool FIRST = true;
   size_t lpnum = gdata->particle_data->elem_count;
   size_t ghostn = gdata->lghostnum;
   size_t ghostnum;
   size_t li;
   int mpirank = gdata->mpirank;
   string filename_t = outputfilename +"/out" + rightFlush(numdigit)+"_";
   string filenameinvisitfile = "out"+rightFlush(numdigit) + "_";
   string filename = filename_t + to_string(mpirank) + ".vtk";
    pdata_t *pad;
    pdata_copy_t * padd;
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
	
	fprintf(outfile,"POINTS %ld double\n",(long int)ghostn);
    
    //pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){
       
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary)
            continue;
        ghostnum = pad->ghostneighbour->elem_count;

        padd = (pdata_copy_t *)pad->ghostneighbour->array; 
        for(size_t j=0;j<ghostnum;j++){
            if(gdata->dimension == 3)
                fprintf(outfile,"%.16g %.16g %.16g\n",padd->xyz[0],padd->xyz[1],padd->xyz[2]);
            else 
                fprintf(outfile,"%.16g %.16g %.16g\n",padd->xyz[0],padd->xyz[1],0);
         padd++;
        }
    }
    
    fprintf(outfile,"POINT_DATA %ld\n",(long int)ghostn);

	fprintf(outfile,"VECTORS Velocity double\n");
    
    //pad = (pdata_t *)sc_array_index_begin(particle_data);
   
    for(li = 0; li<lpnum; li++){
        
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary)
            continue;
        ghostnum = pad->ghostneighbour->elem_count;
        padd = (pdata_copy_t *)pad->ghostneighbour->array; 
        for(size_t j=0; j<ghostnum;j++){
         fprintf(outfile,"%.16g %.16g %.16g\n",padd->v[0],padd->v[1],padd->v[2]);
         padd++;
        }
    }
 

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    
    for(li = 0; li<lpnum; li++){
        
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary)
            continue;
        ghostnum = pad->ghostneighbour->elem_count;
        padd = (pdata_copy_t *)pad->ghostneighbour->array; 
        for(size_t j=0;j<ghostnum;j++){
         fprintf(outfile,"%.16g\n",padd->pressure);
         padd++;
        }
    }
 
    fclose(outfile);
    if(mpirank == 0){
    FILE *visitfile;
    string fname = outputfilename+ "output_alldata.visit";
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
            string name  = filenameinvisitfile + to_string(i) + ".vtk\n";
            fprintf(visitfile,name.c_str());
            
            }
    
    fclose(visitfile);
    }   

}
void ParticleViewer:: writeResult(int step, double time){
   writestep = step; 
   static bool FIRST = true;
   size_t lpnum = gdata->particle_data->elem_count;
   size_t li;
   size_t lfluidnum = gdata->lfluidnum;
   int mpirank = gdata->mpirank;
   string filename_t = outputfilename +"/out_fluid" + rightFlush(numdigit)+"_";
   string filenameinvisitfile = "out_fluid"+rightFlush(numdigit) + "_";
   string filename = filename_t + to_string(mpirank) + ".vtk";
    pdata_t *pad;
    
	FILE *outfile;

    //pad = (pdata_t *)sc_array_index_begin(particle_data);
    
    outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Unable to open file: %s\n",filename.c_str()); 
		return;
	}
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",0.0);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",(long int)lfluidnum);
    
    //pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        if(gdata->dimension == 3)
            fprintf(outfile,"%.16g %.16g %.16g\n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
        else 
            fprintf(outfile,"%.16g %.16g %.16g\n",pad->xyz[0],pad->xyz[1],0);
    }
  
    fprintf(outfile,"POINT_DATA %ld\n",(long int)lfluidnum);

	fprintf(outfile,"VECTORS Velocity double\n");
    
    //pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){

        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        if(gdata->dimension == 3)
            fprintf(outfile,"%.16g %.16g %.16g\n",pad->v[0],pad->v[1],pad->v[2]);
        else
            fprintf(outfile,"%.16g %.16g %.16g\n",pad->v[0],pad->v[1],0);
    }

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",(double)pad->pressure);
    }
    
	fprintf(outfile,"SCALARS density double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",1./(double)pad->volume);
    }
    
	fprintf(outfile,"SCALARS soundspeed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",pad->soundspeed);
    }
	
    fprintf(outfile,"SCALARS localspacing double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",pad->localspacing);
    }
    
    fprintf(outfile,"SCALARS ifhasghostnei int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%d\n",pad->ifhasghostneighbour);
    }
	fprintf(outfile,"SCALARS leftintegral double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",(double)pad->leftintegral);
    }
    
	fprintf(outfile,"SCALARS rightintegral double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    
	fprintf(outfile,"SCALARS deltaq double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",(double)pad->deltaq);
    }
    
	fprintf(outfile,"SCALARS qplusminus double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",(double)pad->qplusminus);
    }
    
    for(li = 0; li<lpnum; li++){
    
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        
        fprintf(outfile,"%.16g\n",(double)pad->rightintegral);
    }
    fclose(outfile);
    if(mpirank == 0){
    FILE *visitfile;
    string fname = outputfilename +"/output_fluiddata.visit";
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
            string name  = filenameinvisitfile + to_string(i) + ".vtk\n";
            fprintf(visitfile,name.c_str());
            
            }
    
    fclose(visitfile);
    }   
    
    if( gdata->pelletnumber){
        size_t li, lnump = gdata->particle_data->elem_count;
        
        double dx = gdata->initlocalspacing;
        double mfr = 0;
        double r = 1;
        double dr = 3*dx;
        double tr, vr;
        double x, y, z;
        double vx, vy, vz;
        double mfr_g = 0;
        for(li = 0; li<lnump; li++){
            pad = (pdata_t *) sc_array_index(gdata->particle_data, li);
            if(pad->ifboundary)
                continue;
            x = pad->xyz[0];
            y = pad->xyz[1];
            z = pad->xyz[2];
            tr = (x*x) + (y*y) + (z*z);
            double dis = fabs(sqrt(tr) - r);
            if(dis < dr){
                vx = pad->v[0];
                vy = pad->v[1];
                vz = pad->v[2];
                vr = (vx*x+vy*y + vz*z)/sqrt(tr);
                mfr += pad->mass*vr;
                
                }
            
            }
        MPI_Allreduce(&mfr,&mfr_g,1,MPI_DOUBLE,MPI_SUM,gdata->mpicomm);

        mfr_g = mfr_g/2/dr;
        
        if(mpirank == 0){
        string mfrfilename = outputfilename + "/massflowrate.txt";
        FILE *mfroutfile;
        
        if(writestep>0)
                mfroutfile = fopen(mfrfilename.c_str(), "a");
        else{
                mfroutfile = fopen(mfrfilename.c_str(), "w");
		fprintf(mfroutfile,"Time  Massflowrate around pellets.\n");
	       }
        if(mfroutfile==nullptr) {
                printf("Unable to open file: %s\n",mfrfilename.c_str());
                return;
        }
        PelletInflowBoundary* b = (PelletInflowBoundary *)gdata->m_vBoundary[0]; 
        fprintf(mfroutfile,"%16g %16g %16g\n",time,b->massflowrate,mfr_g); 
        fclose(mfroutfile); 
        }
     }

}


void ParticleViewer::writeTXTFile(int step){

   string filename = outputfilename +"/out_fluid" + rightFlush(numdigit) + ".txt";
   
	MPI_File outfile;

    MPI_File_open(gdata->mpicomm,filename.c_str(),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&outfile);

   size_t lpnum = gdata->particle_data->elem_count;
   size_t li;
   double x,y,z;
   double vx,vy,vz;
   char buf[400];
   pdata_t *pad;
   for(li =0; li<lpnum; li++){
       pad = (pdata_t *)sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
       x = pad->xyz[0];
       y = pad->xyz[1];
       z = pad->xyz[2];
       vx = pad->v[0];
       vy = pad->v[1];
       vz = pad->v[2];
        
       snprintf(buf,400,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                         x, y, z, vx, vy, vz,pad->pressure,1./pad->volume,pad->soundspeed, pad->leftintegral, pad->deltaq);
       MPI_File_write_shared(outfile,buf,strlen(buf),MPI_CHAR,MPI_STATUS_IGNORE);

       }
    MPI_Barrier(gdata->mpicomm);
    MPI_File_close(&outfile);
    }










