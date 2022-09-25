#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Geometry>
#include <chrono>
#include "mpi.h"
#include <vector>
#include <random>
#include <sstream>

#include "General_Utilities.hpp"
#include "Physical_Parameters.hpp"
#include "IC_Generator.hpp"
#include "Trajectory_Simulation.hpp"
#include "PREM.hpp"
#include "RN_Generators.hpp"

using namespace std;
using namespace std::chrono;

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

	std::vector<double> linspaced;

	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	double num = static_cast<double>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1) 
    {
		linspaced.push_back(start);
		return linspaced;
    }

	double delta = (end - start) / (num - 1);

	for(int i=0; i < num-1; ++i)
    {
		linspaced.push_back(start + delta * i);
    }
	linspaced.push_back(end); // I want to ensure that start and end
	// are exactly the same as the input
	return linspaced;
}

int main(int argc, char *argv[])
{

//INITIALIZATION	
////////////////////////////////////////////////////////////
	//MPI Enviroment
		// Initialize the MPI environment
	    	MPI_Init(NULL, NULL);
	    // Get the number of processes
	    	int numprocs;
	    	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	    // Get the ID number of the process
	    	int myRank;
	    	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

		//Starting time
	    high_resolution_clock::time_point tStart,t1,t2,tEnd;
	    double durationMain=0.0;

	    if(myRank==0)
	    {
	    	tStart = high_resolution_clock::now();
	    }
	//Read in configuration file, set up the detector and define the earth velocity
		Read_Config_File(argv[1]);
	//Initialize Logfile and create folders for inputfiles
		//Copy_Config_File(argv[1]);  if (argc>=2) {mass_dm[0]=atof(argv[1]);printf("DM mass: %le\n",mass_dm[0]);}
	  if (argc>=3) {mChi=atof(argv[2])*MeV;}
	  if (argc>=4) {sigma0=atof(argv[3])*pb;}
	  if (argc>=5) {Global_SampleSize_Initial=atol(argv[4]);}
	  //if (argc>=6) {strcpy(outfname,argv[5]);printf("Output file name: %s\n", outfname);}
		
		if(myRank==0)
		{
			cout << "sample size: " << Global_SampleSize_Initial << endl;
			cout << "DM mass: " << mChi/MeV << " MeV" << endl;
			cout << "Cross section: " << sigma0/pb << " pb" << endl;
			cout << "Number of processors: " << numprocs << endl;
			std::chrono::time_point<std::chrono::system_clock> start;
			start = std::chrono::system_clock::now();
			std::time_t start_time = std::chrono::system_clock::to_time_t(start);
			cout <<"\n##############################"<<endl
			<<"DaMaSCUS"<<version<<" - Simulation"<<endl<<endl
			<<"Starting Time: " <<std::ctime(&start_time)
			<<"Simulation ID: " <<SimID <<endl<<endl;
		}
		//<<"Creating logfile." <<endl;
		//LogFile_Initialize(numprocs);

	//Initialize Random Number Generator:
		std::random_device rd;
		std::mt19937 PRNG(rd());
	//Synchronize processes:
		MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////
	
	//1. Initial MC run without scatterings:
		//Initialize the earth model for initial run
		//Deactivate formfactor, if it is used:

			//computing time
			//t1 = high_resolution_clock::now();
			//durationIni =1e-6*duration_cast<microseconds>( t1 - tStart ).count();
			//cout <<"\tInitial run finished\t(" <<floor(durationIni) <<" s)." <<endl <<endl;

	//2. MC run with scatterings and velocity data collection.
		//Initialize the earth model for the main run
			Initialize_PREM(mChi,sigma0);
		//if (myRank==0) cout <<"Start main MC simulation run with scatterings." <<endl;
		//Particle counter per isodetection ring
			
		//MPI Output File and Offset
			//MPI_Offset offset;
 			//MPI_File   file_Velocity[Isodetection_Rings];
 			//MPI_File   file_Weights[Isodetection_Rings];
 		   	//MPI_Status status;
 		//Velocity Output files including process dependent offset
			
		
		unsigned long long int Local_SampleSize_Initial =ceil((double)Global_SampleSize_Initial/numprocs);
		Global_SampleSize_Initial=Local_SampleSize_Initial*numprocs;

		int NRbins = 100;
		unsigned long long int Ncap_Global = 0;
		unsigned long long int Ncap_Local = 0;
		unsigned long long int Nstats_Global[NRbins];
		unsigned long long int Nstats_Local[NRbins];
		for(int i=0;i<NRbins;i++)
		{
			Nstats_Global[i] = 0;
			Nstats_Local[i] = 0;
		}
		//shells by radius
		vector<double> Rs = linspace(0., 1., NRbins+1);

		//Global_SampleSize_Initial
		for(unsigned long long int i=0; i<Local_SampleSize_Initial; i++)
		{
			Event IC = InitialCondition(0,vEarth,PRNG);
			Eigen::Vector3d vini=IC.Velocity();
			//reject initial velocity smaller than the escape velocity
			while (vini.norm()<=vescape(rEarthatm))
			{
				IC = InitialCondition(0,vEarth,PRNG);
				vini=IC.Velocity();
			}
			Trajectory trajectory=ParticleTrack(mChi,sigma0,IC,vcut,PRNG); 	
			Event Elast=trajectory.TrajectoryEnd();
			Eigen::Vector3d x=Elast.Position();
			Eigen::Vector3d v=Elast.Velocity();
			//cout << x << endl;
			//cout << v << endl;
			//cout << x.norm()/rEarth << " " << v.norm()/km*sec << endl;
			//if (trajectory.Trajectory_Type()>2)
			//{
			//	cout << "Warning: particle does not enter the Earth, exit." << endl;
			//	exit(0);
			//}
			double xnorm = x.norm();
			if (xnorm < rEarthatm)
			{
				int idx = floor(x.norm()/rEarthatm/(1./NRbins));
				//cout << idx << endl;
				++Nstats_Local[idx];
				if (v.norm()<vescape(rEarthatm)) ++Ncap_Local;
			}
			
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&Nstats_Local,&Nstats_Global,NRbins,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&Ncap_Local,&Ncap_Global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		
		if(myRank==0)
		{
			//computing time
			t1 = high_resolution_clock::now();
			durationMain =1e-6*duration_cast<microseconds>( t1 - tStart ).count();
			cout <<"\tMain run finished\t(" <<floor(durationMain) <<" s)." <<endl <<endl;
		}		
		
		if(myRank==0)
		{
			long double fcap = (long double) Ncap_Global / (long double) Global_SampleSize_Initial;
			cout << "capture fraction: " << fcap << endl;
			std::ostringstream streamObj1, streamObj2;
			streamObj1 << mChi/GeV;
			streamObj2 << sigma0/pow(cm,2);
			
			/*
			//write output file			
			string fname = "../results_atm_frac/DM_" + streamObj1.str() + "GeV" + ".txt";
			cout << "Writing to file " << fname << "..." << endl;
			FILE *fp = fopen(fname.c_str(),"a+");
			fprintf(fp, "%1.6e %1.6Le\n", sigma0/pow(cm,2), fcap);			
			fclose(fp);
			cout << "Done." << endl;
			*/
			
		}
		
		MPI_Finalize();
    	return 0;
////////////////////////////////////////////////////////////
}
