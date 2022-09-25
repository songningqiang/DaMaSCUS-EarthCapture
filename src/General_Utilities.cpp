#include "General_Utilities.hpp"

#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <libconfig.h++>

#include <sys/types.h> // required for stat.h
#include <sys/stat.h>	//required to create a folder

#include "Physical_Parameters.hpp"
#include "IC_Generator.hpp"
#include "Trajectory_Simulation.hpp"

using namespace libconfig;

//Input Parameter before reading in the config file. They must be defined at some point outside of a function.
	string SimID="default";
	string FormFactor="default";
	string experiment="default";
	string version = "v1.0";
	double mChi=0,sigma0=0;
	double Detector_Depth=0.0;
	double rhoDM=0.0;
	unsigned long long int Global_SampleSize_Initial = 0;
	unsigned int Global_SampleSize_Velocity=0;
	int Isodetection_Rings=0;

	int SimDate [3] = 	{0,0,0};
	int SimTime [3] = 	{0,0,0};
	double nJ2000 = FractionalDays(SimDate,SimTime);
	Eigen::Vector3d vEarth=EarthVelocity(nJ2000);
	double vcut=0.0;
	int SampleSize=0;

//Read config file and define Input Parameters. Also calculate fractional days, earth velocity, cutoff velocity and set up the detector.
	void Read_Config_File(const char *inputfile)
	{
		Config cfg;
		// Read the file. If there is an error, report it and exit.
			try
			{
			  cfg.readFile(inputfile);
			}
			catch(const FileIOException &fioex)
			{
				std::cerr << "I/O error while reading configuration file." << std::endl;
				exit(EXIT_FAILURE);
			}
			catch(const ParseException &pex)
			{
				std::cerr << "Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
				exit(EXIT_FAILURE);
			}
		//Simulation ID
			try
			{
				SimID = cfg.lookup("simID").c_str();
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'simID' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Experiment
			try
			{
				experiment = cfg.lookup("experiment").c_str();
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'simID' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Form Factor
			try
			{
				FormFactor = cfg.lookup("formfactor").c_str();
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'formfactor' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
			if(FormFactor!="None"&&FormFactor!="HelmApproximation")
			{
				cerr<<"Form factor option \"" <<FormFactor <<"\" is not a valid option." <<endl;
				exit(EXIT_FAILURE);
			}
		//Initial SampleSize
			try
			{
					Global_SampleSize_Initial = cfg.lookup("initialruns");
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'initialruns' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Main Samplesize
			try
			{
				Global_SampleSize_Velocity = cfg.lookup("samplesize");
				SampleSize = Global_SampleSize_Velocity;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'samplesize' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Isodetection rings
			try
			{
				Isodetection_Rings = cfg.lookup("rings");
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'rings' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Velocity cut off
			try
			{
				vcut = cfg.lookup("vcutoff");
				vcut*=cm/sec;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'vcutoff' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Mass
			try
			{
				mChi = cfg.lookup("mass");
				mChi*=MeV;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'mass' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Cross section
			try
			{
				sigma0 = cfg.lookup("sigma");
				sigma0*=pb;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'sigma' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//DM energy density
			try
			{
				rhoDM = cfg.lookup("rho");
				rhoDM*=GeV/cm/cm/cm;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'rho' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//Date
			try
			{
				SimDate[0] = cfg.lookup("date")[0];
				SimDate[1] = cfg.lookup("date")[1];
				SimDate[2] = cfg.lookup("date")[2];
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'date' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//Time
			try
			{
				SimTime[0] = cfg.lookup("time")[0];
				SimTime[1] = cfg.lookup("time")[1];
				SimTime[2] = cfg.lookup("time")[2];
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'time' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//Detector
			try
			{
				Detector_Depth = cfg.lookup("depth");
				Detector_Depth*=meter;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'time' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//Time and Detector
			nJ2000= FractionalDays(SimDate,SimTime);
		//Earth Velocity
			vEarth=EarthVelocity(nJ2000);
		
	}
	//Copy cfg file
	void Copy_Config_File(const char *inputfile)
	{
		std::ifstream    inFile;
		std::ofstream outFile;
		inFile.open(inputfile);
    	outFile.open("../data/"+SimID+".cfg");
		outFile << inFile.rdbuf();
		inFile.close();
		outFile.close();
	}

//For the simulator:
	//Create and initialize Logfile
		std::ofstream LogStream;
		void LogFile_Initialize(int numprocs,std::ofstream& f)
		{	
			std::chrono::time_point<std::chrono::system_clock> start;
			start = std::chrono::system_clock::now();
			std::time_t start_time = std::chrono::system_clock::to_time_t(start);
			// ofstream f;
			f.open("../results/"+SimID+".log");
			f	<<"//DaMaSCUS"<<version<<" Logfile" <<endl<<endl
				<<"//Simulation start: " <<std::ctime(&start_time) <<endl
				<<"//Simulation data:"<<endl
				<<"\tSimID\t\t" <<SimID <<endl
				<<"\tMPI processes\t" <<numprocs<<endl<<endl
				<<"//Number of initial runs and data points per isodetection ring:"<<endl
				<<"\tnInitial\t"<<(double)Global_SampleSize_Initial<<endl
				<<"\tData points\t" <<(double)Global_SampleSize_Velocity<<endl
				<<"\tIsodetection rings\t" <<(double)Isodetection_Rings<<endl <<endl
				<<"//Simulation time:" <<endl
				<<"\tDate\t\t" <<SimDate[0]<<"."<<SimDate[1]<<"."<<SimDate[2]<<endl
				<<"\tTime\t\t"<<SimTime[0]<<":"<<SimTime[1]<<":"<<SimTime[2]<<endl
				<<"\tnJ2000.0\t" <<nJ2000 <<endl
				<<"\tvEarth[km/s]\t" <<InUnits(vEarth.norm(),km/sec) <<endl <<endl
				<<"//DM Data:" <<endl
				<<"\tmDM[GeV]\t" <<mChi <<endl
				<<"\tsig0[pb]\t" <<sigma0/pb <<endl
				<<"\tForm Factor\t" <<FormFactor <<endl <<endl
				<<"//Velocity cut off:"<<endl
				<<"\tvcut[km/sec]\t"<<InUnits(vcut,km/sec)<<endl <<endl
				<<"//Detector depth:" <<endl
				<<"\td[m]\t\t" <<InUnits(Detector_Depth,meter) <<endl
				<<"//Expected depth crossing:\n \t\t\t" <<2.0*pow(rEarth-Detector_Depth,2.0)/pow(rEarth,2.0)<<endl;

			//Create folder for velocity data
				cout <<"Creating folder for velocity data." <<endl;
				std::string sPath = "../data/"+SimID+"_data";
				mode_t nMode = 0733; // UNIX style permissions
				int nError = 0;
				#if defined(_WIN32)
				  nError = _mkdir(sPath.c_str()); // can be used on Windows
				#else 
				  nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
				#endif
				if (nError != 0) 
				{
					cout <<"The folder already exists, data will be overwritten."<<endl<<endl;
				}
		}

	//Log the final simulation parameters
		void LogFile_Finalize(unsigned long long int nParticles,unsigned long long int nScatterings,unsigned long long int nFree,unsigned long long int vcutoff,unsigned long long int nIni,unsigned int nData,unsigned long long int nCrossing0,unsigned long long int nCrossing,double AvDensity,double duration1,double duration2,double duration3,std::ofstream& f)
		{
			double CrossingPerParticle0 = 1.0*nCrossing0/nIni;
			double CrossingPerParticle = 1.0*nCrossing/nParticles;
			std::chrono::time_point<std::chrono::system_clock> end;
			end = std::chrono::system_clock::now();
			std::time_t end_time = std::chrono::system_clock::to_time_t(end);
			f	<<"\n////////////////////////////////////////////////////\n\n"
				<<"//Simulation end: " <<std::ctime(&end_time)<<endl
				<<"//Results" <<endl
				<<"\tData points per iso-ring:\t" <<nData <<endl
				<<"\tSimulated particles:\t\t" <<(double)nParticles <<endl
				<<"\tParticles below vcutoff:\t" <<100.0*vcutoff/nParticles<<" \%"<<endl
				<<"\tFree particles:\t\t\t" <<100.0*nFree/nParticles<<" \%"<<endl
				<<"\tScatterings pp:\t\t\t" <<(double)nScatterings/nParticles<<endl 
				<<"\tDepth crossing pp:\t\t" <<CrossingPerParticle<<"\t("<<CrossingPerParticle0 <<")"<<endl<<endl
				<<"\tAv. Number Density[cm^-3]:\t" <<InUnits(AvDensity,pow(cm,-3.0)) <<endl
				<<"\tAv. Energy Density[GeV cm^-3]:\t" <<InUnits(mChi*AvDensity,GeV/cm/cm/cm) <<endl
				<<"\n//Computation time\n"
				<<"\tInitial run[s]\t\t\t" <<duration1 <<"\t("<<floor(duration1/3600.0)<<":"<<floor(fmod(duration1/60.0,60.0))<<":"<<floor(fmod(duration1,60.0))<<":"<<floor(fmod(1000*duration1,1000.0)) <<")"<<endl
				<<"\tPart/sec\t\t\t" <<nIni/duration1 <<endl <<endl
				<<"\tMain run[s]\t\t\t"<<duration2<<"\t("<<floor(duration2/3600.0)<<":"<<floor(fmod(duration2/60.0,60.0))<<":"<<floor(fmod(duration2,60.0))<<":"<<floor(fmod(1000*duration2,1000.0)) <<")" <<endl
				<<"\tPart/sec\t\t\t" <<nParticles/duration2 <<endl <<endl
				<<"\tTotal[s]\t\t\t"<<duration3<<"\t("<<floor(duration3/3600.0)<<":"<<floor(fmod(duration3/60.0,60.0))<<":"<<floor(fmod(duration3,60.0))<<":"<<floor(fmod(1000*duration3,1000.0)) <<")" <<endl;
				
			f.close();
		}


	//Area of isodetection rings
		double IsoDetectionRing_Area(int theta,int rings,double depth)
		{
			double dTheta = 180.0 / rings;
			return 2*M_PI*pow(rEarth-depth,2)*(cos(theta*deg)-cos((theta+dTheta)*deg));
		}


	std::vector< std::vector<double>> DM_EnergyDensity(long double w0[],int long long unsigned n0total,long double w[],int unsigned long long ntotal,long double w0sq[],long double wsq[],int rings)
	{
		std::vector<std::vector<double>> output;
			for(int i=0;i<rings;i++)
			{
				double dens=1.0*(w[i]/ntotal)/(w0[i]/n0total)*rhoDM;
				double variance = pow(dens/w[i],2.0)*wsq[i] +pow(dens/w0[i],2.0)*w0sq[i]; 
				std::vector <double> buf;
				buf.push_back(dens);
				buf.push_back(sqrt(variance));
				output.push_back(buf);
			}
			return output;
	}


	extern std::vector< std::vector<double>> DM_NumberDensity(double mDM,std::vector<std::vector<double>> density)
	{
		std::vector< std::vector<double>> number_density(density.size(), {0.0,0.0});
		for(unsigned int i=0;i<density.size();i++)
		{
			number_density[i][0] = density[i][0] /mDM;
			number_density[i][1] = density[i][1] / mDM;
		}
		return number_density;
	}



	//Average density over all isodetection rings
		double DM_AverageDensity(std::vector<std::vector<double>> density,int rings,double depth)
		{
			double TotalArea=4*M_PI*pow(rEarth-depth,2);
			double d_theta = 180.0 / rings;
			double sum=0.0;
			for(int i=0;i<rings;i++) 
				sum+= density[i][0]*IsoDetectionRing_Area(i * d_theta,rings,depth)/TotalArea;
			return sum;
		}


//For the analyzer
	void LogFile_Analysis(double duration,int worldsize)
	{
			ofstream f;
			f.open("../results/"+SimID+".log",std::ofstream::app);
			std::chrono::time_point<std::chrono::system_clock> end;
			end = std::chrono::system_clock::now();
			std::time_t end_time = std::chrono::system_clock::to_time_t(end);
			f	<<"\n////////////////////////////////////////////////////\n\n"
			<<"//Data analysis performed:\t" <<std::ctime(&end_time)<<endl
			<<"\tMPI processes:\t\t" <<worldsize<<endl
			<<"\tExperiment:\t\t" <<experiment<<endl
			<<"\tComputation time:\t" <<duration <<"\t("<<floor(duration/3600.0)<<":"<<floor(fmod(duration/60.0,60.0))<<":"<<floor(fmod(duration,60.0))<<":"<<floor(fmod(1000*duration,1000.0)) <<")"<<endl;
			
	}

	std::vector<int> WorkDivision(int WorldSize, int rings)
	{
		std::vector<int> output;
		int overlap = WorldSize*ceil(1.0*rings/WorldSize)-rings;
		for(int i =0 ; i<=WorldSize; i++)
		{
			if(i==0) 			output.push_back(0);
			else if(i<=(WorldSize-overlap))	output.push_back(i*ceil(1.0*rings/WorldSize));
			else				output.push_back(output[i-1]+floor(1.0*rings/WorldSize));
		}
		return output;
	}



