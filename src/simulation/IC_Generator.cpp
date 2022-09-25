#include "IC_Generator.hpp"

#include <iostream>
#include <fstream>

#include "RN_Generators.hpp"
#include "General_Utilities.hpp"

#include "Trajectory_Simulation.hpp"

//Two definitions for the default values
Event InitialCondition(std::mt19937& PRNG){
	return InitialCondition(0,vEarth,PRNG);
}
Event InitialCondition(double t,Eigen::Vector3d &vearth,std::mt19937& PRNG,double R){
	Eigen::Vector3d IniPosi,IniVeli;
	//Initial Velocity
		double vnorm=MaxwellSample(PRNG);
		IniVeli=SphericalCoordinates(vnorm,ThetaSample(PRNG),PhiSample(PRNG));
		//Boost into the geocentric frame:
			IniVeli-=vearth;
		//consider earth gravitational potential
		Eigen::Vector3d vnormdtmp=IniVeli.normalized();
		double vnormtmp=IniVeli.norm();
		double vsunEarth=42.2256*km/sec;//including solar acceleration at the position of the Earth
		double vnormnew=sqrt(vnormtmp*vnormtmp+vescape(rEarthatm)*vescape(rEarthatm)+vsunEarth*vsunEarth);
		IniVeli=vnormnew*vnormdtmp;
	//Random Point in a circle at distance R
		Eigen::Vector3d ez=-IniVeli.normalized();
		Eigen::Vector3d ex(0,ez[2]/sqrt(ez[1]*ez[1]+ez[2]*ez[2]),-ez[1]/sqrt(ez[1]*ez[1]+ez[2]*ez[2]));
		Eigen::Vector3d ey=ez.cross(ex);
		double phi=PhiSample(PRNG);
		double xi=ProbabilitySample(PRNG);
		//IniPosi=R*ez+pow(xi,1.0/2.0)*rEarth*(cos(phi)*ex+sin(phi)*ey);
		//Change the initial position to be on the top of atmosphere
		IniPosi=R*ez+pow(xi,1.0/2.0)*rEarthatm*(cos(phi)*ex+sin(phi)*ey);
	//Projection onto a sphere around earth (optional, but I have to do the projection also in time!)
		// double tProject = 1/pow(IniVeli.norm(),2)*(-IniPosi.dot(IniVeli)-sqrt(pow(IniPosi.dot(IniVeli),2)+pow(IniVeli.norm(),2)*(R*R-pow(IniPosi.norm(),2))));
		// IniPosi+=tProject*IniVeli;
		// t+=tProject;
		
		//if (acos(-IniPosi.normalized().dot(IniVeli.normalized()))>asin(rEarth/IniPosi.norm()))
		if (acos(-IniPosi.normalized().dot(IniVeli.normalized()))>asin(rEarthatm/IniPosi.norm()))
			{
				cout <<"Error in InitialCondition(): Send in particle misses the Earth." <<endl;
			}
	//Return the result

		return Event (t,IniPosi,IniVeli);
}











