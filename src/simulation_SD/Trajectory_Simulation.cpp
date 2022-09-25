#include "Trajectory_Simulation.hpp"
#include "Physical_Parameters.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "PREM.hpp"
#include "RN_Generators.hpp"
#include "IC_Generator.hpp"


//This function calculates the final velocity of a DM particle of mass mX after it scattered on a nucleus A with initial velocity vini.
Eigen::Vector3d Scatter(Eigen::Vector3d& vini,double mX,double A,std::mt19937& PRNG)
{
	Eigen::Vector3d vfinal;
	double mNucleus=NucleusMass(A);
	if(FormFactor=="HelmApproximation")
	{
		//Unit Vector in direction of velocity
			Eigen::Vector3d ev=vini.normalized();
		//Determination of n, the unit vector pointing into the direction of vfinal, see Landau Lifshitz.
			Eigen::Vector3d ep(ev(2),ev(2),-(ev(0)+ev(1)));
			ep.normalize();
		//scattering angle is not uniform anymore and depends on the velocity via the form factor
			//double chi = ThetaSample(PRNG);
			long double qmax2 = 4.0*Mu(mX,mNucleus)*Mu(mX,mNucleus)*vini.dot(vini);
			if(qmax2==0.0)cout<<"Error in Scatter with HelmApproximation. (qmax2==0)" <<endl;
			double xi=ProbabilitySample(PRNG);
			long double rA2 =pow( (0.3+0.9*pow(mNucleus/GeV,1.0/3.0))*fm ,2.0);
			long double q2;
			//If the argument of the exp becomes too small exp() returns just 1. and no random angle is coming out. In this case we use a taylor expanded result
			if(rA2*qmax2>1.0e-10)q2 = -3.0/rA2*log(1.0-xi*(1.0-exp(-1.0/3.0*rA2*qmax2)));
			else 	q2 = xi*qmax2;
			if(q2==0.0)cout<<"Error in Scatter with HelmApproximation. (q2==0)" <<endl;
			
			double chi = acos(1.0-2.0*q2/qmax2);

			
		//construct vector mit scattering angle chi
			Eigen::Vector3d n;
			if(chi<M_PI/2)
			{
				n=ev+tan(chi)*ep;
			}	
			else
			{
				n=-ev-tan(chi)*ep;
			}
			n.normalize();
			//We chose a particular ep, therefore we finally rotate n around ev by a random angle.
				double alpha=PhiSample(PRNG);
				n=(n.dot(ev))*ev+cos(alpha)*(ev.cross(n)).cross(ev)+sin(alpha)*ev.cross(n);
			//Output
				vfinal=mNucleus/(mX+mNucleus)*vini.norm()*n+mX/(mX+mNucleus)*vini;	
	}
	else
	{
		//Unit vector with uniformly distributed random direction:
			Eigen::Vector3d n=SphericalCoordinates(1,ThetaSample(PRNG),PhiSample(PRNG));
		//Output
			vfinal=mNucleus/(mX+mNucleus)*vini.norm()*n+mX/(mX+mNucleus)*vini;
	}

	return vfinal;
}

//This function is no longer used
double VelocityCut(double mX,double Ethreshold,double A,double r)
{
	return (mX+NucleusMass(A))/mX*sqrt(Ethreshold/2/NucleusMass(A));
}

//local escape velocity depending on the position in the Earth
double vescape(double r)
{
	//double x=r/rEarth;
	//double ret=0.0000567511*x + 3.86087e-6*pow(x,2) - 0.0000162102*pow(x,3);
	//return ret*3e5*km/sec;
	//return 11.2*km/sec;
	//escape velocity at the top of atmophere, 100km above ground
	return 11.1161*km/sec;
}

//The above method consume too much memeory, just save the last Event
//This function performs the trajectory simulation for a single DM particle and saves the events of scattering in one list.
Trajectory ParticleTrack(double mX,double sigman0,Event IniCondi, double vcut, std::mt19937& PRNG,double rFinal)
{
	//We start outside the earth
		bool InsideEarth=false;
	//Vector for the output:
		vector<Event> EventList;
	//Initial Conditions
		EventList.push_back(IniCondi);
		double t=IniCondi.Time();
		Eigen::Vector3d x=IniCondi.Position();
		Eigen::Vector3d v=IniCondi.Velocity();
		
		//for test
		//double r0 = x.norm();
		//cout << "r0: " << r0/km << " r0/rEarth " << r0/rEarth << endl;
		//exit(0);		
		
	//First we check, if the particle even passes the surface of the earth. This is usually never called upon, since the IC generator doesn't generate particles pointing away from the earth.
		double alphaCone;
		if(abs(x.norm()-rEarthatm)<1.0*meter) alphaCone=M_PI/2.0;//this is necessary in the case that the particles start on the earth surface. then rearth/x.norm can be >1 for numerical reasons.-> sin(1.00000000000001)=nan.
		else alphaCone=asin(rEarthatm/x.norm());
		if (acos(-x.normalized().dot(v.normalized()))<alphaCone)
		{
			double TimeOfEntry =(-x.dot(v)-sqrt(pow(x.dot(v),2)-v.dot(v)*(x.dot(x)-pow(rEarthatm,2))))/(v.dot(v));
			t+=TimeOfEntry;
			x+=v*TimeOfEntry; 
			//FreePathVector() might interpret sitting on top the surface as being outside the Earth due to numerical imprecision.
			//This would result in a particle track without a single scattering, which is why we move the particle 1mm underground at the start.
				if(x.norm()>rEarthatm)	x+=mm*v.normalized();
			//Save point of entry and continue
				EventList.push_back(Event(t,x,v));
				InsideEarth=true;
		}
		//If the particle doesnt pass the earth surface for whatever reason, we just append an additional point and do not enter the loop.
		else 
		{
			t+=4*rEarthatm/v.norm();
			x+=4*rEarthatm*v.normalized();
			EventList.push_back(Event(t,x,v));
		}
	//Simulation Underground:
		while(InsideEarth)
		{
			//New Position
				Eigen::Vector3d xnew=x+FreePathVector(mX,sigman0,x,v,PRNG);
				//cout << "position: " << xnew.norm()/rEarth << " " << xnew(0) << " " << xnew(1) << " " << xnew(2) << endl;
				//cout << "velocity: " << v.norm()/km*sec << " " << v(0) << " " << v(1) << " " << v(2) << endl;				
				if (xnew.norm()>rEarthatm)
				{
					//InsideEarth=false;
					double TimeOfExit =(-x.dot(v)+sqrt(pow(x.dot(v),2)-v.dot(v)*(x.dot(x)-pow(rEarthatm,2))))/(v.dot(v));
					//Save Point of exit
					t+=TimeOfExit;
					x+=TimeOfExit*v;
					if (v.norm()>vescape(rEarthatm))
					{						
						InsideEarth=false;
						//EventList.push_back(Event(t,x,v));
						//Save final point
						t+=rFinal/v.norm();
						x+=rFinal*v.normalized();
						//EventList.push_back(Event(t,x,v));
					}
					else
					{
						//Eigen::Vector3d vopp=-1*v;
						//if not larger than the escape velocity, put the particle 1 millimeter underground to avoid numerical issues and invert its velocity
						v=-1*v;
						x-=1.0*mm*x/x.norm();
						//EventList.push_back(Event(t,x,v));
					}
				}
				
			//Still in Earth! ->Scattering occurs.
				else
				{	
					t+=(xnew-x).norm()/v.norm();
					x=xnew;
					v=Scatter(v,mX,ScatterNucleus(x,PRNG),PRNG);
					//EventList.push_back(Event(t,x,v));
					//are we still above our velocity cut? if not, we abort the loop.
					//if (v.norm()<vcut)
					//if (v.norm()<vescape(x.norm()))
					if (v.norm()<vescape(rEarthatm))
					{
						InsideEarth=false;
					}
				}
		}
		EventList.push_back(Event(t,x,v));
	return Trajectory(EventList);
}







