#include "PREM.hpp"

#include <cmath>
#include <iostream>
#include <fstream>

#include "General_Utilities.hpp"
#include "RN_Generators.hpp"


//PREM
	//NS: Add atmosphere, up to 100km abovehead
	//Radii
	//							0			1			2			3			4			5			6			7			8			9			10
	const double PREM_Radii[11]={1221.5*km,	3480*km,	5701*km,	5771*km,	5971*km,	6151*km,	6346.6*km,	6356*km,	6368*km,	6371*km, 6471*km};
	//Density Coefficients, we omit the ocean layer with density 1.02, since experiments are located under rocks, not water.
	//NS: Note that the density function for Earth's atmosphere is exponential
	//					0			1			2			3			4			5			6			7		8		9		10		 11	
	//					I.C.		O.C.		L.Mantle 	Trans. 1	Trans. 2	Trans. 3	LVZ&LID		Crust1	Crust2	Ocean	Atmos.   Space	
	const double a[12] ={13.0885,	12.5815,	7.9565,		5.3197,		11.2494,	7.1089,		2.6910,		2.9,	2.6,	2.6,	-2.85,   0};
	const double b[12] ={0,			-1.2638,	-6.4761,	-1.4836,	-8.0298,	-3.8045,	0.6924,		0,		0,		0,		-0.06355,0};
	const double c[12] ={-8.8381,	-3.6426,	5.5283,		0,			0,			0,			0,			0,		0,		0,		0,       0};
	const double d[12] ={0,			-5.5281,	-3.0807,	0,			0,			0,			0,			0,		0,		0,		0,       0};
	//PREM Functions
	int PREM_Layer(double r)
	{
		//The 11 mechanical Layers
		for(int i=10;i>=0;i--)
		{
			if(r>PREM_Radii[i])
			{
				return i+1;
			}
		}
		//Inner Core:
		return 0;
	}
	double PREM_Mass_Density(double r){
		int i = PREM_Layer(r);
		//atmosphere
		if (i==10) {return pow(10.,(a[i]+b[i]*(r-rEarth)/km))*gram*pow(cm,-3);}
		else
		{
			double x = r/rEarth;
			return (a[i]	+b[i]*x+c[i]*pow(x,2)	+d[i]*pow(x,3))*gram*pow(cm,-3);
		}
	}

	
	
//Earth Composition
	//NS: Add atmospheric composition
	const double Elements_Atmos[4][3]=
	{
		{8,16,	0.2318},		//Oxygen		O
		{7,14,	0.7560},		//Nitrogen		N
		{2,4,	7.2472e-7},		//Helium		He
		{18,40,	0.0129}	        //Argon			Ar
	};		
	const double Elements_Crust[9][3]=
	{
		{8,16,	0.467},		//Oxygen		O
		{12,24,	0.021},		//Magnesium		Mg
		{14,28,	0.277},		//Silicon		Si
		{26,56,	0.051},	//Iron			Fe
		{20,40,	0.037},	//Calcium		Ca
		{13,27,	0.081},	//Aluminium		Al
		{11,23,	0.028},	//Natrium		Na
		{19,39,	0.026},	//Potassium		K
		{22,48,	0.006},	//Titanium		Ti		
	};	
	const double Elements_Core[9][3]=
	{
		{26,56,	0.855},		//Iron			Fe
		{14,28,	0.06},		//Silicon		Si
		{28,58,	0.052},		//Nickel		Ni
		{16,32,	0.019},		//Sulfur		S
		{24,52,	0.009},		//Chromium		Cr
		{25,55,	0.003},		//Manganese		Mn
		{15,31,	0.002},		//Phosphorus	P
		{6,12,	0.002},		//Carbon		C
		{1,1,		0.0006}		//Hydrogen		H
	};
	const double Elements_Mantle[14][3]=
	{
		{8,16,	0.440},		//Oxygen		O
		{12,24,	0.228},		//Magnesium		Mg
		{14,28,	0.21},		//Silicon		Si
		{26,56,	0.0626},	//Iron			Fe
		{20,40,	0.0253},	//Calcium		Ca
		{13,27,	0.0235},	//Aluminium		Al
		{11,23,	0.0027},	//Natrium		Na
		{24,52,	0.0026},	//Chromium		Cr
		{28,58,	0.002},		//Nickel		Ni
		{25,55,	0.001},		//Manganese		Mn
		{16,32,	0.0003},	//Sulfur		S
		{6,12,	0.0001},	//Carbon		C
		{1,1,		0.0001},	//Hydrogen		H
		{15,31,	0.00009}		//Phosphorus	P
	};

//Probability to scatter on a certain element. These are calculated in Initialize_PREM().
	double Scatter_Probability_Atmos[4][3];
	double Scatter_Probability_Crust[9][3];
	double Scatter_Probability_Core[9][3];
	double Scatter_Probability_Mantle[14][3];


//calculates the space-independent prefactor and the scatternucleus probabilities for a given parameter point
	//Mean Free Path PreFactor for F(q^2)=1. Including a form factor requires the calculation of these factors at every step of a simulation, because everything becomes velocity depending.
	double g_PREM[4]={0.0,0.0,0.0,0.0};//Core Mantle Crust Atmosphere
	void Initialize_PREM(double mX,double sigma0)
	{
		g_PREM[0]=0.0;
		g_PREM[1]=0.0;
		g_PREM[2]=0.0;
		g_PREM[3]=0.0;
		//1. MFP Prefactor
			//core
				for(int i=0;i<9;i++)
				{
					g_PREM[0]+=Elements_Core[i][2]/NucleusMass(Elements_Core[i][1])*sigmaSI(mX,sigma0,Elements_Core[i][1]);
				}
			//mantle
				for(int i=0;i<14;i++)
				{
					g_PREM[1]+=Elements_Mantle[i][2]/NucleusMass(Elements_Mantle[i][1])*sigmaSI(mX,sigma0,Elements_Mantle[i][1]);
				}
			//crust
				for(int i=0;i<9;i++)
				{
					g_PREM[2]+=Elements_Crust[i][2]/NucleusMass(Elements_Crust[i][1])*sigmaSI(mX,sigma0,Elements_Crust[i][1]);
				}	
			//atmosphere
				for(int i=0;i<4;i++)
				{
					g_PREM[3]+=Elements_Atmos[i][2]/NucleusMass(Elements_Atmos[i][1])*sigmaSI(mX,sigma0,Elements_Atmos[i][1]);
					//cout << "g_PREM: " << g_PREM[3] << endl;
				}					
		//2. Probabilities. The probabilities are velocity depending if form factors are included.
			double total=0.0;
			//core
				total = 0.0;
				for(int i=0;i<9;i++)
				{
					Scatter_Probability_Core[i][0] = Elements_Core[i][0];
					Scatter_Probability_Core[i][1] = Elements_Core[i][1];
					Scatter_Probability_Core[i][2] = Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * sigmaSI(mX,sigma0,Elements_Core[i][1]);
					total+=Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * sigmaSI(mX,sigma0,Elements_Core[i][1]);
				}
				for(int i=0;i<9;i++)
				{
					Scatter_Probability_Core[i][2]/=total;
					//cout <<"Core \t" <<Scatter_Probability_Core[i][0] <<"\t"<<Scatter_Probability_Core[i][1]<<"\t"<<Scatter_Probability_Core[i][2] <<endl;
				}

					
			//mantle
				total = 0.0;
				for(int i=0;i<14;i++)
				{
					Scatter_Probability_Mantle[i][0] = Elements_Mantle[i][0];
					Scatter_Probability_Mantle[i][1] = Elements_Mantle[i][1];
					Scatter_Probability_Mantle[i][2] = Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * sigmaSI(mX,sigma0,Elements_Mantle[i][1]);
					total+=Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * sigmaSI(mX,sigma0,Elements_Mantle[i][1]);
				}
				for(int i=0;i<14;i++)
				{
					Scatter_Probability_Mantle[i][2]/=total;
					//cout <<"Mantle \t" <<Scatter_Probability_Mantle[i][0] <<"\t"<<Scatter_Probability_Mantle[i][1]<<"\t"<<Scatter_Probability_Mantle[i][2] <<endl;
				}


			//crust
				total = 0.0;
				for(int i=0;i<9;i++)
				{
					Scatter_Probability_Crust[i][0] = Elements_Crust[i][0];
					Scatter_Probability_Crust[i][1] = Elements_Crust[i][1];
					Scatter_Probability_Crust[i][2] = Elements_Crust[i][2] / NucleusMass(Elements_Crust[i][1]) * sigmaSI(mX,sigma0,Elements_Crust[i][1]);
					total+=Elements_Crust[i][2] / NucleusMass(Elements_Crust[i][1]) * sigmaSI(mX,sigma0,Elements_Crust[i][1]);
				}
				for(int i=0;i<9;i++)
				{
					Scatter_Probability_Crust[i][2]/=total;
					//cout <<"Crust \t" <<Scatter_Probability_Crust[i][0] <<"\t"<<Scatter_Probability_Crust[i][1]<<"\t"<<Scatter_Probability_Crust[i][2] <<endl;
				}
			//atmosphere
				total = 0.0;
				for(int i=0;i<4;i++)
				{
					Scatter_Probability_Atmos[i][0] = Elements_Atmos[i][0];
					Scatter_Probability_Atmos[i][1] = Elements_Atmos[i][1];
					Scatter_Probability_Atmos[i][2] = Elements_Atmos[i][2] / NucleusMass(Elements_Atmos[i][1]) * sigmaSI(mX,sigma0,Elements_Atmos[i][1]);
					total+=Elements_Atmos[i][2] / NucleusMass(Elements_Atmos[i][1]) * sigmaSI(mX,sigma0,Elements_Atmos[i][1]);
				}
				for(int i=0;i<4;i++)
				{
					Scatter_Probability_Atmos[i][2]/=total;
					//cout <<"Atmosphere \t" <<Scatter_Probability_Atmos[i][0] <<"\t"<<Scatter_Probability_Atmos[i][1]<<"\t"<<Scatter_Probability_Atmos[i][2] <<endl;
				}				

	}

	void Update_PREM(double mX,double sigma0,double velocity)
	{
		if(FormFactor=="HelmApproximation")
		{
			g_PREM[0]=0.0;
			g_PREM[1]=0.0;
			g_PREM[2]=0.0;
			g_PREM[3]=0.0;
			//1. MFP Prefactor
				//core
					for(int i=0;i<9;i++)
					{
						g_PREM[0]+=Elements_Core[i][2]/NucleusMass(Elements_Core[i][1])*TotalsigmaSI(mX,sigma0,Elements_Core[i][1],velocity);
					}
				//mantle
					for(int i=0;i<14;i++)
					{
						g_PREM[1]+=Elements_Mantle[i][2]/NucleusMass(Elements_Mantle[i][1])*TotalsigmaSI(mX,sigma0,Elements_Mantle[i][1],velocity);
					}
				//crust
					for(int i=0;i<9;i++)
					{
						g_PREM[2]+=Elements_Crust[i][2]/NucleusMass(Elements_Crust[i][1])*TotalsigmaSI(mX,sigma0,Elements_Crust[i][1],velocity);
					}
				//atmosphere
					for(int i=0;i<4;i++)
					{
						g_PREM[3]+=Elements_Atmos[i][2]/NucleusMass(Elements_Atmos[i][1])*TotalsigmaSI(mX,sigma0,Elements_Atmos[i][1],velocity);
					}						
			//2. Probabilities. The probabilities are velocity depending if form factors are included.
				double total=0.0;
				//core
					total = 0.0;
					for(int i=0;i<9;i++)
					{
						Scatter_Probability_Core[i][2] = Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Core[i][1],velocity);
						total+=Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Core[i][1],velocity);
					}
					for(int i=0;i<9;i++)
					{
						Scatter_Probability_Core[i][2]/=total;
						//cout <<"Core \t" <<Scatter_Probability_Core[i][0] <<"\t"<<Scatter_Probability_Core[i][1]<<"\t"<<Scatter_Probability_Core[i][2] <<endl;
					}				
				//mantle
					total = 0.0;
					for(int i=0;i<14;i++)
					{
						Scatter_Probability_Mantle[i][2] = Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Mantle[i][1],velocity);
						total+=Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Mantle[i][1],velocity);
					}
					for(int i=0;i<14;i++)
					{
						Scatter_Probability_Mantle[i][2]/=total;
					}
				//crust
					total = 0.0;
					for(int i=0;i<9;i++)
					{
						Scatter_Probability_Crust[i][2] = Elements_Crust[i][2] / NucleusMass(Elements_Crust[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Crust[i][1],velocity);
						total+=Elements_Crust[i][2] / NucleusMass(Elements_Crust[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Crust[i][1],velocity);
					}
					for(int i=0;i<9;i++)
					{
						Scatter_Probability_Crust[i][2]/=total;
					}
				//atmosphere
					total = 0.0;
					for(int i=0;i<4;i++)
					{
						Scatter_Probability_Atmos[i][2] = Elements_Atmos[i][2] / NucleusMass(Elements_Atmos[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Atmos[i][1],velocity);
						total+=Elements_Atmos[i][2] / NucleusMass(Elements_Atmos[i][1]) * TotalsigmaSI(mX,sigma0,Elements_Atmos[i][1],velocity);
					}
					for(int i=0;i<4;i++)
					{
						Scatter_Probability_Atmos[i][2]/=total;
					}						
		}
		else
		{
			cout <<"Error in Update_PREM." <<endl;
		}
	}


	//Return Prefactor:	
	double g_Factor(int layer)
	{
		if(layer<=1) return g_PREM[0];
		else if (layer<=6) return g_PREM[1];
		else if (layer<=9) return g_PREM[2];
		else if (layer<=10) return g_PREM[3];
		else return 0.0;
	}


//Free Path Vector
	//When does a particle leave its current layer?
	double tExit(Eigen::Vector3d& position, Eigen::Vector3d& vel)
	{
			double r=position.norm();
			double v=vel.norm();
			double rv=position.dot(vel);
			int layer = PREM_Layer(r);
			//Check if the particle is well inside the earth.
			//if(layer==10)
			if(layer==11)
			{
				cout <<"Error in tExit(): Position argument lies outside the earth."<<endl;
				return 0;
			}
			//If not we have to distinguish two cases.
			double R[2];
			//Case A:Most often position lies directly on the border between two layers.
				if(r==PREM_Radii[layer])
				{
					//Position directly at the Earth surface:
					//if(layer==9)
					if(layer==10)
					{	
						//Leaving the Earth
						if(rv/r/v>=0) return 1*km/v;
						//Entering the Earth.
						//R[0]=PREM_Radii[9];
						//R[1]=PREM_Radii[8];
						R[0]=PREM_Radii[10];
						R[1]=PREM_Radii[9];						
					}
					//Position directly at core border:
					else if(layer==0)
					{
						R[0]=PREM_Radii[1];
						R[1]=PREM_Radii[0];
					}
					else
					{
						R[0]=PREM_Radii[layer+1];
						R[1]=PREM_Radii[layer-1];
					}
				}
			//Case B: Position is not at the layer border:
				else
				{
					//Core
					if(layer==0) return (-rv+sqrt( pow(rv,2)-pow(v*r,2)+pow(v*PREM_Radii[layer],2) ))/pow(v,2)+10*cm/v;
					else
					{
						R[0]=PREM_Radii[layer];
						R[1]=PREM_Radii[layer-1];
					} 
				}
			//With the two radii we can find tExit.
			double texit=0.0;
			//Outer radius:
				texit = (-rv+sqrt(pow(rv,2)-pow(v*r,2)+pow(v*R[0],2)))/pow(v,2);
			//Inner Radius
				double radic=pow(rv,2)-pow(v*r,2)+pow(v*R[1],2);
				if(radic>=0)
				{
					double t1 = (-rv-sqrt(radic))/pow(v,2);
					double t2 = (-rv+sqrt(radic))/pow(v,2);
					if(t1>0&&t1<texit) texit=t1;
					if(t2>0&&t2<texit) texit=t2;
				}
		return texit+10*cm/v;
	}

	//Exponent of the scattering probability
	//The function form of atmosphere is different from the layers below
	double ProbabilityExponent(Eigen::Vector3d& position,Eigen::Vector3d& vel,double L)
	{
		double r0=position.norm();
		//for test
		//cout << "r0: " << position/km << " r0/km: " << r0/km << endl;
		//cout << "v0: " << vel/km*sec << " v0/km/s: " << vel.norm()/km*sec << endl;
		//exit(0);		
		double v=vel.norm();
		double cosalpha=position.dot(vel)/r0/v;

		int i=PREM_Layer(r0);
		//cout << "layer: " << i << endl;
		double Ltilde=sqrt(L*L+2*L*r0*cosalpha+r0*r0);
		//Coefficients
		if (i==10) //atmosphere, needs numerical integral
		{
			auto f = [i,r0,cosalpha](double x) {return pow(10.,(b[i]*(sqrt(x*x+2*x*r0*cosalpha+r0*r0)-rEarth)/km));};
			double tol = 1e-3;
			int max_refinements = 12;
			//test this function, NewtonMethod does not seem to find a solution
			double fint;
			//avoid trapzoidal error in earlier versions of boost
			if (L==0.) {fint=0.;}
			else {fint = trapezoidal(f, 0.0, L, tol, max_refinements);}
			//double fint = trapezoidal(f, 0.0, L, tol, max_refinements);	
			//cout << "fint: " << fint << endl;
			//exit(0);
			//cout << "gfactor: " << g_Factor(i) << endl;
			return g_Factor(i)*pow(10.,(a[i]))*fint*gram*pow(cm,-3);		
		}
		else
		{
			double C1,C2,C3,C4;
			C1 = L*rEarth*rEarth;
			C3 = L*(r0*r0+r0*L*cosalpha+L*L/3);
			double epsilon=1e-14;
			if(abs(cosalpha+1) > epsilon)
			{
				C2 = rEarth/2*(Ltilde*(L+r0*cosalpha)-r0*r0*cosalpha+(1-cosalpha*cosalpha)*r0*r0*log((L+Ltilde+r0*cosalpha)/((1+cosalpha)*r0)));
				C4 = 1.0/8/rEarth*((5-3*cosalpha*cosalpha)*(cosalpha*pow(r0,3)*Ltilde-pow(r0,4)*cosalpha)+2*L*L*Ltilde*(L+3*r0*cosalpha)+L*Ltilde*r0*r0*(5+cosalpha*cosalpha)+3*pow(r0,4)*pow(1-cosalpha*cosalpha,2)*log((L+Ltilde+r0*cosalpha)/((1+cosalpha)*r0)));
			}
			else
			{
				C2 = (L*(L - 2*r0)*sqrt(pow(L - r0,2))*rEarth)/(2.*(L - r0));
				C4 = (L*(L - 2*r0)*sqrt(pow(L - r0,2))*(pow(L,2) - 2*L*r0 + 2*pow(r0,2)))/(4.*(L - r0)*rEarth);
			}
			return g_Factor(i)/rEarth/rEarth*(a[i]*C1+b[i]*C2+c[i]*C3+d[i]*C4)*gram*pow(cm,-3);
		}
	}
	double ScatterProbability(Eigen::Vector3d& position,Eigen::Vector3d& vel,double L)
	{
		return 1-exp(-ProbabilityExponent(position,vel,L));
	}
	
	//modify newtonmethod, use bisection for atmosphere
	double NewtonMethod(Eigen::Vector3d& position,Eigen::Vector3d& vel,double xi,double Lambda_Total,double NewtonStart=500*km)
	{	
		double r0=position.norm();
		double v=vel.norm();
		double cosalpha=position.dot(vel)/r0/v;
		int i=PREM_Layer(r0);
		double delta=1e-10;
		double L=NewtonStart;
		double f,dfdL,Ltilde;
		//for test
		/*
		cout << "r0: " << r0/km << " layer: " << i << endl;
		cout << " xi: " << xi << " log(1-xi): " << log(1-xi) << " Lambda_tot: " << Lambda_Total << endl;
		cout << "cosalpha: " << cosalpha << endl;
		double sinalpha = sqrt(1-cosalpha*cosalpha);
		double ratio = PREM_Radii[i-1]/PREM_Radii[i];
		cout << "sinalpha: " << sinalpha << "ratio: " << ratio << endl;
		if (sinalpha<ratio)
		{
			double Lmax = r0*cosalpha+sqrt(PREM_Radii[i-1]*PREM_Radii[i-1]-r0*r0*sinalpha*sinalpha);
			double t=tExit(position,vel);
			double Lmax2=t*v;
			cout << "Lmax1: " << Lmax << " Lmax2: " << Lmax2 << endl;
		}
		exit(0);
		*/
		
		if (i==10)
		{
			double t=tExit(position,vel);
			double Lmax=t*v;
			auto f=[&](double x) {return Lambda_Total+ProbabilityExponent(position,vel,x)+log(1-xi);};
			auto tol=[&](double l, double r) {return abs(l-r)<1e-4*km;};
			std::pair<double, double> Lfind=bisect(f,0.,Lmax,tol);
			L=0.5*(Lfind.first+Lfind.second);
			//cout << "bisect L: " << L/km << " Lmax: " << Lmax << " Lfirst: " << Lfind.first/km << " Lsecond: " << Lfind.second/km << endl;
			//cout << "fsol: " << Lambda_Total+ProbabilityExponent(position,vel,L)+log(1-xi) << endl;
			//exit(0);
		}
		else
		{	
			//int k=0;
			while(abs(Lambda_Total+ProbabilityExponent(position,vel,L)+log(1-xi))>delta)
			{
				/*
				k=k+1;
				if (k<=100) {cout << "L: " << L << " difference: " << Lambda_Total+ProbabilityExponent(position,vel,L)+log(1-xi) << endl;}
				else {
					cout << "int L=0: " << Lambda_Total+ProbabilityExponent(position,vel,0) << endl;
					cout << "int L=Inf: " << Lambda_Total+ProbabilityExponent(position,vel,2*rEarth) << endl;
					cout << "int L=-1e20: " << Lambda_Total+ProbabilityExponent(position,vel,-1e20) << endl;
					exit(0);
				}
				*/
				
				Ltilde=sqrt(L*L+2*L*r0*cosalpha+r0*r0);
				f=Lambda_Total+ProbabilityExponent(position,vel,L)+log(1-xi);
				//if (i==10) {dfdL=g_Factor(i)*pow(10.,(a[i]+b[i]*(Ltilde-rEarth)/km))*gram/pow(cm,3)*;}
				//else {dfdL=g_Factor(i)*(a[i]+Ltilde/rEarth*b[i]+pow(Ltilde,2)/(rEarth*rEarth)*c[i]+pow(Ltilde,3)/pow(rEarth,3)*d[i])*gram/pow(cm,3);}
				dfdL=g_Factor(i)*(a[i]+Ltilde/rEarth*b[i]+pow(Ltilde,2)/(rEarth*rEarth)*c[i]+pow(Ltilde,3)/pow(rEarth,3)*d[i])*gram/pow(cm,3);
				L-=f/dfdL;
				
				/*
				double dL=f/dfdL;
				while (L<0.)
				{
					L=L+dL;
					dL=dL/2.;
					L-=dL;
				}
				*/
				//cout << "L: " << L << " dL: " << dL << " f: " << f << " dfdL: " << dfdL << endl;
			}
		}
			
		return L;
	}
	
	//Find the free path vector
	Eigen::Vector3d FreePathVector(double mX,double sigma,Eigen::Vector3d& position,Eigen::Vector3d& vel,std::mt19937& PRNG)
	{
		if(FormFactor!="None") Update_PREM(mX,sigma,vel.norm());
		Eigen::Vector3d r=position;
		double v=vel.norm();
		double xi=ProbabilitySample(PRNG);
		double logxi=-log(1-xi);
		//Add up terms in the exponent of P(L)
		double Lambda = 0;
		int layer=PREM_Layer(r.norm());
		double t,L,Lambda_l,LambdaNew;
		//while(layer<10)
		while(layer<11)
		{
			t=tExit(r,vel);
			L=t*v;
			Lambda_l=ProbabilityExponent(r,vel,L);
			LambdaNew=Lambda+Lambda_l;
			//cout << "Lambdanew: " << LambdaNew << " logxi: " << logxi << " layer:" << layer << endl;
			//No scattering in this layer
			if(LambdaNew<logxi)
			{
				Lambda=LambdaNew;
				r+=t*vel;
				layer=PREM_Layer(r.norm());
			}
			//Scattering in this layer
			else
			{
				//With Newton Method:
				double l =NewtonMethod(r,vel,xi,Lambda, (logxi-Lambda)/Lambda_l*L );
				//double l =NewtonMethod(r,vel,xi,Lambda, 2*rEarthatm );
				r+=l/v*vel;
				//layer=10; //exit while loop
				layer=11; //exit while loop
			}
		}
		return (r-position);
	}



	



//Which nucleus does the DM scatter on?
	 int ScatterNucleus(Eigen::Vector3d& position,std::mt19937& PRNG)
	{
		int layer = PREM_Layer(position.norm());
		double xi = ProbabilitySample(PRNG);
		double sum=0.0;
		//The two compositional layers
		//Mantle
		if(layer>1&&layer<7)
		{
			for(int i=0;i<14;i++)
			{
				sum+=Scatter_Probability_Mantle[i][2];
				if(sum>xi)
				{
					return Scatter_Probability_Mantle[i][1];
				}
			}
		}
		//Crust
		if(layer>6&&layer<10)
		{
			for(int i=0;i<9;i++)
			{
				sum+=Scatter_Probability_Crust[i][2];
				if(sum>xi)
				{
					return Scatter_Probability_Crust[i][1];
				}
			}
		}
		//Atmosphere
		if(layer>9&&layer<11)
		{
			for(int i=0;i<4;i++)
			{
				sum+=Scatter_Probability_Atmos[i][2];
				if(sum>xi)
				{
					return Scatter_Probability_Atmos[i][1];
				}
			}
		}		
		//Core
		else if(layer<2)
		{
			for(int i=0;i<9;i++)
			{
				sum+=Scatter_Probability_Core[i][2];
				if(sum>xi)
				{
					return Scatter_Probability_Core[i][1];
				}
			}
		}
		//Space
		return 0;
	}
