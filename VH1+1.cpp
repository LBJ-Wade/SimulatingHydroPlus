/*
Causal Viscous Hydro Code for Central Heavy Ion Collisions
with back-reaction

by

Paul Romatschke
and G-Nasty

v0.0 January 2007
v0.1 February 2007

The causal viscous hydro equations have been
(re-) derived in 

R.~Baier, P.~Romatschke and U.A.~Wiedemann, arxiv:hep-ph/0602249
(Phys.Rev.C73:064903,2006).

The main algorithm of this code was explained in

R.~Baier and P.~Romatschke, arxiv:nucl-th/0610108

and this version has the updates used in

P.~Romatschke, arxiv:nucl-th-0701032.

If you use this code or some part of it, be sure to reference these
articles.

Copyright of this code is granted provided you keep this disclaimer.

!!! represents a line that I don't understand (most likely terminology or physics related)
//Troubleshoot means that its a hack that should eventually be done away with.

//Questions
1. what is chi_M?
2. How do I define Q_i when using radial coordinates? Bessel-Fourier transform?
3. Do I turn on phi evolution abruptly?
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <Eigen/Dense>

using namespace std;

//typedef Eigen::ArrayXd vec;

#define SQRT3		sqrt(3.0)
#define SQRT2		sqrt(2.0)


// these global vars are initialized from parameters file
// defaults set here are overridden by that file
int    NUM=20, UPDATE=100,SNAPUPDATE=1000;
long int STEPS=40000;
double A=0.05,EPS=0.001,TINIT=1,ETAOS=0.0,TF=0.15,TSTART=0.2;
double EF = .015;
int whichdisc=0;
long int counter = 0;

// controls value of tau_Pi
double COEFF=3.0;

// midpoint and half-width of the critical point
double TC = .2, DELTA = .03;

// anulus width past which critical eos is turned on
double DELTA_R_CRIT = 2;

// Set to .1973 to use proper units
double fac = 1;

// controls rate of phi relaxation
double LAMBDA_M = 1/fac/fac;

// number of phi modes to include
int NUM_MODES = 20;

// these hold the current values
double **u,*e,*pirr,*piee;
double **phi;//, **phi_eq

// these hold the updated values
double **U,*E,*Pirr,*Piee;
double **Phi;//, **Phi_eq

//overall time
double t = 0;

//these are global for convenience; used in doInc
double **dtmat, **rhs;

//also for convenience, defined in loadeos;
long int length,globali;
double globalx;

double *eoT4,*cs2i,*poT4,*Ti;
double *xi,*dXi,*d2Xi,*c_v,*Q,*dQ;

//radius of nucleus in fm
double Rnuc=6.4;
//wood-saxon parameter in fm;
double anuc=0.54;

//flags
int wflag=0;
int reachedTf=0;
bool crit_switch1 = false;
bool crit_switch2 = false;
bool back_react = true;

// output files
fstream T_out;
fstream u_out;
fstream freeze_out;

// initialize global arrays
void allocateMemory() {

    	cout << "==> Allocating memory\n";


	//Fluid velocities u[0]=u^t, u[1]=u^r
	u = new double*[2];
	for (int i=0;i<2;i++) u[i] = new double[NUM+2];

	////Equilibrium phi values
	//phi_eq = new double*[NUM_MODES];
	//for (int i=0;i<NUM_MODES;i++) phi_eq[i] = new double[NUM+2];

	//Energy density
	e = new double[NUM+2];

	//Pi^r_r
	pirr = new double[NUM+2];
	//Pi^\eta_\eta
	piee = new double[NUM+2];

	//critical mode
	phi = new double*[NUM_MODES];
	for(int i=0; i<NUM_MODES; ++i) phi[i] = new double[NUM+2];

	//same as above, at new timestep

	U = new double*[2];
	for (int i=0;i<2;i++) U[i] = new double[NUM+2];

	E = new double[NUM+2];
	Pirr = new double[NUM+2];
	Piee = new double[NUM+2];

	Phi = new double*[NUM_MODES];
	for(int i=0; i<NUM_MODES; ++i) Phi[i] = new double[NUM+2];

	//for convenience only

	dtmat = new double*[3];
	for (int i=0;i<3;i++) dtmat[i] = new double[3];

	double **b;
	rhs = new double*[3];
	for (int i=0;i<3;i++) rhs[i] = new double[1];

	return;
}


void cleanup()
{

  delete [] u;
  delete [] U;

  delete [] e;
  delete [] E;
  delete [] pirr;
  delete [] Pirr;
  delete [] piee;
  delete [] Piee;

  delete [] phi;
  delete [] Phi;
}

// enforces vN boundary conditions by explicitly copying values into the padding elements !!! why von Neumann?
void enforceVNBCs()
{
  for (int i=0;i<2;i++)
    {
      u[i][0]=u[i][1];
      u[i][NUM+1]=u[i][NUM];     
    }

  e[0]=e[1];
  e[NUM+1]=e[NUM];

  pirr[0]=pirr[1];
  pirr[NUM+1]=pirr[NUM];

  piee[0]=piee[1];
  piee[NUM+1]=piee[NUM];

  piee[0]=piee[1];
  piee[NUM+1]=piee[NUM];

  for(int i=0; i<NUM_MODES; ++i)
	{
		phi[i][0]=phi[i][1];
  	phi[i][NUM+1]=phi[i][NUM];
	}
}


void enforceBCs()
{
  enforceVNBCs();
}


// copies uppercase to lowercase letters
void copyDown() {
for (int s=1;s<=NUM;s++) 
  {
    for (int i=0;i<2;i++) u[i][s]=U[i][s];
    e[s]=E[s];
    pirr[s]=Pirr[s];
    piee[s]=Piee[s];
    for(int i=0; i<NUM_MODES; ++i) phi[i][s]=Phi[i][s];
  }
 enforceBCs();
}


//loads the equation of state
void loadeos()
{
  fstream eosf;

  //ideal equation of state
  eosf.open("idEOS.dat", ios::in);

  //Critical EOS
  //eosf.open("gregEOS.dat", ios::in);
  //eosf.open("gregRyanEOS.dat", ios::in);

  //qcd equation of state from Mikko Laine and York Schroeder, 
  //hep-ph/0603048
  //eosf.open("qcdEOS.dat",ios::in);
  //interpolated equation of state -- numerically more stable
  //eosf.open("qcdIEOS.dat",ios::in);

  if (eosf.is_open())
    {
      eosf >> length;
      cout << "Read length of " << length << endl;
      
      eoT4 = new double[length];
      cs2i = new double[length];
      poT4 = new double[length];
      Ti = new double[length];
			xi = new double[length];

			Q = new double[NUM_MODES];
			dQ = new double[NUM_MODES];

			if(back_react)
			{
				dXi = new double[length];
				d2Xi = new double[length];
				c_v = new double[length];
				//chi_Mi = new double[length];
			}


      for (int i=1;i<=length;i++)
					{
						eosf >> Ti[i-1];
						eosf >> eoT4[i-1];
						eosf >> poT4[i-1];
						eosf >> cs2i[i-1];
					}

      eosf.close();
    }
  else
    cout << "Could not open EOS file" << endl; 
}

//Wood-Saxon routine
double WS(double radius, double z)
{
  
  //cout <<"R= " << R << endl;

  double temp;
  temp=radius*radius+z*z;
  temp=sqrt(temp);
  temp-=Rnuc;
  temp/=anuc;
  return 1/(1+exp(temp));
}


//overlap function T_A !!!
double TA(double radius)
{
  double temp=0;
  for (int i=0;i<1200;i++)
    temp+=WS(radius,(i+0.5)*Rnuc/400.)*Rnuc/400.;

  //return result normalized to a gold nucleus
  return 2*temp*197./1175.22;
}

//return interpolated temperature in lattice units
double getintT(int i,double x)
{
	return (Ti[i]+x*(Ti[i+1]-Ti[i]))*A;
}

//provides Temperature*lattice spacing
double T(int site)
{
  //ideal
  double temp;

  long int j;
  //j=geti(e[site]);
  j=globali;

  if (j!=-1)
    return getintT(globali,globalx);
  else
	{
		//Troubleshoot
		//if(e[site] < 0) e[site] = -e[site];
		//!!! Are these the correct units?
    return sqrt(sqrt(e[site]/eoT4[0]));
	}
}

//number of wounded nucleons !!!
double getwnuc(double rr)
{
  double mTA=TA(rr);
  //return mTA*2*(1.-pow(1-mTA/197.*4.,197.));
  
  return mTA*2*(1.-exp(-mTA*4.));
}

//gets index associated with energy density -- internal use

long int geti(double mye)
{
	//Make sure mye is in lattice units at this point
  long int i;
  mye/=A*A*A*A;
  for(i=0;i<length;i++){if(eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i]>mye) break;}
	return (i-1);
}

//get fraction associated with index i !!!
//geti and getx are complementary. give geti an energy, find bin in eos file that contains it.
//give getx an i and energy, find exactly where in bin it resides. 
double getx(long int i,double mye)
{
  mye/=A*A*A*A;
  
  double temp;

	if(i!=-1)
	{
  	temp=(eoT4[i+1]-eoT4[i])/(Ti[i+1]-Ti[i])*(Ti[i+1]+Ti[i])/2;
  	temp+=2*(eoT4[i+1]+eoT4[i]);
  	temp*=(Ti[i+1]+Ti[i])/2;
  	temp*=(Ti[i+1]+Ti[i])/2;
  	temp*=(Ti[i+1]+Ti[i])/2;
  
		return (mye-eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i])/temp/(Ti[i+1]-Ti[i]);
	}
	else return 0;

}

void setInitialConditions()
{
/* initialize u[0] = 1, u[1]=pi=0 everywhere, initialize eps to something like
 * a Wood-Saxon with eps_0 = eps(TSTART) */

  double sig;

  //freeze-out temperature, converted to lattice units
  TF=TF*A;

  cout << "TF=" << TF/A << endl;
  cout << "TSTART=" << TSTART << endl;

  //load equation of state
  loadeos();

  //convert fm/c to lattice units
  t=TINIT*5.06842/A;

  //Saxon-Wood
  double s0;
  
	//find index at which TSTART is closest
  long int i;
  for (i=0;i<length;i++)
  {
    if (Ti[i]>TSTART)
		break;
  }
  
  cout << "found temperature " << Ti[i] << " at " << i << endl;

  //entropy at that temperature in physical units
  s0=(eoT4[i]+poT4[i])*Ti[i]*Ti[i]*Ti[i];

  //energy density at that temperature in physical units

  double e0;

  e0=(eoT4[i])*Ti[i]*A;
  e0*=Ti[i]*A;
  e0*=Ti[i]*A;
  e0*=Ti[i]*A;

  cout << "central entropy " << s0 << endl;

  //normalizing

  double temp;

  temp=Rnuc/anuc;
  temp=exp(temp);
  temp=anuc*log(1+temp);
  temp*=temp;
  s0=e0/temp;
  
	double xiInv;
  for (int s=1;s<=NUM;s++)
  {
    u[0][s]=1.0;
    u[1][s]=0.0;
      
    //energy density WN scaling
    e[s]=e0*getwnuc(s*A/5.06842)/4.29048;

    pirr[s]=0;
    piee[s]=0;
  }

  enforceBCs();
}

//Speed of Sound squared = dp/depsilon
double cs2(double ed)
{

  long int j;
  //j=geti(ed);

  j=globali;

  if (j==-1)
    j=0;

  //return cs2i[j]+globalx*(cs2i[j+1]-cs2i[j]);
  return cs2i[j];
}


//Spatial derivatives
//I'm using a naive version which might be enough 
//for finite viscosity. Would have to cross-check with e.g. RHLLE
//in the future, though...

double Dru(int i,int site)
{
  double temp=0;
  if (site!=1)
    {
      if (site==NUM)
	 {temp=u[i][site]-u[i][site-1];}
       else
	 temp=(u[i][site+1]-u[i][site-1])/2;
    }
  else
    //temp=(2*u[i][site+1]-0.5*u[i][site+2]-1.5*u[i][site])/2.;
    temp=u[i][site+1]-u[i][site];
  return temp;
}

double Dre(int site)
{
  double temp=0;
  if (site!=1)
    {
       if (site==NUM)
	 {temp=e[site]-e[site-1];}
       else
	 temp=(e[site+1]-e[site-1])/2;
    }
  else
    //temp=(2*e[site+1]-0.5*e[site+2]-1.5*e[site])/2.;
    temp=e[site+1]-e[site];
  return temp;
}

double Drpirr(int site)
{
  double temp=0;
  if (site!=1)
    {
      if (site==NUM)
	{temp=(pirr[site]-pirr[site-1]);}
      else
	temp=(pirr[site+1]-pirr[site-1])/2;
    }
  else
    //temp=(2*pirr[site+1]-0.5*pirr[site+2]-1.5*pirr[site])/2.;
    temp=pirr[site+1]-pirr[site];
  return temp;
}

double Drpiee(int site)
{
  double temp=0;
  if (site!=1)
    {
      if (site==NUM) temp=(piee[site]-piee[site-1]);
      else temp=(piee[site+1]-piee[site-1])/2;
    }
  else temp=piee[site+1]-piee[site];
  return temp;
}

//Equation of State; returns p(e)
double eos(double mye, int s)
{
  double temp=0;

  long int j=globali;
  double x = globalx;

  if (j!=-1)
    {
			//temp = (Ti[j]*Ti[j]*Ti[j]*Ti[j]*poT4[j] + x*Ti[j+1]*Ti[j+1]*Ti[j+1]*Ti[j+1]*poT4[j+1])*A*A*A*A;
			//temp = (Ti[j]*Ti[j]*Ti[j]*Ti[j]*poT4[j]+0*x*Ti[j+1]*Ti[j+1]*Ti[j+1]*Ti[j+1]*poT4[j+1])*A*A*A*A;
      temp=Ti[j]*A;
      temp*=temp;
      temp*=temp;
      temp*=poT4[j];
    }
  else
    {
      temp=mye*cs2i[0];
    }

  return temp;
}

//Expansion scalar; returns nabla_mu u^mu
void theta(double *result, int site)
{
  result[0]=1.0;
  result[1]=0.0;
  result[2]=0.0;
  result[3]=u[0][site]/t+u[1][site]/(site)+Dru(1,site);
}

//Provides <\nabla_r u_r>
void nablaurr(double *result, int site)
{
  double th[4];
  theta(th,site);
  result[0]=2.0/3.0*(1+u[1][site]*u[1][site])*th[0];
  result[1]=-2*u[1][site]*u[0][site];
  result[2]=0.0;
  result[3]=2.0/3.0*(1+u[1][site]*u[1][site])*th[3]
    -2*(1+u[1][site]*u[1][site])*Dru(1,site);
}

//Provides t^4 * <\nabla_eta u_eta>
void t4_nablauee(double *result, int site)
{
  double th[4];
  theta(th,site); //I guess t is a global variable
  result[0]=2.0/3.0*th[0]*t*t;
  result[1]=0.0;
  result[2]=0.0;
  result[3]=2.0/3.0*th[3]*t*t
    -2*u[0][site]*t;
}

//Provides r^4 * <\nabla_phi u_phi>
void r4_nablaupp(double *result, int site)
{
  double th[4];
  theta(th,site);
  result[0]=2.0/3.0*th[0]*(site)*(site);
  result[1]=0.0;
  result[2]=0.0;
  result[3]=2.0/3.0*th[3]*(site)*(site)
    -2*u[1][site]*(site);
}

//Provides \beta_2 = tau_Pi/eta/2 \sim 3/4/p
double beta2(int site)
{
	//Muller-Israel-Stewart
	double p;
	//if(back_react) p = crit_eos(e[site], site);
	//else p = eos(e[site]);
	p = eos(e[site], site);
  return COEFF/p/4;
}

//just some debug routine
double getfinals()
{
  double tots=0;
  for (int s=1;s<=NUM;s++)
    {
      globali=geti(e[s]);
      globalx=getx(globali,e[s]);
       //Very approximate entropy:

      tots+=(eoT4[globali]+poT4[globali])*T(s)*T(s)*T(s)*s;
    }
  tots*=t;
  return tots;
}


//Provides tau_Pi/lattice spacing;
double taupi(int site)
{
  double temp=COEFF*2.0*ETAOS/T(site);
  return temp;
}

//Provides \partial_tau pirr
void Dtpirr(double *result,int site,bool verbose=false)
{
	//In this function, collect all coefficients on the RHS
  double vr=u[1][site]/u[0][site];
  double nbrr[4];
  nablaurr(nbrr,site);
  double b2=beta2(site);
  result[0]=-nbrr[0]/b2/2/u[0][site]-2*u[1][site]*vr*pirr[site]; //NOTE: uses 1/2*b2 = eta/tau_PI
  result[1]=-nbrr[1]/b2/2/u[0][site]+2*u[1][site]*pirr[site];
  result[2]=0.0;
  result[3]=-nbrr[3]/b2/2/u[0][site];
  result[3]+=2*u[1][site]*pirr[site]*vr*Dru(1,site);
  result[3]-=2*u[1][site]*vr*pirr[site]*vr*Dru(0,site); 
  result[3]-=vr*Drpirr(site);
  result[3]-=pirr[site]/taupi(site)/u[0][site];
	//Troubleshoot
	if(verbose) 
	{
		printf("nablaurr, taupi, COEFF, ETAOS, T(site): %e, %e, %e, %e, %e \n",nbrr[3], taupi(site), COEFF, ETAOS, getintT(globali,globalx));
    printf("globali, globalx, e, eoT4, %ld %f %e %e \n",globali, globalx, e[site],eoT4[0]);
	}
}

//Provides d_\mu \Pi^\mu_tau
void Dmupimutau(double *result,int site)
{
	//coefficients wrt the RHS
  double vr=u[1][site]/u[0][site];
  double intres[4];
  Dtpirr(intres,site);
  result[0]=intres[0]*vr*vr-2*pirr[site]*vr*vr/u[0][site];
  result[1]=intres[1]*vr*vr+2*pirr[site]*vr/u[0][site];
  result[2]=0.0;
  result[3]=intres[3]*vr*vr;
  result[3]+=pirr[site]*Dru(1,site)/u[0][site];
  result[3]-=pirr[site]*Dru(0,site)/u[0][site]*vr;
  result[3]+=vr*Drpirr(site);
  result[3]+=vr*vr*pirr[site]/t;
  result[3]+=vr*pirr[site]/(site);
  result[3]+=piee[site]/t;

  result[0]*=-1;
  result[1]*=-1;
  result[2]*=-1;
  result[3]*=-1;

}

//Provides d_\mu \Pi^\mu_r
void Dmupimur(double *result,int site,bool verbose=false)
{
  double vr=u[1][site]/u[0][site];
  double intres[4];
  Dtpirr(intres,site);
  result[0]=intres[0]*vr-pirr[site]*vr/u[0][site];
  result[1]=intres[1]*vr+pirr[site]/u[0][site];
  result[2]=0.0;
  result[3]=intres[3]*vr;
  result[3]+=Drpirr(site)+vr*pirr[site]/t;
  result[3]+=(piee[site]+(2-vr*vr)*pirr[site])/(site);
	//Troubleshoot
	if(verbose) printf("Dtpirr, Drpirr, pirr, piee, 1-v^2: %e, %e, %e, %e, %e\n", intres[3], Drpirr(site), pirr[site], piee[site], 1-vr*vr);
}

//provides \Delta^{nu r} d_\mu \Pi^{\mu}_{\nu} 
//(negative of the PI terms in the Du equations in Paul's notes)
void Deltardmupi(double *result,int site)
{
	//coefficients wrt RHS
  double dmupir[4],dmupit[4];
  Dmupimur(dmupir,site);
  Dmupimutau(dmupit,site);
	for(int cc=0; cc<4; ++cc)
	{
  	result[cc]  = -dmupit[cc]*u[1][site]*u[0][site];
  	result[cc] -= dmupir[cc]*(1.0+u[1][site]*u[1][site]);
	}
}

//provides \Delta^{tau \nu} d_\mu \Pi^{\mu}_\nu 
void Deltatdmupi(double *result,int site,bool verbose)
{
	//coefficients wrt RHS
  double dmupir[4],dmupit[4];
  Dmupimur(dmupir,site);
  Dmupimutau(dmupit,site);
	for(int cc=0; cc<4; ++cc)
	{
  	result[cc]  = dmupit[cc]*(1-u[0][site]*u[0][site]);
  	result[cc] -= dmupir[cc]*u[1][site]*u[0][site];
	}
	//Troubleshoot
	if(verbose) printf("dmupit, dmupir, 1-u0^2, u0: %e, %e, %e, %e\n", dmupit[3], dmupir[3], 1-u[0][site]*u[0][site], u[1][site]);
}



//gauss-jordan elimination procedure
//provides 
//void gaussj(double **a, int n,double **b, int m)
#include "GJE.cpp"

//Hydro+ functions
#include "crit.cpp"

//here some diagnostic routines are defined
#include "diags.cpp"

double prev_dpit,prev_dpir;
//main update routine -- this is the core of the code
inline void doInc(double eps, long int iter)
{
  double th[4],nbrr[4],t4_nbee[4],r4_nbpp[4],dpit[4],dpir[4],dtpirr[4],crit_dtp[4];
	double pi_phph = 0;
  double vr=0.0;

  int check=0;

  fstream lout;
  if((iter-1)%SNAPUPDATE == 0)
	{
    char fname[255];
    sprintf(fname,"../data/snapshot/LinAlg_%.2f.dat",t/5.06842*A);
    lout.open(fname, ios::out);
	}

  for (int s=1;s<=NUM;s++)
    {
      globali=geti(e[s]);
      globalx=getx(globali,e[s]);

      vr=u[1][s]/u[0][s];

      //the idea here is that I identify partial_tau u^tau, partial_tau u^r,
      //partial_tau p and constant parts as 0,1,2,3, components
      //of my equation. This can then be rewritten into equations
      //which are explicit in \partial_\tau u^tau,... after inversion
      //of the linear equation system and thus solved
      //see nucl-th/0610108 for details

      theta(th,s);
      nablaurr(nbrr,s);
      t4_nablauee(t4_nbee,s);
      r4_nablaupp(r4_nbpp,s);
			//Troubleshoot
      Deltatdmupi(dpit,s,false);
      Deltardmupi(dpir,s);
      Dtpirr(dtpirr,s);

      //time derivatives in Depsilon=-(e+p)theta...
      //d_t u^t
      dtmat[0][0]=0.0;
      //d_t u^r
      dtmat[0][1]=0.0;
      //d_t eps
      dtmat[0][2]=u[0][s];

      //rhs vector, only first index can differ from 0
      rhs[0][0]=u[1][s]*Dre(s)*(-1.0);

			double p = 0;
			double Drp = 0;
			if(crit_switch2 && back_react)
			{
				p = crit_eos(e[s], s);
				//if(crit_switch2) printf("(p_plus - p)/p: %e\n", (p - eos(e[s]))/eos(e[s]));
				//p = eos(e[s], s);
				//Drp = crit_Drp(e[s], s);
				crit_dp(Drp, crit_dtp, e[s], s);
				//Drp = cs2(e[s])*Dre(s);
				//printf("(Drp_plus - Drp)/Drp: %e\n", (crit_Drp(e[s],s) - Drp)/Drp);
			}
			else 
			{
				p = eos(e[s], s);
				Drp = cs2(e[s])*Dre(s);
			}

			/*matrix coefficients evaluated wrt LHS*/

			//D eps
			pi_phph = -(piee[s] + (1-vr*vr)*pirr[s]);

      rhs[0][0]-=(p+e[s])*th[3];
      rhs[0][0]-=0.5*pirr[s]*(1-vr*vr)*(1-vr*vr)*nbrr[3];
      rhs[0][0]-=0.5*piee[s]*t4_nbee[3]/t/t;
      rhs[0][0]-=0.5*pi_phph/(s)/(s)*r4_nbpp[3];

      for (int cc=0;cc<3;cc++)
			{
				dtmat[0][cc]+=(p+e[s])*th[cc];
				dtmat[0][cc]+=0.5*pirr[s]*(1-vr*vr)*(1-vr*vr)*nbrr[cc];
				dtmat[0][cc]+=0.5*piee[s]*t4_nbee[cc]/t/t;
				dtmat[0][cc]+=0.5*pi_phph/(s)/(s)*r4_nbpp[cc];
			}

      //Du^t=...
      rhs[1][0]=-(p+e[s])*u[1][s]*Dru(0,s);
      rhs[1][0]-=u[0][s]*u[1][s]*Drp;

      dtmat[1][0]=(p+e[s])*u[0][s]; //Derivative of u_t
      dtmat[1][1]=0.0;              //Derivative of u_r
			if(crit_switch2 && back_react && false)
			{
				dtmat[1][2] = -(1-u[0][s]*u[0][s])*crit_dtp[2];
				rhs[1][0] += (1-u[0][s]*u[0][s])*crit_dtp[3];
			}
			else dtmat[1][2]=-cs2(e[s])*(1-u[0][s]*u[0][s]);

      dtmat[1][0]+=dpit[0];
      dtmat[1][1]+=dpit[1];
      dtmat[1][2]+=dpit[2];
      rhs[1][0]-=dpit[3];


      //Du^r=...
      rhs[2][0]=-(p+e[s])*u[1][s]*Dru(1,s);
      rhs[2][0]-=(1.0+u[1][s]*u[1][s])*Drp;

      dtmat[2][0]=0.0;                 //Derivative of u_t
      dtmat[2][1]=(p+e[s])*u[0][s];    //Derivative of u_r
			if(crit_switch2 && back_react && false)
			{
				dtmat[2][2] = u[0][s]*u[1][s]*crit_dtp[2];
				rhs[2][0] -= u[0][s]*u[1][s]*crit_dtp[3];
			}
			else dtmat[2][2]=u[0][s]*u[1][s]*cs2(e[s]);


      dtmat[2][0]+=dpir[0];
      dtmat[2][1]+=dpir[1];
      dtmat[2][2]+=dpir[2];
      rhs[2][0]-=dpir[3];
      
			check=gaussj(dtmat,3,rhs,1);

			if(check==0)
			{
				for (int i=0;i<2;i++) U[i][s]=u[i][s]+eps*rhs[i][0];
				E[s]=e[s]+eps*rhs[2][0];
				Pirr[s]=pirr[s]+eps*(dtpirr[0]*rhs[0][0]+dtpirr[1]*rhs[1][0]+dtpirr[2]*rhs[2][0]+dtpirr[3]);
				Piee[s]=piee[s]+eps/u[0][s]*(-0.5/beta2(s)*(t4_nbee[0]*rhs[0][0]+t4_nbee[1]*rhs[1][0]+t4_nbee[2]*rhs[2][0]+t4_nbee[3])/t/t
												-piee[s]/taupi(s)-u[1][s]*Drpiee(s));
			}
			else 
			{
				if(wflag==0)
				{
					//Troubleshoot
					cout << "Warning: nan at r = "<< s/5.06842*A << "encountered!" << endl;

					char fname[255];
					sprintf(fname,"../data/snapshot/BAD_SPECS_%.2f.dat",t/5.06842*A);
					lout.open(fname, ios::out);
					lout << s/5.06842*A << "\t" << endl;
					lout << dtmat[0][0] << "\t" << dtmat[0][1] <<"\t" << dtmat[0][2] << "\t" << rhs[0][0] << endl;
					lout << dtmat[1][0] << "\t" << dtmat[1][1] <<"\t" << dtmat[1][2] << "\t" << rhs[1][0] << endl;
					lout << dtmat[2][0] << "\t" << dtmat[2][1] <<"\t" << dtmat[2][2] << "\t" << rhs[2][0] << endl;

					lout << "u, Dru, vr: " << u[0][s] << "\t" << u[1][s] << "\t" << Dru(0,s) << "\t" << Dru(1,s) << "\t" << vr << endl;
					lout << "p, e, Dre, cs2: " << p << "\t" << e[s] << "\t" << Dre(s) << "\t" << cs2(e[s]) << endl;
					lout << "pirr, piee, dpit, dpir, prev_dpit, prev_dpir: " << pirr[s] << "\t" << piee[s] << "\t" << dpit[3] << "\t" << dpir[3] << "\t" << prev_dpit << "\t" << prev_dpir << endl;
					lout << "nbrr, nbee, nbpp, theta: " << nbrr[3] << "\t" << t4_nbee[3] << "\t" << r4_nbpp[3] << "\t" << th[3] << endl;
					Deltatdmupi(dpit,s,true);
					Dmupimur(dpit,s,true);
					Dmupimur(dpit,s,true);
					Dtpirr(dpit,s,true);
					lout << endl;
					lout.close();
					snapshot(t);
					wflag=1;
				}
				for (int i=0;i<2;i++) U[i][s]=u[i][s];
				E[s]=e[s];
				Pirr[s]=pirr[s];
				Piee[s]=piee[s];
				if((iter-1)%SNAPUPDATE == 0) lout.close();
				prev_dpit = dpit[3];
				prev_dpir = dpir[3];
			}
			
			double xiInv = 1/getintXi(s);
			double phi_eq;
			//if the eos is critical, evolve the phi derivatives
			if(crit_switch2) for(int i=0; i<NUM_MODES; ++i){
				phi_eq = 1/(Q[i]*Q[i]+xiInv*xiInv);
				Phi[i][s] = phi[i][s] + eps*Dtphi(i,s,phi_eq);
			}
		}
	counter++;
}

//main driver routine
void Evolve()
{
  //setting step sizes to maximum step size
  double eps = EPS;
  long int i=0;
  //for (long int i=1;i<=STEPS;i++) {
  while((reachedTf==0)&&(wflag==0)) {
	  if(i==0)
		{
			outputMeasurements(t); 
      snapshot(t); 
		}
    i++;
    // evolve fields eps forward in time storing updated fields in captial vars

    doInc(eps,i);

    // measurements and data dump
    if ( (i>1 && (i-1)%UPDATE==0)) {
      outputMeasurements(t);
			crit_switch1 = calcDeltaR();
    }
    if ((i-1)%SNAPUPDATE==0) {
      
      snapshot(t); 
    }
	
		//copy fields from capital vars to lowercase vars
		copyDown();
	
		// increment time
		t += eps;

		if(crit_switch1 && !crit_switch2)
		{
			load_crit_eos();
			crit_switch2 = true;
			cout << "Critical eos being used now!!!" << endl;
			SNAPUPDATE=400;
		}
  }
}

void printDivider()
{
	int dwidth = 13;
	for (int i=0;i<8*dwidth;i++) cout << "-"; cout << endl;
	return;
}

int main()
{
	
	extern void readParameters(char*);

	// open files to put the data in first	
	//clengths_out.open("data/clengths.dat", ios::out);
	
	printDivider();
	
	readParameters("params.txt");
	

	printDivider();

	allocateMemory();
		
	setInitialConditions();
	    
	cout << "Initial system entropy" << getfinals() << endl;
       
	//cout << "UV viol init: " << uviol() << endl;

	cout << "check viol init: " << consistcheck1(1) << endl;

	//snapshot(otime);
	
	printDivider();

	cout << "Spatial Step Size (A): " << A << endl;
	
	//cout << "Adaptive stepsize algorithm" << endl;
	
	cout << "Maximal Temporal Step Size (EPS): " << EPS << endl;
	
	cout << "Longitudinal Box Size (A*" << NUM << "): " << A*NUM << endl;

	cout << "Tc, Critical temperature: " << TC << endl;
	cout << "halfwidth of critical region: " << DELTA << endl;
	
	printDivider();
	

	T_out.open("../data/T.dat", ios::out);
	freeze_out.open("../data/freezeout.dat", ios::out);

	printDivider();
	
	int dwidth = 13;
	cout.width(dwidth); cout << "Time";
	cout.width(dwidth); cout << "T(site 1)";
	cout.width(dwidth); cout << "c2(site 1)";
	cout.width(dwidth); cout << "e(site 1)";
	cout.width(dwidth); cout << "T(site NUM)";
	cout.width(dwidth); cout << "u(0, site 1)";
	cout.width(dwidth); cout << "u(1, site 1)";
	cout.width(dwidth); cout << "check";
	cout.width(dwidth); cout << "<UV>";
	cout.width(dwidth); cout << "check2";//"Texp";
	//cout.width(dwidth); cout << "Phi1";
	//cout.width(dwidth); cout << "Phi2";
	cout << endl;
	
	printDivider();
	
	Evolve();
	    
	printDivider();

	cout << "Final system entropy" << getfinals() << endl;

	if (wflag==0)
	  cout << "Done.\n";
	else
	  cout << "Aborted because encountered nan.\n";

	printDivider();
	
	// finally close files
	//clengths_out.close();
	
	T_out.close();
	freeze_out.close();

	cleanup();

	return 0;
}

/*
Change-Log: v0.1: (15.2.2007): Bug in Dtpirr routine fixed
*/
