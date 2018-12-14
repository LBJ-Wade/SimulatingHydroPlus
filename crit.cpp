/* Extra inputs required by hydro+ */

//return interpolated correlation length in lattice units
double getintXi()
{
	//xi is stored in fm.  So we must multiply by 1 = (fac GeV fm)^-1, i.e. divide by fac
  double result;
  long int i = globali;
	double x = globalx;
  if(i!=-1) result = (xi[i]+x*(xi[i+1]-xi[i]))/A;
  else result = 1/A;
  return result/fac;
} 

void load_crit_eos()
{
  fstream eosf2;
  eosf2.open("gregRyanEOS_backReact.dat", ios::in);

	double Qinv[NUM_MODES];
	//for(int i = 0; i<NUM_MODES; ++i) Qinv[i] = 512.*fac*A - i*(512.*fac*A - 0.)/(NUM_MODES);
	for(int i = 0; i<DEL1; ++i) Qinv[i] = 512.*fac*A - i*(512.*fac*A - 6.)/DEL1;
	for(int i = 0; i<DEL2; ++i) Qinv[i+DEL1] = 6. - i*(6.-0.5)/DEL2;
	for(int i = 0; i<DEL1; ++i) Qinv[i+DEL1+DEL2] = 1. - i/DEL1;

	Q = new double[NUM_MODES];
	dQ = new double[NUM_MODES];

  phi = new double*[NUM_MODES];
  for(int i=0; i<NUM_MODES; ++i) phi[i] = new double[NUM+2];

	Phi = new double*[NUM_MODES];
  for(int i=0; i<NUM_MODES; ++i) Phi[i] = new double[NUM+2];
	
  //phi momenta, midpoint interpolation
  for(int j=0; j<NUM_MODES; ++j)
  {
    //Q[j] = 2*M_PI*(j+0)/NUM;//*fac;//!!!
		Q[j] = fac*A/Qinv[j];
    if(j!=0) dQ[j] = Q[j] - Q[j-1];
		//printf("%e %e\n", Q[j], dQ[j]);
  }
  dQ[0] = dQ[1];

  if (eosf2.is_open())
	{
  	eosf2 >> length;

  	for (int i=1;i<=length;i++)
  	{
    	eosf2 >> Ti[i-1];
    	eosf2 >> eoT4[i-1];
    	eosf2 >> poT4[i-1];
    	eosf2 >> cs2i[i-1];
    	eosf2 >> xi[i-1];
			eosf2 >> dXi[i-1];
			eosf2 >> d2Xi[i-1];
			// c_v is useless and should be taken out at some point
			eosf2 >> c_v[i-1];
    }

		//set initial phi to equilibrium value
    double xiInv;
    for (int s=1;s<=NUM;s++)
    {
    	globali = geti(e[s]);
    	globalx = getx(globali, e[s]);
    	xiInv = 1/getintXi();
    	for(int j=0; j<NUM_MODES; ++j) phi[j][s] = 1./(Q[j]*Q[j]+xiInv*xiInv);
    }
    eosf2.close();
  }
  else cout << "Could not open EOS file" << endl;
}

//return interpolated d xi/d eps in lattice units
double getint_dXi(int site)
{
  double result;
	double a = 113.425, c = 14.3269;
	long int i = globali;
	double x = globalx;
  if(i!=-1) result = (dXi[i]+x*(dXi[i+1]-dXi[i]))/(A*A*A*A*A);
	else result = exp(a*T(site)/A - c);
  return result/fac;
}

//return interpolated d^2 xi / d eps^2 in lattice units
double getint_d2Xi(int site)
{
  double result;
	double a = 57.0826, c = 6.29069;
  double norm_fac = 1/(A*A*A*A);
	long int i = globali;
	double x = globalx;
  norm_fac *= norm_fac/A;

  if(i!=-1) result = (d2Xi[i]+x*(d2Xi[i+1]-d2Xi[i]))*norm_fac;
	else result = exp(a*T(site)/A - c);
  return result/fac;
}

/* p_(+) and derivatives of p_(+) and critical modes phi[i] */

double Drphi(int i,int site)
{
  double temp=0;
  if (site!=1)
    {
      if (site==NUM) temp=phi[i][site]-phi[i][site-1];
      else temp=(phi[i][site+1]-phi[i][site-1])/2;
    }
  else temp=phi[i][site+1]-phi[i][site];
  return temp;
}

//Provides d\phi_j for each phi mode of momentum Q_j
double Dtphi(int i, int site, double phi_eq)
{
  return -u[1][site]/u[0][site]*Drphi(i,site) - LAMBDA_M/u[0][site] * phi_eq/phi[i][site] * (1 - phi_eq/phi[i][site]);
}

//returns p_(+)(e)
double crit_eos(double mye, int s)
{
  double result=0, measure=0;
  double ds=0, db=0; //contributions from s_(+) and \beta_(+)
  double xiInv = 1/getintXi();
  double phi_eq, phi_p; //d phi_eq / d \epsilon

  long int j=globali;
  if(j == -1) j = 0;

  for(int i=0; i<NUM_MODES; ++i)
  {
    phi_eq = 1/(Q[i]*Q[i] + xiInv*xiInv);
    phi_p = 2 * phi_eq*phi_eq * xiInv*xiInv*xiInv * getint_dXi(s);
    measure = Q[i]*Q[i]*dQ[i]/(2*M_PI*M_PI);

    ds += .5*measure * (log(phi[i][s]/phi_eq) - phi[i][s]/phi_eq + 1.);
    db += .5*measure * phi_p/phi_eq * (phi[i][s]/phi_eq - 1);
  }

	//printf("%e %e\n", eos(mye,s)+mye, (eos(mye, s) +mye +T(s)*ds)/(1+T(s)*db));
  return (eos(mye, s) + mye + T(s)*ds)/(1. + T(s)*db) - mye;
}


//returns dp_(+)/r, dp_(+)/dtau
void crit_dp(double &dpdr, double *result, double mye, int s)
{
  /* dpdr = dp_(+)/dr,
   * result[0] = coefficient of du^tau/dtau in dp_(+)/dtau
   * result[1] = coefficient of du^r/dtau
   * result[2] = coefficient of de/dtau
   * result[3] = everything else
  */
	result[0]=0;
	result[1]=0;

  double measure=0, dpde=0;
  double ds=0, db=0; //contributions from s_(+) and \beta_(+)
  double xiInv = 1/getintXi(), xi_p = getint_dXi(s), xi_pp = getint_d2Xi(s);
  double phi_eq, phi_p, phi_pp; //d phi_eq / d \epsilon

  for(int i=0; i<NUM_MODES; ++i)
  {
    phi_eq = 1/(Q[i]*Q[i] + xiInv*xiInv);
    phi_p = 2 * phi_eq*phi_eq * xiInv*xiInv*xiInv * xi_p;
    phi_pp = phi_p * (2*phi_p/phi_eq - 3*xi_p*xiInv) + 2*xi_pp * phi_eq*phi_eq * xiInv*xiInv*xiInv;
    measure = Q[i]*Q[i]*dQ[i]/(2*M_PI*M_PI);

    ds += .5*measure * (log(phi[i][s]/phi_eq) - phi[i][s]/phi_eq + 1.);
    db += .5*measure * phi_p/phi_eq * (phi[i][s]/phi_eq - 1);

    dpde += -.5*measure * (phi_p*phi_p*(phi_eq - 2*phi[i][s]) + phi_eq*phi_pp*(phi[i][s] - phi_eq)) / (phi_eq*phi_eq*phi_eq);
  }

  double wplus = (eos(mye, s) + mye + T(s)*ds)/(1. + T(s)*db);
  double bplus = 1/T(s) + db;

  dpde += 1/T(s)/(mye + eos(mye, s)) * cs2(mye);
  dpde *= wplus/bplus;
	result[2] = dpde;


	dpdr = dpde * Dre(s);
	double dpdphi;
	result[3]=0;
  for(int i = 0; i<NUM_MODES; ++i)
	{
    phi_eq = 1/(Q[i]*Q[i] + xiInv*xiInv);
    phi_p = 2 * phi_eq*phi_eq * xiInv*xiInv*xiInv * xi_p;
		measure = Q[i]*Q[i]*dQ[i]/(2*M_PI*M_PI);
		dpdphi = .5*measure/bplus/phi_eq * (-wplus * phi_p/phi_eq + phi_eq/phi[i][s]-1);

		dpdr += dpdphi * Drphi(i,s);
		result[3] += dpdphi * Dtphi(i,s,phi_eq);
	}
	//printf("%e %e\n", dpdr, result[3]);
	//if(counter%10==0) abort();
}
