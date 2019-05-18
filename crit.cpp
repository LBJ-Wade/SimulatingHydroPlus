/* Extra inputs required by hydro+ */

//return interpolated correlation length in lattice units
double getintA()
{
  double result;
  long int i = globali;
	double x = globalx;
  if(i!=-1) result = (Aa[i]+x*(Aa[i+1]-Aa[i]))/A;
  else result = xi[0]/A;
  return result;
} 

string label(string title, string element, double num)
{
  int expo = floor(log10(num));
  double dig = floor(num/pow(10.,expo));
  double remainder = num/pow(10.,expo) - dig;
	//control for numerical precision problems
	if(abs(remainder-1) <= .00001)
	{
		dig = dig+1;
		remainder = 0;
	}
	int int_dig = dig;

  if(remainder == 0) return title.append("_").append(element).append("_").append(to_string(int_dig)).append("e").append(to_string(expo));
  else
  {
    string remain = to_string(remainder);
    remain.erase(remain.find_last_not_of('0') + 1, string::npos);
    remain.erase(0,1);
    return title.append("_").append(element).append("_").append(to_string(int_dig)).append(remain).append("e").append(to_string(expo));
  }
}

void set_output_string(string &str, string &dir)
{
	str = string("");
	str = label(str, "Tc", TC);
	str = label(str, "aL", AL);
	str = label(str, "aH", AH);
	str = label(str, "dT", DT);
	str = label(str, "Xm", XM);
	str = label(str, "TL", TL);
	str = label(str, "TH", TH);
	dir = str;
	if(back_react) dir = dir.append("/with_backreaction");
	else dir = dir.append("/without_backreaction");
	cout << "outputing to " << dir << endl;
}

void load_crit_eos()
{
  fstream eosf2;
  //eosf2.open("gregRyanEOS_backReact.dat", ios::in);
	
	cout << string("EOS/richEOSdimensionless").append(param_str).append(".dat") << endl;
  eosf2.open(string("EOS/richEOSdimensionless").append(param_str).append(".dat"), ios::in);

	double Qinv[NUM_MODES];
	double fm_to_lat = 1/.1973/A; //conversion from units of fm to lattice spacing
	for(int i = 0; i<DEL1; ++i) Qinv[i] = (512.*.1973*A - ((double)i)*(512.*.1973*A - 6.)/(1.* DEL1)) * fm_to_lat;
	for(int i = 0; i<DEL2; ++i) Qinv[i+DEL1] = (6. - ((double)i)*(6.-0.5)/(1.* DEL2)) * fm_to_lat;
	for(int i = 0; i<DEL1; ++i) Qinv[i+DEL1+DEL2] = (0.5 - i*0.5/(1.*DEL1)) * fm_to_lat;

	Q = new double[NUM_MODES];
	dQ = new double[NUM_MODES];

  phi = new double*[NUM_MODES];
  for(int i=0; i<NUM_MODES; ++i) phi[i] = new double[NUM+2];

	Phi = new double*[NUM_MODES];
  for(int i=0; i<NUM_MODES; ++i) Phi[i] = new double[NUM+2];
	
  //phi momenta, midpoint interpolation
  for(int j=0; j<NUM_MODES; ++j)
  {
    Q[j] = 1/Qinv[j];
    if(j!=0) dQ[j] = Q[j] - Q[j-1];
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
    	eosf2 >> Aa[i-1];
			eosf2 >> dlogAde[i-1];
			eosf2 >> d2logAde2[i-1];
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
  double norm_fac = 1/(A*A*A*A*A);
  if(i!=-1) result = (dXi[i]+x*(dXi[i+1]-dXi[i]))*norm_fac;
  else result = dXi[0];
  return result;
}

//return interpolated d^2 xi / d eps^2 in lattice units
double getint_d2Xi(int site)
{
  double result;
  double norm_fac = 1/(A*A*A*A);
  long int i = globali;
  double x = globalx;
  norm_fac *= norm_fac/A;

  if(i!=-1) result = (d2Xi[i]+x*(d2Xi[i+1]-d2Xi[i]))*norm_fac;
  else result = d2Xi[0]*norm_fac;
  return result;
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
  return -u[1][site]/u[0][site]*Drphi(i,site) + LAMBDA_M/u[0][site] * 1./phi_eq * (1 - phi[i][site]/phi_eq);
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


//returns p_(+) dp_(+)/r, dp_(+)/dtau
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
  // aa is what's usually called A (xi^-2), Atp is A^tilde prime, Atpp is A^tilde double-prime
  double aa = 1/getintXi(), Atp = getint_Atp(s), Atpp = getint_Atpp(s);
  double phi_eq, phi_p, phi_pp; //d phi_eq / d \epsilon

  for(int i=0; i<NUM_MODES; ++i)
  {
    phi_eq = 1/(Q[i]*Q[i] + aa);
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
