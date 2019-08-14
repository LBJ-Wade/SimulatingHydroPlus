/* Extra inputs required by hydro+ */

//return interpolated correlation length in lattice units
double getint_A()
{
  double result;
  long int i = globali;
	double x = globalx;
  if(i!=-1) result = (Aa[i]+x*(Aa[i+1]-Aa[i]));
  else result = Aa[0];
  return result;
} 

double getint_Atp(int site)
{
	/* A tilde prime: w d logA / de 
	 * 
	*/
  double result;
  long int i = globali;
	double x = globalx;
  if(i!=-1) result = (Atp[i]+x*(Atp[i+1]-Atp[i]));
  else
	{
		//linear interpolation assuming Atp at T=0 is 0.
		double tmp = T(site);
		result = 0 + (T(site) - 0) * (Atp[0] - 0)/(Ti[0]*A - 0);
	}
  return result;
} 

double getint_Atpp(int site)
{
	/* A tilde double prime: w^2 d^2 logA / de^2 
	 * 
	*/
  double result;
  long int i = globali;
	double x = globalx;
  if(i!=-1) result = (Atpp[i]+x*(Atpp[i+1]-Atpp[i]));
  else
	{
		//linear interpolation assuming Atpp at T=0 is 0.
		double tmp = T(site);
		result = 0 + (T(site) - 0) * (Atpp[0] - 0)/(Ti[0]*A - 0);
	}
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
	str = label(str, "TL", TL*TC);
	str = label(str, "TH", TH*TC);
	dir = str;
	if(SIMTYPE==2)      dir = dir.append("/with_backreaction");
	else if(SIMTYPE==1) dir = dir.append("/without_backreaction");
	else                dir = dir.append("/noCP");
	cout << "outputing to " << dir << endl;
}

void load_crit_eos()
{
  fstream eosf2;
  //eosf2.open("gregRyanEOS_backReact.dat", ios::in);
	
	cout << string("EOS/richEOSdimensionless").append(param_str).append(".dat") << endl;
  eosf2.open(string("EOS/richEOSdimensionless").append(param_str).append(".dat"), ios::in);

	double Qinv[NUM_MODES];
	double invGeV_to_fm = .1973;
	double xi0_to_fm = XI_0;
	double fm_to_lat = 1/.1973/A; //conversion from units of fm to lattice spacing
	for(int i = 0; i<DEL1; ++i) Qinv[i] = (NUM*A*invGeV_to_fm - ((double)i)*(NUM*A*invGeV_to_fm - 6.)/(1.* DEL1)) * fm_to_lat;
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
    //eosf2 >> c1;
    //eosf2 >> c2;
    //eosf2 >> c3;
    //eosf2 >> c4;
		//eosf2 >> c5;
		//cout << length << "\t" << c1 << "\t" << c2 << "\t" << c3 << "\t" << c4 << endl;

  	for (int i=1;i<=length;i++)
  	{
    	eosf2 >> Ti[i-1];
    	eosf2 >> eoT4[i-1];
    	eosf2 >> poT4[i-1];
    	eosf2 >> cs2i[i-1];
    	eosf2 >> Aa[i-1];
			Aa[i-1] = Aa[i-1]/(XI_0*XI_0*fm_to_lat*fm_to_lat); //convert to A in units of lattice spacing
			eosf2 >> Atp[i-1];
			eosf2 >> Atpp[i-1];
    }

		//set initial phi to equilibrium value
    double A_; // A_ = xi^-2
    for (int s=1;s<=NUM;s++)
    {
    	globali = geti(e[s]);
    	globalx = getx(globali, e[s]);
    	A_ = getint_A();
    	for(int j=0; j<NUM_MODES; ++j) phi[j][s] = 1./(Q[j]*Q[j]+A_);
    }
    eosf2.close();
  }
  else cout << "Could not open EOS file" << endl;
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
	//double lam = GAMMA_0*.1973*.1973*.1973*A*A*A;
	double fm_to_lat = 1/.1973/A; //conversion from units of fm to lattice spacing
	double gamma = GAMMA_0/fm_to_lat;
	double tmpA = getint_A();
	double ratio = tmpA * (XI_0*XI_0*fm_to_lat*fm_to_lat);
	double v_r = u[1][site]/u[0][site];
  return -v_r*Drphi(i,site) - gamma/u[0][site] * ratio * (1+Q[i]*Q[i]/tmpA) * (phi[i][site] - phi_eq);
}

//returns p_(+) dp_(+)/r, dp_(+)/dtau
void crit_dp(double &pplus, double &dpdr, double *result, double mye, int s)
{
  /* dpdr = dp_(+)/dr,
   * result[0] = coefficient of du^tau/dtau in dp_(+)/dtau
   * result[1] = coefficient of du^r/dtau
   * result[2] = coefficient of de/dtau
   * result[3] = everything else
  */
	result[0]=0;
	result[1]=0;
	double myT = T(s), myp = eos(mye, s), w = mye + myp, cs2_ = cs2(mye); //!!!!!

  double measure=0, q=0, q_sq=0, Yi=0;
	double dpde_sum=0, dpde;
  double ds=0, db=0;
  double A_ = getint_A(), Adot = getint_Atp(s), A2dot = getint_Atpp(s); //A = xi^-2, Atp = \tilde{A}', Atpp = \tilde{A}''
  double phi_eq, phi_dot, phi_2dot;
	double dpdr_term1 = 0, dpdr_term2 = 0, dpdt_term1=0, dpdt_term2=0;

  for(int i=0; i<NUM_MODES-3; ++i)
  {
		q = Q[i]*Q[i];
		q_sq = q/A_;
		phi_eq = 1/(q + A_);
		Yi = phi[i][s]/phi_eq - 1;
		phi_dot  = -Adot/(1+q_sq);														   			 // w * d log phi_eq / de
		phi_2dot = -A2dot/(1+q_sq) - q_sq*Adot*Adot/(1+q_sq)/(1+q_sq);      // w^2 * d^2 log phi_eq / de^2
		measure = q*dQ[i]/(2*M_PI*M_PI);
		ds           +=  .5*measure * (log(1+Yi) - Yi);
		db           +=  .5*measure * phi_dot * Yi / w;
		dpde_sum     +=  .5*measure * (-Yi*phi_2dot + (1+Yi) * phi_dot*phi_dot);

		dpdt_term1   += -.5*measure * phi_dot*(1+Yi) / phi[i][s] * Dtphi(i, s, phi_eq);
		dpdr_term1   += -.5*measure * phi_dot*(1+Yi) / phi[i][s] * Drphi(i, s);

		dpdt_term2   += -.5*measure * Yi * Dtphi(i, s, phi_eq)/phi[i][s];
		dpdr_term2   += -.5*measure * Yi * Drphi(i, s)/phi[i][s];

  }
	double Tbplus = (1 + myT*db);
	double dp    = myT*(ds - w*db)/Tbplus;
	       pplus = myp + dp;
	double wplus = mye + pplus;

	dpdr       = myT/Tbplus * (wplus/w * dpdr_term1 + dpdr_term2);
	result[3]  = myT/Tbplus * (wplus/w * dpdt_term1 + dpdt_term2);

	dpde       = wplus/w/Tbplus * cs2_ + myT*wplus/w/w/Tbplus * dpde_sum;
	result[2]  = dpde;
	dpdr      += dpde * Dre(s);
}
