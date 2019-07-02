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
	double GeV_to_fm = .1973;
	double xi0_to_fm = XI_0;
	double fm_to_lat = 1/.1973/A; //conversion from units of fm to lattice spacing
	for(int i = 0; i<DEL1; ++i) Qinv[i] = (NUM*A*GeV_to_fm - ((double)i)*(NUM*A*GeV_to_fm - 6.)/(1.* DEL1)) * fm_to_lat;
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
	double gamma = GAMMA_0;
	double fm_to_lat = 1/.1973/A; //conversion from units of fm to lattice spacing
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
	double myT = T(s), myp = eos(mye, s), w = mye + myp;
	double eeoT4 = mye/T(s)/T(s)/T(s)/T(s), ppoT4 = getint_poT4(), ccs2 = cs2(mye); //!!!!!
	double woT4 = eeoT4 + ppoT4;

  double logmeasure_o2T3, measure_o2=0, q=0, q_sq=0, Yi=0, Xi=0, wplus_oT4=0, tmp=0, sign=1;
	double dpde_sum=0, dpde, dpdlogphi;
  double ds=0, db=0;
  double A_ = getint_A(), Atp_ = getint_Atp(s), Atpp_ = getint_Atpp(s);
  double phi_eq, phi_dot, phi_2dot;
	double dpdr_term = 0, dpdt_term=0;

  for(int i=0; i<NUM_MODES-2; ++i)
  {
		q = Q[i]*Q[i];
		q_sq = q/A_;
		phi_eq = 1/(q + A_);
		Yi = phi[i][s]/phi_eq - 1;
		phi_dot  = -Atp_/(1+q_sq);														   			 // w * d log phi_eq / de
		phi_2dot = -Atpp_/(1+q_sq) + Atp_*Atp_/(1+q_sq)/(1+q_sq);      // w^2 * d^2 log phi_eq / de^2
		//cout << i << "\t" << q << "\t" << Yi << "\t" << phi_dot << "\t" << phi_2dot << endl;
    //logmeasure_o2T3 = log(q*dQ[i]/(2*M_PI*M_PI)/2) - 3*log(T(s)); // log of integration measure divided by 2 T^3
		measure_o2 = q*dQ[i]/(2*M_PI*M_PI)/2;
		ds           +=  measure_o2 * (log(1+Yi) - Yi);
		//if((i == 119 or i == 118) and s==1) 
		if(Yi+1 <= 0)
		{
			cout << t*.1973*A << "\t" << i << "\t" << s << "\t" << phi[i][s] << "\t" << phi_eq << endl;
			abort();
		}
		db           +=  measure_o2 * phi_dot * Yi / w;
		// integral in expression for dp/de, to be multiplied by w(+)/beta(+)w
		dpde_sum     +=  measure_o2 * (Yi*phi_2dot + (1-Yi) * phi_dot*phi_dot);

		tmp           = -measure_o2 * phi_dot*(1-Yi) / phi[i][s];
		dpdt_term    += tmp * Dtphi(i, s, phi_eq);
		dpdr_term    += tmp * Drphi(i, s);

		result[3]    += measure_o2 * Yi * Dtphi(i, s, phi_eq)/phi[i][s];
		dpdr         += measure_o2 * Yi * Drphi(i, s)/phi[i][s];

		/*
		//To avoid overflow
		tmp = log(Xi) - Yi;
		if(tmp != 0)
		{
			tmp = logmeasure_o2T3 + log(-tmp);
			if(tmp > -30) Tm3ds -= exp(tmp);
		}
		tmp = phi_p*Yi;
		if(tmp != 0)
		{
			sign = tmp/abs(tmp);
			tmp = logmeasure_o2T3 + log(sign*tmp);
			if(tmp>-30) Tdb += sign*exp(tmp);
		}
		tmp = phi_pp*Yi + phi_p*phi_p*(1-Yi);
		//tmp = Yi * (phi_pp * woT4 + phi_p);
		if(tmp!=0)
		{
			sign = tmp/abs(tmp); 
			tmp = logmeasure_o2T3 + log(sign*tmp);
			if(tmp>-30) dpde_sum += sign*exp(tmp);
		}
		*/
  }
	double Tbplus = (1 + myT*db);
	double dp    = myT*(ds - w*db)/Tbplus;
	       pplus = myp + dp;
	double wplus = mye + pplus;

	//if(t*.1973*A > 4.157 and isnan(dp)) cout << s << "\t" << Tbplus << "\t" << ds << "\t" << db << "\t" << myT <<  endl;

	dpdr      += dpdr_term * wplus/w;
	dpdr      *= myT/Tbplus;

	result[3] += dpdt_term * wplus/w;
	result[3] *= myT/Tbplus;

	dpde       = ccs2 + myT*wplus/w/w/Tbplus * dpde_sum;
	result[2]  = dpde;
	dpdr      += dpde * Dre(s);
	//cout << s << "\t" << ds << "\t" << db << "\t" << dp << endl;
	//cout << s << "\t" << mye << "\t" << myT << "\t" << w << "\t" << pplus << "\t"  << Tbplus << "\t" << ccs2 << "\t" << dpde << "\t" << Dre(s) << "\t" << dpdr << endl;


	/*
	wplus_oT4 = (woT4 + Tm3ds)/(1+Tdb);
	pplus = (wplus_oT4)*T(s)*T(s)*T(s)*T(s) - mye;
	//if((s>200) && (s<250)) cout << "s, T (GeV), e/T^4, p/T^4, w/T^4, T^-3 ds, T db: " << s << "\t" << T(s)/A << "\t" << eeoT4 << "\t" << ppoT4 << "\t" << woT4 << "\t" << Tm3ds << "\t" << Tdb << endl;
	if(isnan(pplus)){cout << s*.1973*A << " pplus is nan " << endl; abort();}
	dpde = wplus_oT4/woT4/(1+Tdb)*(ccs2 - woT4 * dpde_sum);
	//dpde = wplus_oT4/woT4/(1+Tdb)*(ccs2 - dpde_sum);
	result[2] = dpde;
	//cout << s << "\t" << ccs2 << "\t" << dpde_sum << "\t" << phi_p*phi_p << endl;

	dpdr = dpde * Dre(s);
	result[3]=0;
	
  for(int i = 0; i<NUM_MODES; ++i)
	{
		q = Q[i]*Q[i];
		q_sq = q/A_;
		phi_eq = 1/(q + A_);
		Yi = phi[i][s]/phi_eq - 1;
		phi_p = Atp_/(1+q_sq);
    measure = q*dQ[i]/(2*M_PI*M_PI);

		dpdlogphi = -.5*T(s)*measure/(1+Tdb) * (Xi*phi_p*wplus_oT4 + Yi);
		dpdr += dpdlogphi * (Drphi(i,s)/phi[i][s]); // !!!!! perhaps make a d log(phi) / dr and d log(phi) / dt
		result[3] += dpdlogphi * Dtphi(i,s,phi_eq)/phi[i][s];
		//cout << s << "\t" << dpdlogphi << "\t" << Drphi(i,s)/phi[i][s] << endl;
	}
	//cout << s << "\t" << dpde << "\t" << dpdlogphi << "\t" << dpdr << "\t" << result[3] << endl;
	*/
	
}
