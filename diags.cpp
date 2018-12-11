double uviol()
{
  double tmp=0;
  for (int s=1;s<=NUM;s++)
    { 
      tmp+=u[0][s]*u[0][s]-u[1][s]*u[1][s];
    }
  tmp/=NUM;
  return sqrt(tmp)-1;
}

double consistcheck1(int site)
{
  return 2*pirr[site]+piee[site];
}

double check()
{
  double temp=0;
  for (int s=1;s<=NUM;s++)
    temp+=fabs(consistcheck1(s));
  return temp/NUM;
}


//constant temperature freeze-out
int freezesite()
{
  int temp=1;

  //search for freeze-out point
  for (int s=2;s<=NUM;s++)
    {
      globali=geti(e[s]);
      globalx=getx(globali,e[s]);
      if (T(s)>TF)
	  temp++;
    }
  if (temp==1)
      reachedTf=1;
  return temp;
}

//mean-free path freeze-out
int freezenew()
{
  int temp=1;

  double limiter=t;

  if (Rnuc*5.06842/A<t)
    limiter=Rnuc*5.06842/A;

  //search for freeze-out point
  for (int s=2;s<=NUM;s++)
    {
      globali=geti(e[s]);
      globalx=getx(globali,e[s]);
      if (T(s)>ETAOS*20/limiter)
	temp++;
    }
  if (temp==1)
    reachedTf=1;
  return temp;
}

void outputMeasurements(double t) 
{

  double uv=uviol();

 
  int fsite=freezesite();

  // output to screen
  cout.precision(5);
  int dwidth = 13;
  cout.width(dwidth); cout << t/5.06842*A;
  globali=geti(e[1]);
  globalx=getx(globali,e[1]);
  cout.width(dwidth); cout << T(1)/A;
  cout.width(dwidth); cout << cs2(e[1]);
  globali=geti(e[NUM]);
  globalx=getx(globali,e[NUM]);
  cout.width(dwidth); cout << T(NUM)/A;
  globali=geti(e[1]);
  globalx=getx(globali,e[1]);
  cout.width(dwidth); cout << u[0][1];
  cout.width(dwidth); cout << u[1][1];
  cout.width(dwidth); cout << check();
  cout.width(dwidth); cout << uv;
  cout.width(dwidth); cout << consistcheck1(1);
  //cout.width(dwidth); cout << 0.27453*exp(log(0.1973/t)/3);
  //cout.width(dwidth); cout << -piee[1]*t*t/4/p[1];
  //cout.width(dwidth); cout << Phi[1]/4/p[1];
  //cout.width(dwidth); cout << p[1]/T(1)/T(1)/T(1)/T(1);
  //cout.width(dwidth); cout << taupi(1)/2/beta2(1)/4/1.7546/T(1)/T(1)/T(1);
  cout << endl;

  T_out << t/5.06842*A << "\t" << T(1)/A << endl;

  double vr=u[1][fsite]/u[0][fsite];
  globali=geti(e[fsite]);
  globalx=getx(globali,e[fsite]);
  freeze_out << fsite/5.06842*A << "\t" << t/5.06842*A << "\t";
  freeze_out << u[0][fsite] << "\t";
  freeze_out << u[1][fsite] << "\t";
  freeze_out << pirr[fsite]/(e[fsite]+eos(e[fsite], fsite)) << "\t";
  freeze_out << piee[fsite]/(e[fsite]+eos(e[fsite], fsite)) << "\t";
  freeze_out << eos(e[fsite], fsite)/(A*A*A*A) << "\t";
  freeze_out << fabs(pirr[fsite]*3*vr*vr-piee[fsite])/(eos(e[fsite], fsite)+e[fsite])/4 << "\t";
  freeze_out << 1/sqrt(fabs(pirr[fsite]*3*vr*vr-piee[fsite])/(eos(e[fsite], fsite)+e[fsite])/4)*T(fsite)/A << "\t";
  freeze_out << T(fsite)/A << "\t";
  freeze_out << (eoT4[globali]+poT4[globali])*T(fsite)*T(fsite)*T(fsite)*fsite*t/A << endl;

  //gausslaw_out << time << "\t" << log10(gausslaw) << endl; 
}


void snapTprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"../data/snapshot/Tprofile_%.2f.dat",time/5.06842*A);
  out.open(fname, ios::out);
  for (int s=1;s<=NUM;s++)
  {
    globali=geti(e[s]);
    globalx=getx(globali,e[s]);
    out << s/5.06842*A << "\t";
    out << T(s)/A << endl;
  }
  out.close();
}

void snapEDprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"../data/snapshot/edprofile_%.2f.dat",time/5.06842*A);
  out.open(fname, ios::out);
  for (int s=1;s<=NUM;s++)
  {  
    out << s/5.06842*A << "\t";
    out << e[s]/A/A/A/A << endl;
  }
  out.close();
}

void snapvprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"../data/snapshot/vprofile_%.2f.dat",time/5.06842*A);
  out.open(fname, ios::out);
  for (int s=1;s<=NUM;s++)
  {
    out << s/5.06842*A << "\t";
    out << u[1][s]/u[0][s] << endl;
  }
  out.close();
}

void snapRinvprofile(double time)
{
  fstream out;
  char fname[255];
  double vr=0;
  sprintf(fname,"../data/snapshot/Rinvprofile_%.2f.dat",time/5.06842*A);
  out.open(fname, ios::out);
  for (int s=1;s<=NUM;s++)
  {
    globali=geti(e[s]);
    globalx=getx(globali,e[s]);
    vr=u[1][s]/u[0][s];
    out << s/5.06842*A << "\t";
    out << fabs(pirr[s]*3*vr*vr-piee[s])/(eos(e[s], s)+e[s])/4 << endl;
  }
  out.close();
}


void snappirrprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"../data/snapshot/Pirrprofile_%.2f.dat",time/5.06842*A);
  out.open(fname, ios::out);
  for (int s=1;s<=NUM;s++)
  {
    globali=geti(e[s]);
    globalx=getx(globali,e[s]);
    out << s/5.06842*A << "\t";
    out << pirr[s]/(eos(e[s], s)+e[s]) << "\t";
    out	<< pirr[s] << endl;
  }
  out.close();
}

void snappieeprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"../data/snapshot/Pieeprofile_%.2f.dat",time/5.06842*A);
  out.open(fname, ios::out);
  for (int s=1;s<=NUM;s++)
  {
    globali=geti(e[s]);
    globalx=getx(globali,e[s]);
    out << s/5.06842*A << "\t";
    out << piee[s]/(eos(e[s], s)+e[s]) << "\t";
    out	<< piee[s] << endl;
  }
  out.close();
}

void snapPhiprofile(double time)
{
  //if(time>10.77) printf("\n Inside Man!!! \n");
  fstream out;
  char fname[255];
  sprintf(fname,"../data/snapshot/Lam%.1f_Phiprofile_%.3f.dat", LAMBDA_M, time/5.06842*A);
  out.open(fname, ios::out);

  fstream out2;
  char fname2[255];
  sprintf(fname2,"../data/snapshot/entropyCorrection_%.3f.dat", time/5.06842*A);
  out2.open(fname2, ios::out);

  double phi_eq, xiInv, entropy, measure;
  for (int s=1;s<=NUM;s++)
  {
    globali=geti(e[s]);
    globalx=getx(globali,e[s]);
    xiInv = 1/getintXi();
		double ds=0;
    out << s/5.06842*A << "\t";

    for(int i=0; i<NUM_MODES; ++i)
    {
      phi_eq = 1/(Q[i] * Q[i] + xiInv * xiInv);
			measure = Q[i]*Q[i]*dQ[i]/(2*M_PI*M_PI);
      ds += .5*measure * (log(phi[i][s]/phi_eq) - phi[i][s]/phi_eq + 1.);
      //if(i != 0) phi_eq = 1/(Q[i] * Q[i] + xiInv * xiInv);
      //else phi_eq = globalx*1.;
      out << Q[i] << "\t" << phi_eq << "\t" << phi[i][s] << "\t";
			if(i<5 && false)
			{
				printf("contribution %e, %e, %e\n", measure, log(phi[i][s]/phi_eq), phi[i][s]);
			}
    }
		entropy = (e[s] + eos(e[s], s))/T(s);
		out2 << s*fac*A << "\t" << ds/entropy << endl;
    out << endl;
  }
  out.close();
	out2.close();
}

void snapPplusProfile(double time)
{
  if(crit_switch && back_react)
  {
    fstream out;
    char fname[255];
    sprintf(fname,"../data/snapshot/pplusprofile_%.3f.dat",time/5.06842*A);
    out.open(fname, ios::out);
    double fac = 1/A/A/A/A/A;
    for (int s=1;s<=NUM;s++)
    {
			globali = geti(e[s]);
			globalx = getx(globali, e[s]);
      out << s/5.06842*A << "\t";
      out << eos(e[s], s)/A/A/A/A << "\t" << crit_eos(e[s], s)/A/A/A/A << "\t" << globalx << endl;//<< "\t" << crit_Drp(e[s],s)*fac << endl;
    }
    out.close();
  }
}

void snapshot(double time)
{
  snapTprofile(time);
  snapvprofile(time);
  snapRinvprofile(time);
  snapEDprofile(time);
  if(crit_switch) 
	{
		snapPhiprofile(time);
		snapPplusProfile(time);
	}

  //snappirrprofile(time);
  //snappieeprofile(time);
}
