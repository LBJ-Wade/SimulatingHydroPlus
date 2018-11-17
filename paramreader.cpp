#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>


// external vars defined in solve.cpp which are loaded here
extern int    NUM,UPDATE,SNAPUPDATE;
extern long int STEPS;
extern double A,EPS,TINIT,ETAOS,TF,EF,TSTART,COEFF,TC,DELTA,DELTA_R_CRIT,LAMBDA_M,NUM_MODES;

using namespace std;

// this workhorse examines a key to see if corresponds to a var we are setting
// and then attempts to set the var corresponding to key by converting value to the
// appropriate type.  lots of hardcoding here
void setParameter(char *key, char *value) {
	// integer params
	if (strcmp(key,"NUM")==0)          NUM=atoi(value);
	if (strcmp(key,"TINIT")==0)        TINIT=atof(value);
	if (strcmp(key,"STEPS")==0)        STEPS=atoi(value);
	if (strcmp(key,"UPDATE")==0)       UPDATE=atoi(value);
	if (strcmp(key,"SNAPUPDATE")==0)   SNAPUPDATE=atoi(value);
	if (strcmp(key,"A")==0)            A=atof(value);
	if (strcmp(key,"COEFF")==0)        COEFF=atof(value);
	if (strcmp(key,"EPS")==0)          EPS=atof(value);
	if (strcmp(key,"ETAOS")==0)        ETAOS=atof(value);
	if (strcmp(key,"TF")==0)           TF=atof(value);
	if (strcmp(key,"EF")==0)           EF=atof(value);
	if (strcmp(key,"TSTART")==0)       TSTART=atof(value);
	if (strcmp(key,"TC")==0)           TC=atof(value);
	if (strcmp(key,"DELTA")==0)        DELTA=atof(value);
	if (strcmp(key,"DELTA_R_CRIT")==0) DELTA_R_CRIT=atof(value);
	//if (strcmp(key,"LAMBDA_M")==0)     LAMBDA_M=atof(value);
	//if (strcmp(key,"NUM_MODES")==0)    NUM_MODES=atof(value);
	return;
}

//
// This routine assumes that paramters are in text file with
// each parameter on a new line in the format 
//
// PARAMKEY	PARAMVALUE
//
// The PARAMKEY must begin the line and only tabs and spaces
// can appear between the PARAMKEY and PARAMVALUE.
// 
// Lines which begin with 'commentmarker' defined below are ignored
//
void readParameters(char *filename) {
		
	string commentmarker = "//"; 
	char space = ' '; 
	char tab = '\t';

	int maxline = 128; // maximum line length used in the buffer for reading
	char buffer[maxline];
	ifstream paramFile(filename);
	
	while(!paramFile.eof()) {
		paramFile.getline(buffer,maxline,'\n');
		string line = buffer; int length = strlen(buffer);
		if (line.substr(0,commentmarker.length())!=commentmarker && line.length()>0) {
			char key[32]="",value[32]=""; int founddelim=0;
			for (int i=0;i<length;i++) {
				if (buffer[i]==space || buffer[i]==tab) founddelim=1;
				else {
					if (founddelim==0) key[strlen(key)] = buffer[i];
					else value[strlen(value)] = buffer[i];
				}
			}
			if (strlen(key)>0 && strlen(value)>0) {
				setParameter(key,value);
				cout << key << " = " << value << endl;
			}
		}
	}
	
	return;	
}

