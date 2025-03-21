// February 28, 2025: MCMC_arctrace01d.cpp: Use an MCMC to fit
// an orbit to input data in hldet format. Position-only distance metric
// Try to calculate the covariance matrix of the state vectors.
// Output RA, Dec, and
// angular velocities rather than state vectors.
//
// Uses fitting based on arctrace01_6Dh.cpp, which is
// like arctrace01_6D.cpp, but takes its input file in heliolinc hldet
// format, and requires an observatory code file.
//
// Description of ancestor program arctrace01_6D.cpp
// Given an input file giving observer positions and measured RA and Dec
// for an object at various times, given as MJD, fit the object's orbit
// based on a starting approximate state vector, and the MJD at which this
// matched the object's position. The input MJD must fall exactly on
// on step of the planet files read from the configuration file.
// Perform a full, 6-dimensional orbit fit
// including planetary perturbations. 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define MINCOLS 9
#define DEBUG 0
#define RESETNUM 2000
#define MINACCEPT 0.2
#define MAXACCEPT 0.5
#define POSSTEP 1000.0
#define PRINTNUM 1000
#define COVARNUM 10000
#define VELTIMESCALE 1.0e7l
#define EIGENOFFMAX 1e-100
#define EIGENITMAX 1000000
#define REPNUM 10

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain
// MJD, observer X, Y, Z, VX, VY, VZ, RA, and Dec. The observer coordinates and velocities must
// be in km and km/sec, relative to the Sun. The RA and Dec must be in decimal degrees.
static void show_usage()
{
  cerr << "Usage: MCMC_arctrace01e -cfg configfile -ranseed random_number_seed -observations obsfile -mjd starting_mjd -statevec x y z vx vy vz -maxnum maxnum -repnum repnum -obscode obscodefile -outfile outfile -repfile repfile\n";
}

vector <long double> posvel2statevec01(point3LD pos, point3LD vel);
int statevec2posvel01(const vector <long double> &invec, point3LD &pos, point3LD &vel);
vector <long double> posvel2statevec02(point3LD pos, point3LD vel, long double timescale);
int statevec2posvel02(const vector <long double> &invec, point3LD &pos, point3LD &vel, long double timescale);

int main(int argc, char *argv[])
{
  ifstream instream1;
  int status, planetnum, planetct, pctSun, pctEarth, configread;
  long maxnum,itct,resetct,i, j, k;
  maxnum = 10000;
  itct = resetct = 0;
  int badread,reachedeof,obsnum,obsct;
  i=j=k=status=configread=badread=reachedeof=obsnum=obsct=0;
  int polyorder = 3;
  long covarnum = COVARNUM;
  long covarct = 0;
  vector <hldet> detvec = {};
  vector <observatory> observatory_list = {};
  string obscodefile;
  double obslon = 0.0l;
  double plxcos = 0.0l;
  double plxsin = 0.0l;
  string stest;
  string configfile;
  string obsfile;
  string planetfile;
  string outfile = "outjunk01a.txt";
  string repfile = "repjunk01a.txt";
  double kepspan=15.0;
  long double ldval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> observer_barypos;
  vector <long double> obsMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> sigastrom;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> resid;
  string lnfromfile;
  vector <string> linestringvec;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  point3LD obspos = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  point3LD testpos = point3LD(0,0,0);
  point3LD testvel = point3LD(0,0,0);
  vector <long double> beststatevec;
  double chisq = 0.0;
  double newchi = 0.0;
  double bestchi = LARGERR;
  vector <double> bestRA;
  vector <double> bestDec;
  vector <double> bestresid;
  long double astromrms = 1.0;
  long double mjdstart = 0.0l;
  int planetfile_refpoint=0;
  int verbose = 0;
  int changepar = 0;
  string seedstring;
  vector <double> stepvec;
  vector <double> keepvec;
  vector <double> normvec;
  make_dvec(6, stepvec);
  make_dvec(6, keepvec);
  make_dvec(6, normvec);
  stepvec[0] = stepvec[1] = stepvec[2] = POSSTEP;
  stepvec[3] = stepvec[4] = stepvec[5] = POSSTEP;
  vector <point3LD> orbit05pos;
  vector <point3LD> orbit05vel;
  vector <long double> orbit05MJD;
  long double RA1,Dec1,RAvel,Decvel;
  double bestephRA,bestephDec,bestephRAvel,bestephDecvel;
  bestephRA = bestephDec = bestephRAvel = bestephDecvel = 0;
  vector <long double> startstate;
  vector <long double> tempstate;
  vector <long double> beststate;
  vector <long double> onestate;
  vector <vector <long double>> accepted_states;
  vector <vector <long double>> representative_states;
  vector <vector <long double>> accepted_all;
  vector <vector <long double>> covariance_matrix;
  long double timescale = VELTIMESCALE;
  vector <vector <long double>> eigenvecs;
  vector <long double> eigenvals;
  vector <long double> statedist;
  vector <long double> statechi;
  long double maxdist;
  long distantstate = 0;
  long repct = 0;
  long repnum = REPNUM;
 
  make_LDvec(6,onestate);
  make_LDmat(6,6,covariance_matrix);
  make_LDmat(6,6,eigenvecs);


  if(argc<9) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-c" || string(argv[i]) == "-cfg" || string(argv[i]) == "-config" || string(argv[i]) == "--config" || string(argv[i]) == "--configfile" || string(argv[i]) == "--configuration" || string(argv[i]) == "--ConfigFile") {
      if(i+1 < argc) {
	//There is still something to read;
	configfile=argv[++i];
	i++;
	// Read the configuration file. This must happen here, so that
	// default values supplied in this file can be overwritten later
	// if the user desires.
	// Read configuration file.
	ifstream instream1 {configfile};
	if(!instream1) {
	  cerr << "ERROR: can't open input config file " << configfile << "\n";
	  return(1);
	}
	// Read polyorder
	status=readconfigint(instream1,&polyorder);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&polyorder);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Polynomial order for orbit integration read as " << polyorder << "\n";
	// Read the number of planets.
	status=readconfigint(instream1,&planetnum);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&planetnum);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Number of planets read as " << planetnum << "\n";
	// Read the index of the Sun within the planet vectors
	status=readconfigint(instream1,&pctSun);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctSun);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "The Sun is planet number " << pctSun << "\n";
	// Read the index of Earth within the planet vectors
	status=readconfigint(instream1,&pctEarth);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctEarth);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Earth is planet number " << pctEarth << "\n";
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the mass for planet number planetct
	  status=readconfigLD(instream1,&ldval);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigLD(instream1,&ldval);
	  }
	  planetmasses.push_back(ldval);
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
	}
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the ephemeris file for planet number planetct
	  status=readconfigstring(instream1,planetfile);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigstring(instream1,planetfile);
	  }
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "Ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
	  mjdtest={};
	  temppos={};
	  tempvel={};
	  read_horizons_fileLD(planetfile,mjdtest,temppos,tempvel);
	  if(planetct==0) planetmjd=mjdtest;
	  else {
	    for(j=0;j<long(planetmjd.size());j++) {
	      if(mjdtest[j]!=planetmjd[j]) {
		cout << "ERROR: time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(1);
	      }
	    }
	  }
	  for(j=0;j<long(temppos.size());j++) {
	    planetpos.push_back(temppos[j]);
	  }
	  if(planetct == pctEarth) {
	    Earthpos = temppos;
	    Earthvel = tempvel;
	  }
	  if(planetct == pctSun) {
	    Sunpos = temppos;
	    Sunvel = tempvel;
	  }
	  cout << "Finished reading ephemeris file " << planetfile << "\n";
	}
	// Close input stream that was reading the config file.
	instream1.close();
	configread=1;
      } else {
	cerr << "Configuration file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }
    // Stop execution if we haven't read a config file successfully by now.
    if(configread!=1) {
      cerr << "ERROR: configuration file must be supplied and successfully";
      cerr << "read before any of the other parameters. This is because\n";
      cerr << "it contains defaults that can optionally be overridden by\n";
      cerr << "user-supplied values in later arguments.\n";
      return(2);
    } else cout << "Configuration file read successfully\n";
    if(string(argv[i]) == "-rs" || string(argv[i]) == "-ranseed" || string(argv[i]) == "-rseed" || string(argv[i]) == "-seed" || string(argv[i]) == "--ranseed" || string(argv[i]) == "--randomseed" || string(argv[i]) == "--randomnumberseed" || string(argv[i]) == "--random_number_seed") {
      if(i+1 < argc) {
	// There is still something to read;
	seedstring=argv[++i];
	i++;
      } else {
	cerr << "Random number seed keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-o" || string(argv[i]) == "-obs" || string(argv[i]) == "-obsfile" || string(argv[i]) == "-observations" || string(argv[i]) == "--obsfile" || string(argv[i]) == "--observation" || string(argv[i]) == "--observationfile" || string(argv[i]) == "--observations") {
      if(i+1 < argc) {
	// There is still something to read;
	obsfile=argv[++i];
	i++;
      } else {
	cerr << "Observation file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjd") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdstart=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Starting MJD keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-statevec" || string(argv[i]) == "-sv" || string(argv[i]) == "-svec" || string(argv[i]) == "--statevector" || string(argv[i]) == "--statevec" || string(argv[i]) == "-statevector" || string(argv[i]) == "--sv") {
      if(i+1 < argc) {
	//There is still something to read;
	startpos.x=stold(argv[++i]);
	cout << "read sv entry as " << startpos.x << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startpos.y=stold(argv[++i]);
	cout << "read sv entry as " << startpos.y << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startpos.z=stold(argv[++i]);
	cout << "read sv entry as " << startpos.z << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
 	startvel.x=stold(argv[++i]);
	cout << "read sv entry as " << startvel.x << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startvel.y=stold(argv[++i]);
	cout << "read sv entry as " << startvel.y << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startvel.z=stold(argv[++i]);
	cout << "read sv entry as " << startvel.z << "\n";
	i++;
      }
      else {
	cerr << "State vector keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxnum") {
      if(i+1 < argc) {
	//There is still something to read;
	maxnum=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Max interation number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-repnum") {
      if(i+1 < argc) {
	//There is still something to read;
	repnum=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Representative statevec number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obscode" || string(argv[i]) == "-obs" || string(argv[i]) == "-oc" || string(argv[i]) == "-obscodes" || string(argv[i]) == "--obscode" || string(argv[i]) == "--obscodes" || string(argv[i]) == "--observatorycodes") {
      if(i+1 < argc) {
	//There is still something to read;
	obscodefile=argv[++i];
	i++;
      }
      else {
	cerr << "Observatory code file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outorb" || string(argv[i]) == "--outorbits") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rep" || string(argv[i]) == "-repfile") {
      if(i+1 < argc) {
	//There is still something to read;
	repfile=argv[++i];
	i++;
      }
      else {
	cerr << "Representative state vector file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  
  cout.precision(17);  
  cout << "randum number seed string " << seedstring << "\n";
  cout << "input configuration file " << configfile << "\n";
  cout << "input observation file " << obsfile << "\n";
  cout << "maxnum = " << maxnum << "\n";
  cout << "repnum = " << repnum << "\n";
  cout << "output file " << outfile << "\n";

  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, verbose);
  if(status!=0) {
    cerr << "ERROR reading observatory code file " << obscodefile << "\n";
    return(1);
  }
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";
  if(DEBUG>=2) {
    for(i=0;i<long(observatory_list.size());i++) {
      cout << observatory_list[i].obscode << " " << observatory_list[i].obslon << " " << observatory_list[i].plxcos << " " << observatory_list[i].plxsin << "\n";
    }
  }
  
  // Read input observation file.
  detvec={};
  status=read_hldet_file(obsfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << obsfile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  obsnum = detvec.size();
  cout << "Read " << obsnum << " data lines from observation file " << obsfile << "\n";

  // Load the MJD, RA, and Dec vectors needed by arctrace02.
  obsMJD = {};
  obsRA = obsDec = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD.push_back(detvec[obsct].MJD);
    obsRA.push_back(detvec[obsct].RA);
    obsDec.push_back(detvec[obsct].Dec);
    sigastrom.push_back(detvec[obsct].sig_across);
  }
  // Calculate the observer's heliocentric position at the time of each observation.
  observer_barypos = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    status = obscode_lookup(observatory_list,detvec[obsct].obscode,obslon,plxcos,plxsin);
    if(status>0) {
      cerr << "ERROR: obscode_lookup failed for observatory code " << detvec[obsct].obscode << "\n";
      return(3);
    }
    // Calculate observer's exact barycentric position and velocity.
    observer_barycoords01LD(obsMJD[obsct], 5, (long double)obslon, (long double)plxcos, (long double)plxsin, planetmjd, Earthpos, obspos);
    observer_barypos.push_back(obspos);
  }
  cout << obsMJD.size() << " " << obsMJD[obsMJD.size()-1] << " " << obsRA[obsMJD.size()-1] << " " << obsDec[obsMJD.size()-1] << "\n";
  if(obsnum<2 || long(obsRA.size())!=obsnum || long(obsDec.size())!=obsnum || long(observer_barypos.size())!=obsnum) {
    cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << obsnum << " " << obsRA.size() << " " << obsDec.size() << " " << observer_barypos.size() << "\n";
    return(1);
  }

  planetfile_refpoint=-99;
  j=0;
  while(j<long(planetmjd.size()) && planetfile_refpoint==-99) {
    j++;
    if(fabs(planetmjd[j]-mjdstart) <= IMAGETIMETOL/SOLARDAY) planetfile_refpoint=j;
  }
  if(planetfile_refpoint>=0) {
    cout << "Input mjd " << mjdstart << " corresponds to planet file entry number " << planetfile_refpoint << "\n";
  } else {
    cerr << "ERROR: input mjd " << mjdstart << " does not match any entry in the planet file\n";
    return(1);
  }

  chisq = tortoisechi03(polyorder,planetnum,planetmjd,planetmasses,planetpos,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,startpos,startvel,planetfile_refpoint,fitRA,fitDec,resid);
  // Load starting position and velocity into state vector.
  startstate = posvel2statevec02(startpos, startvel, timescale);
  // Initialize eigenvector matrix to the identity matrix
  for(i=0;i<6;i++) {
    for(j=0;j<6;j++) {
      if(i==j) eigenvecs[i][j] = 1.0l;
      else eigenvecs[i][j] = 0.0l;
    }
  }
  
  cout << "Initial chi-square value is " << chisq << "\n";
  ofstream outstream1 {outfile};
  changepar=0;
  for(i=0;i<6;i++) normvec[i] = keepvec[i] = 0.0;
  resetct=covarct=0;
  accepted_states={};
  accepted_all={};
  statechi={};
  for(itct=0;itct<maxnum;itct++) {
    tempstate = startstate;
    ldval = gaussian_deviate_mt(generator);
    for(i=0;i<6;i++) {
      tempstate[i] += eigenvecs[i][changepar]*stepvec[changepar]*ldval;
    }
    statevec2posvel02(tempstate, testpos, testvel, timescale);
    newchi = tortoisechi03(polyorder,planetnum,planetmjd,planetmasses,planetpos,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,testpos,testvel,planetfile_refpoint,fitRA,fitDec,resid);
    if(itct%PRINTNUM==0) cout << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << testpos.x << " " << testpos.y << " " << testpos.z << " "  << fixed << setprecision(10) << testvel.x << " " << testvel.y << " " << testvel.z << " " << newchi << " ";
    if(newchi<=chisq) {
      // The fit has become better: accept the point
      chisq = newchi;
      startstate = tempstate;
      keepvec[changepar]+=1.0;
      normvec[changepar]+=1.0;
      accepted_states.push_back(tempstate);
      accepted_all.push_back(tempstate);
      statechi.push_back(newchi);
      covarct++;
      outstream1 << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << testpos.x << " " << testpos.y << " " << testpos.z << " "  << fixed << setprecision(10) << testvel.x << " " << testvel.y << " " << testvel.z << " " << newchi << "\n";
      if(itct%PRINTNUM==0) cout << "ACCEPTED" << " itct " << itct << "\n";
      if(newchi<=bestchi) {
	bestchi = newchi;
	beststate = startstate;
	bestRA = fitRA;
	bestDec = fitDec;
	bestresid = resid;
	bestephRA = RA1;
	bestephDec = Dec1;
	bestephRAvel = RAvel;
	bestephDecvel = Decvel;
	ldval = 0.0;
	for(obsct=0;obsct<obsnum;obsct++) ldval += bestresid[obsct]*bestresid[obsct];
	astromrms = sqrt(ldval/double(obsnum));
      }
    } else {
      double testpar = unitvar(generator);
      if(testpar<exp(0.5*(chisq-newchi))) {
	// The fit has become only slightly worse: accept the point
	chisq = newchi;
	startstate = tempstate;
	keepvec[changepar]+=1.0;
	normvec[changepar]+=1.0;
	accepted_states.push_back(tempstate);
	accepted_all.push_back(tempstate);
	statechi.push_back(newchi);
	covarct++;
	outstream1 << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << testpos.x << " " << testpos.y << " " << testpos.z << " "  << fixed << setprecision(10) << testvel.x << " " << testvel.y << " " << testvel.z << " " << newchi << "\n";
	if(itct%PRINTNUM==0) cout << "ACCEPTED" << " itct " << itct << "\n";
      } else {
	// Reject the point
	normvec[changepar]+=1.0;
	if(itct%PRINTNUM==0) cout << "REJECTED" << " itct " << itct << "\n";
      }
    }
    resetct++;
    if(resetct>=RESETNUM) {
      // Test to see if the step sizes are appropriate
      for(i=0;i<6;i++) {
	if(normvec[i]>0) keepvec[i]/=normvec[i];
	else {
	  cerr << "ERROR: zero normalization for changepar " << i << "\n";
	  return(1);
	}
      }
      for(i=0;i<6;i++) {
	if(keepvec[i]<MINACCEPT) {
	  // Step size for parameter i is too coarse. Make it smaller.
	  stepvec[i]/=2.0;
	} else if(keepvec[i]>MAXACCEPT) {
	  // Step size for parameter i is too fine. Make it larger.
	  stepvec[i]*=2.0;
	}
      }
      for(i=0;i<6;i++) {
	cout << "Parameter " << i << ", keepvec was " << keepvec[i] << ", step size " << stepvec[i] << "\n";
      }
      // Reset the reset count, and zero the acceptance-tracking vectors
      resetct=0;
      for(i=0;i<6;i++) normvec[i] = keepvec[i] = 0.0;
    }
    if(covarct>=covarnum) {
      // Calculate the covariance matrix
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) {
	  covariance_matrix[i][j] = 0.0l;
	  for(k=0;k<long(accepted_states.size());k++) {
	    covariance_matrix[i][j] += (accepted_states[k][i]-beststate[i])*(accepted_states[k][j]-beststate[j]);
	  }
	  covariance_matrix[i][j]/=(long double)(accepted_states.size());
	}
      }
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) {
	  cout << covariance_matrix[i][j] << " ";
	}
	cout << "\n";
      }
      // Calculate eigenvectors of the covariance matrix
      status = eigensolve02(covariance_matrix,eigenvecs,eigenvals,EIGENOFFMAX,EIGENITMAX);
      if(status!=0) {
	cerr << "ERROR: eigensolve01 returned error value " << status << "\n";
      }
      for(i=0;i<6;i++) {
	cout << "Eigenvalue " << i << " is " << eigenvals[i] << ", with eigenvector " << "\n";
	for(j=0;j<6;j++) cout << eigenvecs[j][i] << " ";
	cout << "\n";
      }
      covarct=0;
      accepted_states={};
    }
    // Cycle changepar	  
    if(changepar<5 && changepar>=0) changepar++;
    else changepar=0;
  }
  outstream1.close();

  statevec2posvel02(beststate, outpos, outvel, timescale);
  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromrms << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  cout << "State vectors correspond to reference point " << planetfile_refpoint << " in the input planet files\n";

  // Loop over all accepted states to find representative ones that span the space
  outstream1.open(repfile);
  statedist={};
  representative_states={};
  representative_states.push_back(beststate);
  statevec2posvel02(beststate, outpos, outvel, timescale);
  maxdist=0.0l;
  cout << fixed << setprecision(10) << "Representative statevec at MJD: " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << " chisq = " << bestchi << " dist = " << maxdist << "\n";
  outstream1 << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << " " << bestchi << " " << "\n";
  // Calculate the parameter-space distance from each accepted state to the best state
  for(i=0;i<long(accepted_all.size());i++) {
    statedist.push_back(sqrt(LDSQUARE(accepted_all[i][0]-beststate[0]) + LDSQUARE(accepted_all[i][1]-beststate[1]) + LDSQUARE(accepted_all[i][2]-beststate[2])));
  }
  if(statedist.size() != accepted_all.size()) {
    cerr << "ERROR: statedist vector does not have the same length as accepted_all!\n";
    cerr << "Difference is " << statedist.size() << " vs. " << accepted_all.size() << "\n";
    return(3);
  }
  repct=1;
  while(repct<repnum) {
    // Identify the state that is furthest from the nearest representative state
    maxdist = statedist[0];
    distantstate = 0;
    for(i=0;i<long(statedist.size());i++) {
      if(statedist[i] > maxdist) {
	maxdist = statedist[i];
	distantstate = i;
      }
    }
    // Load this newly-identified most isolated state to representative_states
    onestate = accepted_all[distantstate];
    representative_states.push_back(onestate);
    statevec2posvel02(onestate, testpos, testvel, timescale);
    cout << fixed << setprecision(10) << "Representative statevec at MJD: " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << testpos.x << " " << testpos.y << " " << testpos.z << " "  << fixed << setprecision(10) << testvel.x << " " << testvel.y << " " << testvel.z << " chisq = " << statechi[distantstate] << " dist = " << maxdist << "\n";
    outstream1 << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << testpos.x << " " << testpos.y << " " << testpos.z << " "  << fixed << setprecision(10) << testvel.x << " " << testvel.y << " " << testvel.z << " " << statechi[distantstate] << "\n";

    // Re-calculate all the distances
    for(i=0;i<long(accepted_all.size());i++) {
      ldval = sqrt(LDSQUARE(accepted_all[i][0]-onestate[0]) + LDSQUARE(accepted_all[i][1]-onestate[1]) + LDSQUARE(accepted_all[i][2]-onestate[2]));
      if(ldval<statedist[i]) statedist[i] = ldval;
    }
    repct++;
  }
  outstream1.close();

  return(0);
}

vector <long double> posvel2statevec01(point3LD pos, point3LD vel)
{
  vector <long double> outvec;
  outvec.push_back(pos.x);
  outvec.push_back(pos.y);
  outvec.push_back(pos.z);
  outvec.push_back(vel.x);
  outvec.push_back(vel.y);
  outvec.push_back(vel.z);
  return(outvec);
}

int statevec2posvel01(const vector <long double> &invec, point3LD &pos, point3LD &vel)
{
  if(long(invec.size())<6) {
    cerr << "ERROR: statevec2posvel01 given input vector of insufficient length\n";
    return(1);
  } else if(long(invec.size())>6) {
    cerr << "WARNING: statevec2posvel01 supplied input vector with more than 6 entries (total is " << invec.size() << ")\n";
  }
  pos.x = invec[0];
  pos.y = invec[1];
  pos.z = invec[2];
  vel.x = invec[3];
  vel.y = invec[4];
  vel.z = invec[5];
  return(0);
}

vector <long double> posvel2statevec02(point3LD pos, point3LD vel, long double timescale)
{
  vector <long double> outvec;
  outvec.push_back(pos.x);
  outvec.push_back(pos.y);
  outvec.push_back(pos.z);
  outvec.push_back(vel.x*timescale);
  outvec.push_back(vel.y*timescale);
  outvec.push_back(vel.z*timescale);
  return(outvec);
}

int statevec2posvel02(const vector <long double> &invec, point3LD &pos, point3LD &vel, long double timescale)
{
  if(long(invec.size())<6) {
    cerr << "ERROR: statevec2posvel02 given input vector of insufficient length\n";
    return(1);
  } else if(long(invec.size())>6) {
    cerr << "WARNING: statevec2posvel02 supplied input vector with more than 6 entries (total is " << invec.size() << ")\n";
  }
  pos.x = invec[0];
  pos.y = invec[1];
  pos.z = invec[2];
  vel.x = invec[3]/timescale;
  vel.y = invec[4]/timescale;
  vel.z = invec[5]/timescale;
  return(0);
}
