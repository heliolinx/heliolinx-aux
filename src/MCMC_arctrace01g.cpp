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
#define DEBUG 1
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
  cerr << "Usage: MCMC_arctrace01g -cfg configfile -ranseed random_number_seed -observations obsfile -mjd starting_mjd -statevec x y z vx vy vz -maxnum maxnum -repnum repnum -obscode obscodefile -outfile outfile -repfile repfile -verbose verbosity\n";
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
  itct = resetct = planetnum = planetct = pctSun = pctEarth = 0;
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
  long repstate = 0;
  long repnum = REPNUM;
  long double main_sample_step=0.0l;
  long double eigenrange0 = 0.0l;
  long double eigenrange1 = 0.0l;
  long double eigenrange2 = 0.0l;
  long double eigenrange3 = 0.0l;
  long double eigenrange4 = 0.0l;
  long double eigenrange5 = 0.0l;

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
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--VERB" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--verbose") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "keyword for verbosity supplied with no corresponding argument";
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
    if(itct%PRINTNUM==0 || verbose>0) cout << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << testpos.x << " " << testpos.y << " " << testpos.z << " "  << fixed << setprecision(10) << testvel.x << " " << testvel.y << " " << testvel.z << " " << newchi << " ";
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
      if(itct%PRINTNUM==0 || verbose>0) cout << "ACCEPTED" << " itct " << itct << "\n";
      if(newchi<=bestchi) {
	bestchi = newchi;
	beststate = startstate;
	bestRA = fitRA;
	bestDec = fitDec;
	bestresid = resid;
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
	if(itct%PRINTNUM==0 || verbose>0) cout << "ACCEPTED" << " itct " << itct << "\n";
      } else {
	// Reject the point
	normvec[changepar]+=1.0;
	if(itct%PRINTNUM==0 || verbose>0) cout << "REJECTED" << " itct " << itct << "\n";
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

  // Calculate the covariance matrix of the entire accepted sample set
  for(i=0;i<6;i++) {
    for(j=0;j<6;j++) {
      covariance_matrix[i][j] = 0.0l;
      for(k=0;k<long(accepted_all.size());k++) {
	covariance_matrix[i][j] += (accepted_all[k][i]-beststate[i])*(accepted_all[k][j]-beststate[j]);
      }
      covariance_matrix[i][j]/=(long double)(accepted_all.size());
    }
  }
  cout << "Covariance matrix over all accepted states:\n";
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

  // Load best state to the ouput matrix
  outstream1.open(repfile);
  statedist={};
  representative_states={};
  representative_states.push_back(beststate);
  
  // Loop over all accepted states to find representative ones that span the space
  // Start with eigenvector 0, which has the largest eigenvalue
  long double ldelem = 0.0l;
  long index = 0;
  ldouble_index ldi = ldouble_index(ldelem,index);
  vector <ldouble_index> ldivec0 ={};
  
  // Take the dot-product of eigenvector 0 with every accepted state,
  // and store in the vector ldivec0.
  for(i=0;i<long(accepted_all.size());i++) {
    // Take the dot product of (accepted_all[i]-beststate) with eigenvector 0
    ldelem = 0.0l;
    for(j=0;j<6;j++) ldelem += (accepted_all[i][j]-beststate[j])*eigenvecs[j][0];
    index = i;
    ldi = ldouble_index(ldelem,index);
    ldivec0.push_back(ldi);
  }

  // Sort ldivec0 by increasing value of the dot-product
  sort(ldivec0.begin(),ldivec0.end(),lower_ldouble_index());
  // Calculate the full range in the dot-product with this eigenvector
  eigenrange0 = ldivec0[ldivec0.size()-1].ldelem - ldivec0[0].ldelem;
  main_sample_step = eigenrange0/double(repnum);
  long eigenrepnum0 = repnum;
  long eigenrepct0 = 0;
  if(DEBUG>0) cout << "eigenrange0 = " << eigenrange0 << ", main_sample_step = " << main_sample_step << "\n";
  // Split ldivec0 into eigenrepnum0 equally spaced parts
  long i0=0;
  for(eigenrepct0=1;eigenrepct0<=eigenrepnum0;eigenrepct0++) {
    vector <ldouble_index> ldivec1 ={};
    while(i0<long(ldivec0.size()) && ldivec0[i0].ldelem < ldivec0[0].ldelem + eigenrange0*double(eigenrepct0)/double(eigenrepnum0)) {
      // Take the dot-product of this state relative to eigenvector 1
      k = ldivec0[i0].index;
      ldelem = 0.0l;
      for(j=0;j<6;j++) ldelem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][1];
      index = k;
      ldi = ldouble_index(ldelem,index);
      ldivec1.push_back(ldi);
      i0++;
    }
    // We are working with a subset of points, split by the direction defined by eigenvector 0
    if(DEBUG>0) cout << "Loaded ldivec1 with " << ldivec1.size() << " points\n";
    eigenrange1 = 0.0l;
    // Sort ldivec1 by increasing value of the dot-product with eigenvector 1
    if(ldivec1.size()>0) {
      sort(ldivec1.begin(),ldivec1.end(),lower_ldouble_index());
      eigenrange1 = ldivec1[ldivec1.size()-1].ldelem - ldivec1[0].ldelem;
      if(DEBUG>0) cout << "eigenrepct0 = " << eigenrepct0 << ", eigenrange1 = " << eigenrange1 << ", main_sample_step = " << main_sample_step << "\n";
    }
    if(ldivec1.size()>0 && eigenrange1>main_sample_step) {
      // Must divide this subset further
      long eigenrepnum1 = eigenrange1/main_sample_step + 1.0;
      long eigenrepct1 = 0;
      long i1=0;
      if(DEBUG>0) cout << "eigenrepnum1 = " << eigenrepnum1 << "\n";
      for(eigenrepct1=1;eigenrepct1<=eigenrepnum1;eigenrepct1++) {
	vector <ldouble_index> ldivec2 ={};
	while(i1<long(ldivec1.size()) && ldivec1[i1].ldelem < ldivec1[0].ldelem + eigenrange1*double(eigenrepct1)/double(eigenrepnum1)) {
	  // Take the dot-product of this state relative to eigenvector 2
	  k = ldivec1[i1].index;
	  ldelem = 0.0l;
	  for(j=0;j<6;j++) ldelem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][2];
	  index = k;
	  ldi = ldouble_index(ldelem,index);
	  ldivec2.push_back(ldi);
	  i1++;
	}
	// We are working with a subset of points, split by the direction defined by eigenvector 1
	// Sort ldivec2 by increasing value of the dot-product with eigenvector 2
	eigenrange2 = 0.0l;
	if(ldivec2.size()>0) {
	  sort(ldivec2.begin(),ldivec2.end(),lower_ldouble_index());
	  eigenrange2 = ldivec2[ldivec2.size()-1].ldelem - ldivec2[0].ldelem;
	  if(DEBUG>0) cout << "eigenrepct1 = " << eigenrepct1 << ", eigenrange2 = " << eigenrange2 << ", main_sample_step = " << main_sample_step << "\n";
	}
	if(ldivec2.size()>0 && eigenrange2>main_sample_step) {
	  // Must divide this subset further
	  long eigenrepnum2 = eigenrange2/main_sample_step + 1.0;
	  long eigenrepct2 = 0;
	  long i2=0;
	  if(DEBUG>0) cout << "eigenrepnum2 = " << eigenrepnum2 << "\n";
	  for(eigenrepct2=1;eigenrepct2<=eigenrepnum2;eigenrepct2++) {
	    vector <ldouble_index> ldivec3 ={};
	    while(i2<long(ldivec2.size()) && ldivec2[i2].ldelem < ldivec2[0].ldelem + eigenrange2*double(eigenrepct2)/double(eigenrepnum2)) {
	      // Take the dot-product of this state relative to eigenvector 3
	      k = ldivec2[i2].index;
	      ldelem = 0.0l;
	      for(j=0;j<6;j++) ldelem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][3];
	      index = k;
	      ldi = ldouble_index(ldelem,index);
	      ldivec3.push_back(ldi);
	      i2++;
	    }
	    // We are working with a subset of points, split by the direction defined by eigenvector 2
	    // Sort ldivec3 by increasing value of the dot-product with eigenvector 3
	    eigenrange3 = 0.0l;
	    if(ldivec3.size()>0) {
	      sort(ldivec3.begin(),ldivec3.end(),lower_ldouble_index());
	      eigenrange3 = ldivec3[ldivec3.size()-1].ldelem - ldivec3[0].ldelem;
	      if(DEBUG>0) cout << "eigenrepct2 = " << eigenrepct2 << ", eigenrange3 = " << eigenrange3 << ", main_sample_step = " << main_sample_step << "\n";
	    }
	    if(ldivec3.size()>0 && eigenrange3>main_sample_step) {
	      // Must divide this subset further
	      long eigenrepnum3 = eigenrange3/main_sample_step + 1.0;
	      long eigenrepct3 = 0;
	      long i3=0;
	      if(DEBUG>0) cout << "eigenrepnum3 = " << eigenrepnum3 << "\n";
	      for(eigenrepct3=1;eigenrepct3<=eigenrepnum3;eigenrepct3++) {
		vector <ldouble_index> ldivec4 ={};
		while(i3<long(ldivec3.size()) && ldivec3[i3].ldelem < ldivec3[0].ldelem + eigenrange3*double(eigenrepct3)/double(eigenrepnum3)) {
		  // Take the dot-product of this state relative to eigenvector 4
		  k = ldivec2[i3].index;
		  ldelem = 0.0l;
		  for(j=0;j<6;j++) ldelem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][4];
		  index = k;
		  ldi = ldouble_index(ldelem,index);
		  ldivec4.push_back(ldi);
		  i3++;
		}
		// We are working with a subset of points, split by the direction defined by eigenvector 3
		eigenrange4 = 0.0l;
		// Sort ldivec4 by increasing value of the dot-product with eigenvector 4
		if(ldivec4.size()>0) {
		  sort(ldivec4.begin(),ldivec4.end(),lower_ldouble_index());
		  eigenrange4 = ldivec4[ldivec4.size()-1].ldelem - ldivec4[0].ldelem;
		  if(DEBUG>0) cout << "eigenrepct3 = " << eigenrepct3 << ", eigenrange4 = " << eigenrange4 << ", main_sample_step = " << main_sample_step << "\n";
		}
		if(ldivec4.size()>0 && eigenrange4>main_sample_step) {
		  // Must divide this subset further
		  long eigenrepnum4 = eigenrange4/main_sample_step + 1.0;
		  long eigenrepct4 = 0;
		  long i4=0;
		  if(DEBUG>0) cout << "eigenrepnum4 = " << eigenrepnum4 << "\n";
		  for(eigenrepct4=1;eigenrepct4<=eigenrepnum4;eigenrepct4++) {
		    vector <ldouble_index> ldivec5 ={};
		    while(i4<long(ldivec4.size()) && ldivec4[i4].ldelem < ldivec4[0].ldelem + eigenrange4*double(eigenrepct4)/double(eigenrepnum4)) {
		      // Take the dot-product of this state relative to eigenvector 5
		      k = ldivec4[i4].index;
		      ldelem = 0.0l;
		      for(j=0;j<6;j++) ldelem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][5];
		      index = k;
		      ldi = ldouble_index(ldelem,index);
		      ldivec5.push_back(ldi);
		      i4++;
		    }
		    // We are working with a subset of points, split by the direction defined by eigenvector 4
		    eigenrange5 = 0.0l;
		    // Sort ldivec5 by increasing value of the dot-product with eigenvector 5
		    if(ldivec5.size()>0) {
		      sort(ldivec5.begin(),ldivec5.end(),lower_ldouble_index());
		      eigenrange5 = ldivec5[ldivec5.size()-1].ldelem - ldivec5[0].ldelem;
		      if(DEBUG>0) cout << "eigenrepct4 = " << eigenrepct4 << ", eigenrange5 = " << eigenrange5 << ", main_sample_step = " << main_sample_step << "\n";
		    }
		    if(ldivec5.size()>0 && eigenrange5>main_sample_step) {
		      // Must divide this subset further
		      long eigenrepnum5 = eigenrange5/main_sample_step + 1.0;
		      long eigenrepct5 = 0;
		      long i5=0;
		      if(DEBUG>0) cout << "eigenrepnum5 = " << eigenrepnum5 << "\n";
		      for(eigenrepct5=1;eigenrepct5<=eigenrepnum5;eigenrepct5++) {
			while(i5<long(ldivec5.size()) && ldivec5[i5].ldelem < ldivec5[0].ldelem + eigenrange5*(double(eigenrepct5)-0.5l)/double(eigenrepnum5)) i5++;
			if(i5<long(ldivec5.size())) {
			  // i5 should be the middle of a sample step
			  repstate = ldivec5[i5].index;
			  representative_states.push_back(accepted_all[repstate]);
			}
		      }
		      // End case in which it was necessary to subdivide along the dimension defined by eigenvector 5
		    } else if(ldivec5.size()>0) {
		      if(DEBUG>0) cout << "eigenrange5 was small enough that no further division was necessary\n";
		      // Just output the median point.
		      k=ldivec5.size()/2;
		      repstate = ldivec5[k].index;
		      representative_states.push_back(accepted_all[repstate]);
		    } else {
		      cout << "Split on eigenvector 5 had zero points\n";
		    }
		  }
		  // End case in which it was necessary to subdivide along the dimension defined by eigenvector 4
		} else if(ldivec4.size()>0) {
		  if(DEBUG>0) cout << "eigenrange4 was small enough that no further division was necessary\n";
		  // Just output the median point.
		  k=ldivec4.size()/2;
		  repstate = ldivec4[k].index;
		  representative_states.push_back(accepted_all[repstate]);
		} else {
		  cout << "Split on eigenvector 4 had zero points\n";
		}
	      }
	      // End case in which it was necessary to subdivide along the dimension defined by eigenvector 3
	    } else if(ldivec3.size()>0) {
	      if(DEBUG>0) cout << "eigenrange3 was small enough that no further division was necessary\n";
	      // Just output the median point.
	      k=ldivec3.size()/2;
	      repstate = ldivec3[k].index;
	      representative_states.push_back(accepted_all[repstate]);
	    } else {
	      cout << "Split on eigenvector 3 had zero points\n";
	    }
	  }
	  // End case in which it was necessary to subdivide along the dimension defined by eigenvector 2
	} else if(ldivec2.size()>0) {
	  if(DEBUG>0) cout << "eigenrange2 was small enough that no further division was necessary\n";
	  // Just output the median point.
	  k=ldivec2.size()/2;
	  repstate = ldivec2[k].index;
	  representative_states.push_back(accepted_all[repstate]);
	} else {
	  cout << "Split on eigenvector 2 had zero points\n";
	}
      }
      // End case in which it was necessary to subdivide along the dimension defined by eigenvector 1
    } else if(ldivec1.size()>0) {
      if(DEBUG>0) cout << "eigenrange1 was small enough that no further division was necessary\n";
      // Just output the median point.
      k=ldivec1.size()/2;
      if(DEBUG>0) cout << "size = " << ldivec1.size() << ", k = " << k << "\n";
      repstate = ldivec1[k].index;
      if(DEBUG>0) cout << "repstate = " << repstate << ", allsize = " << accepted_all.size() << "\n";
      representative_states.push_back(accepted_all[repstate]);
    } else {
      cout << "Split had zero points\n";
    }
  }
  // Write all of the representative states to the output file
  cout << "Writing " << representative_states.size() << " lines to output file " << repfile << "\n";
  for(i=0;i<long(representative_states.size());i++) {
    onestate = representative_states[i];
    statevec2posvel02(onestate, outpos, outvel, timescale);
    cout << fixed << setprecision(10) << "Representative statevec at MJD: " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
    outstream1 << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << " " "\n";
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
