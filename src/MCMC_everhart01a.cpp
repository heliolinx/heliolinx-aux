// October 21, 2025: MCMC_everhart01a.cpp
// Like MCMC_arctrace01g.cpp, but integrates orbits with the newer functions
// based on Everhart integrators.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define DEBUG 0
#define HNUM 4
#define RESETNUM 2000
#define MINACCEPT 0.2
#define MAXACCEPT 0.5
#define POSSTEP 1000.0
#define PRINTNUM 1000
#define COVARNUM 10000
#define VELTIMESCALE 1.0e7l
#define EIGENOFFMAX 1e-50
#define EIGENITMAX 1000000
#define REPNUM 10

vector <double> scale_statevel(const vector <double> &statevec, double timescale);
vector <double> unscale_statevel(const vector <double> &statevec, double timescale);

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain observations in heliolinc's hldet format.
static void show_usage()
{
  cerr << "Usage: MCMC_everhart01a -cfg configfile -ranseed random_number_seed -observations obsfile -mjdref mjdref -maxnum maxnum -repnum repnum -obscode obscodefile -rmsfloor min astrometric uncertainty -outfile outfile -repfile repfile -verbose verbosity\n";
}

int main(int argc, char *argv[])
{
  long maxnum=10000;
  long i,j,k;
  i=j=k=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string obsfile,obscodefile;
  vector <hldet> detvec = {};
  vector <observatory> observatory_list = {};
  string imfile;
  string Earthfile,Sunfile,planetfile;
  string spaceobsfile;
  string outfile = "";
  string repfile = "repjunk01a.txt";
  double mjdref = 0.0L;
  int configread=0;
  int planetnum=0;
  int planetct=0;
  double ldval=0.0;
  vector <vector <double>> Earth_statevecs;
  vector <vector <double>> Sun_statevecs;
  vector <double> Earth_mjd;
  vector <double> Sun_mjd;
  vector <vector <double>> statevecs;
  vector <double> planetmasses;
  vector <double> mjdtest;
  vector <double> planet_backward_mjd;
  vector <vector <double>> planet_backward_statevecs;
  vector <double> planet_forward_mjd;
  vector <vector <double>> planet_forward_statevecs;
  vector <double> starting_statevec;
  vector <double> tempstate;
  vector <double> startstate;
  vector <double> beststate;
  vector <double> realstate;
  make_dvec(6, starting_statevec);
  vector <double> targMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> bestRA;
  vector <double> bestDec;
  vector <double> fitresid;
  vector <double> obsMJD;
  vector <double> obsTDB;
  vector <double> sigastrom;
  string lnfromfile;
  vector <string> linestringvec;
  ofstream outstream1;
  double timestep = 5.0;
  int hnum = HNUM;
  vector <double> hspace;
  int verbose = 0;
  vector <double> obsstate;
  vector <double> sunstate;
  vector <vector <double>> observer_statevecs;
  vector <vector <double>> observer_heliostate;
  long obsnum,obsct,refpoint;
  double obslon,plxcos,plxsin;
  double astromrms=0.0;
  double timescale = VELTIMESCALE;

  double chisq=0.0l;
  double newchi=0.0l;
  double bestchi = LARGERR;
  vector <double> onestate;
  vector <vector <double>> accepted_states;
  vector <vector <double>> representative_states;
  vector <vector <double>> accepted_all;
  vector <vector <double>> covariance_matrix;
  vector <vector <double>> eigenvecs;
  vector <double> eigenvals;
  vector <double> statedist;
  vector <double> statechi;
  long repstate = 0;
  long repnum = REPNUM;
  make_dvec(6,onestate);
  make_dmat(6,6,covariance_matrix);
  make_dmat(6,6,eigenvecs);
  string seedstring;
  int changepar=0;
  vector <double> stepvec;
  vector <double> keepvec;
  vector <double> normvec;
  make_dvec(6, stepvec);
  make_dvec(6, keepvec);
  make_dvec(6, normvec);
  stepvec[0] = stepvec[1] = stepvec[2] = POSSTEP;
  stepvec[3] = stepvec[4] = stepvec[5] = POSSTEP;
  long covarnum = COVARNUM;
  long covarct = 0;
  long resetct = 0;
  long itct = 0;
  double main_sample_step=0.0l;
  double eigenrange0 = 0.0l;
  double eigenrange1 = 0.0l;
  double eigenrange2 = 0.0l;
  double eigenrange3 = 0.0l;
  double eigenrange4 = 0.0l;
  double eigenrange5 = 0.0l;
  double rmsfloor=0.0;

  if(argc<11) {
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
	// Read HNUM for the Everhart integrator (that is, number of substeps per big step).
	status=readconfigint(instream1,&hnum);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&hnum);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "HNUM (number of substeps per big step) for Everhart integrator read as " << hnum << "\n";
	// Read the timestep for the Everhart integrator in days
	status=readconfigd(instream1,&timestep);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&timestep);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Timestep for Everhart integrator (big step in days) read as " << timestep << "\n";
	// Read the evenly-sampled ephemeris file for the sun
	status=readconfigstring(instream1,Sunfile);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigstring(instream1,Sunfile);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Evenly-sampled ephemeris file for the sun is named " << Sunfile << "\n";
	mjdtest={};
	statevecs={};
	read_horizons_statevec(Sunfile,Sun_mjd,Sun_statevecs);
	if(Sun_statevecs.size()!=Sun_mjd.size()) {
	  cerr << "ERROR: number of statevecs " << Sun_statevecs.size() << " for the sun did not match number of mjds " << Sun_mjd.size() << "\n";
	  return(1);
	}
      	// Read the evenly-sampled ephemeris file for Earth
	status=readconfigstring(instream1,Earthfile);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigstring(instream1,Earthfile);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Evenly-sampled ephemeris file for Earth is named " << Earthfile << "\n";
	mjdtest={};
	statevecs={};
	read_horizons_statevec(Earthfile,Earth_mjd,Earth_statevecs);
	if(Earth_statevecs.size()!=Earth_mjd.size()) {
	  cerr << "ERROR: number of statevecs " << Earth_statevecs.size() << " for Earth did not match number of mjds " << Earth_mjd.size() << "\n";
	  return(1);
	}
	if(Earth_statevecs.size()!=Sun_statevecs.size()) {
	  cerr << "ERROR: number of statevecs " << Earth_statevecs.size() << " for Earth did not match number for the sun " << Sun_statevecs.size() << "\n";
	  return(2);
	}
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
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the mass for planet number planetct
	  status=readconfigd(instream1,&ldval);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigd(instream1,&ldval);
	  }
	  planetmasses.push_back(ldval);
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else {
	    cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
	  }
	}
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the backward ephemeris file for planet number planetct
	  status=readconfigstring(instream1,planetfile);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigstring(instream1,planetfile);
	  }
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "Backward ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
	  mjdtest={};
	  statevecs={};
	  read_t6text_statevec(planetfile,mjdtest,statevecs);
	  if(planetct==0) {
	    planet_backward_mjd=mjdtest;
	    planet_backward_statevecs = statevecs;
	  }
	  else {
	    if(statevecs.size()!=planet_backward_mjd.size()) {
	      cerr << "ERROR: number of backward statevecs " << statevecs.size() << " for planet " << planetct << " did not match number of mjds " << planet_backward_mjd.size() << "\n";
	      return(1);
	    }
	    if(statevecs.size()!=planet_backward_statevecs.size()) {
	      cerr << "ERROR: number of backward statevecs " << statevecs.size() << " for planet " << planetct << " did not match master array " << planet_backward_statevecs.size() << "\n";
	      return(2);
	    }
	    for(j=0;j<long(planet_backward_mjd.size());j++) {
	      if(mjdtest[j]!=planet_backward_mjd[j]) {
		cout << "ERROR: backward time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(3);
	      }
	    }
	    // Add state vectors for the latest planet to the master array
	    for(j=0;j<long(planet_backward_statevecs.size());j++) {
	      for(k=0;k<6;k++) planet_backward_statevecs[j].push_back(statevecs[j][k]);
	    }
	  }
	  // Read the forward ephemeris file for planet number planetct
	  status=readconfigstring(instream1,planetfile);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigstring(instream1,planetfile);
	  }
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "Forward ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
	  mjdtest={};
	  statevecs={};
	  read_t6text_statevec(planetfile,mjdtest,statevecs);
	  if(planetct==0) {
	    planet_forward_mjd=mjdtest;
	    planet_forward_statevecs = statevecs;
	  }
	  else {
	    if(statevecs.size()!=planet_forward_mjd.size()) {
	      cerr << "ERROR: number of forward statevecs " << statevecs.size() << " for planet " << planetct << " did not match number of mjds " << planet_forward_mjd.size() << "\n";
	      return(1);
	    }
	    if(statevecs.size()!=planet_forward_statevecs.size()) {
	      cerr << "ERROR: number of forward statevecs " << statevecs.size() << " for planet " << planetct << " did not match master array " << planet_forward_statevecs.size() << "\n";
	      return(2);
	    }
	    for(j=0;j<long(planet_forward_mjd.size());j++) {
	      if(mjdtest[j]!=planet_forward_mjd[j]) {
		cout << "ERROR: forward time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(3);
	      }
	    }
	    // Add state vectors for the latest planet to the master array
	    for(j=0;j<long(planet_forward_statevecs.size());j++) {
	      for(k=0;k<6;k++) planet_forward_statevecs[j].push_back(statevecs[j][k]);
	    }
	  }
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
    } else if(string(argv[i]) == "-mjdref" || string(argv[i]) == "-mr" || string(argv[i]) == "-MJDref" || string(argv[i]) == "--mjdref" || string(argv[i]) == "--MJDref" || string(argv[i]) == "--ModifiedJulianDayref" || string(argv[i]) == "--modifiedjuliandayref") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD ref keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-statevec" || string(argv[i]) == "-sv" || string(argv[i]) == "-svec" || string(argv[i]) == "--statevector" || string(argv[i]) == "--statevec" || string(argv[i]) == "-statevector" || string(argv[i]) == "--sv") {
      for(j=0;j<6;j++) {
	if(i+1 < argc) {
	  // There is still something to read;
	  starting_statevec[j]=stod(argv[++i]);
	  cout << "read sv entry as " << starting_statevec[j] << "\n";
	} else {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
      }
      i++;
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
    } else if(string(argv[i]) == "-rmsfloor" || string(argv[i]) == "-minsigma" || string(argv[i]) == "-minrms" || string(argv[i]) == "-min_astrom_sigma" || string(argv[i]) == "-min_astrometric_sigma" || string(argv[i]) == "-rms_floor" || string(argv[i]) == "--min_astrometric_sigma" ) {
      if(i+1 < argc) {
	//There is still something to read;
	rmsfloor=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Min astrometric uncertainty keyword supplied with no corresponding argument\n";
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
  if(outfile.size()>0) cout << "output file " << outfile << "\n";
  else cout << "No output file name has been specified, so none will be written\n";
  
  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Load hspace vector based on selected value of hnum
  if(hnum<3 || hnum>8) {
    cerr << "Warning: hnum set to out-of-range value (" << hnum << "), will be reset to default of " << HNUM << "\n";
    hnum=HNUM;
  }
  hspace={};
  if(hnum==3) {
    hspace.push_back(h3_1);
    hspace.push_back(h3_2);
    hspace.push_back(h3_3);
  } else if(hnum==4) {
    hspace.push_back(h4_1);
    hspace.push_back(h4_2);
    hspace.push_back(h4_3);
    hspace.push_back(h4_4);
  } else if(hnum==5) {
    hspace.push_back(h5_1);
    hspace.push_back(h5_2);
    hspace.push_back(h5_3);
    hspace.push_back(h5_4);
    hspace.push_back(h5_5);
  } else if(hnum==6) {
    hspace.push_back(h6_1);
    hspace.push_back(h6_2);
    hspace.push_back(h6_3);
    hspace.push_back(h6_4);
    hspace.push_back(h6_5);
    hspace.push_back(h6_6);
  } else if(hnum==7) {
    hspace.push_back(h7_1);
    hspace.push_back(h7_2);
    hspace.push_back(h7_3);
    hspace.push_back(h7_4);
    hspace.push_back(h7_5);
    hspace.push_back(h7_6);
    hspace.push_back(h7_7);
  } else if(hnum==8) {
    hspace.push_back(h8_1);
    hspace.push_back(h8_2);
    hspace.push_back(h8_3);
    hspace.push_back(h8_4);
    hspace.push_back(h8_5);
    hspace.push_back(h8_6);
    hspace.push_back(h8_7);
    hspace.push_back(h8_8);
  }
  
  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, verbose);
  if(status!=0) {
    cerr << "ERROR reading observatory code file " << obscodefile << "\n";
    return(1);
  }
  
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";
  if(verbose>=2) {
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

  // Load the MJD, RA, and Dec vectors.
  obsMJD = obsTDB = obsRA = obsDec = sigastrom = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD.push_back(detvec[obsct].MJD);
    obsTDB.push_back(detvec[obsct].MJD + TTDELTAT/SOLARDAY);
    obsRA.push_back(detvec[obsct].RA);
    obsDec.push_back(detvec[obsct].Dec);
    sigastrom.push_back(sqrt(DSQUARE(detvec[obsct].sig_across)+DSQUARE(detvec[obsct].sig_along)+DSQUARE(rmsfloor)));
    cout << "sigastrom = " << sigastrom[obsct] << "\n";
  }
  // Check that MJD values are time-ordered.
  for(obsct=1;obsct<obsnum;obsct++) {
    if(obsMJD[obsct]<obsMJD[obsct-1]) {
      cerr << "ERROR: input observations are not properly sorted by time.\n";
      return(2);
    }
    if(verbose>0) cout << "obsct, MJD, RA, Dec: " << obsct << " " << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << "\n";
  }

  // Calculate the observer's heliocentric position at the time of each observation.
  observer_statevecs = observer_heliostate = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    status = obscode_lookup(observatory_list,detvec[obsct].obscode,obslon,plxcos,plxsin);
    if(verbose>0) cout << "Detection " << obsct << " is from obscode " << detvec[obsct].obscode << ", coords " << obslon << " " << plxcos << " " << plxsin << "\n";
    if(status>0) {
      cerr << "ERROR: obscode_lookup failed for observatory code " << detvec[obsct].obscode << "\n";
      return(3);
    }
    // Calculate observer's exact barycentric position and velocity.
    // Note that observer_barystate01 assumes input is UTC, so we use obsMJD here rather than obsTDB
    if(verbose>0) cout << "About to call observer_barystate01 " << Earth_mjd.size() << " " << Earth_statevecs.size() << "\n";
    observer_barystate01(obsMJD[obsct], 5, obslon, plxcos, plxsin, Earth_mjd, Earth_statevecs, obsstate, verbose);
    observer_statevecs.push_back(obsstate);
    // Calculate the sun's position at the same time.
    // Note that planetpos02 assumes input is TDB
    planetpos02(obsTDB[obsct],5,Sun_mjd,Sun_statevecs,sunstate);
    // Convert obsstate from barycentric to heliocentric coords.
    for(k=0;k<6;k++) obsstate[k] -= sunstate[k];
    observer_heliostate.push_back(obsstate);
    if(verbose>0) {
      cout << "Observation " << obsct << ", barystate : " << observer_statevecs[obsct][0] << " " << observer_statevecs[obsct][1] << " " << observer_statevecs[obsct][2] << " " << observer_statevecs[obsct][3] << " " << observer_statevecs[obsct][4] << " " << observer_statevecs[obsct][5] << "\n";
       cout << "heliostate : " << observer_heliostate[obsct][0] << " " << observer_heliostate[obsct][1] << " " << observer_heliostate[obsct][2] << " " << observer_heliostate[obsct][3] << " " << observer_heliostate[obsct][4] << " " << observer_heliostate[obsct][5] << "\n";
    }
  }
  cout << obsMJD.size() << " " << obsMJD[obsMJD.size()-1] << " " << obsRA[obsMJD.size()-1] << " " << obsDec[obsMJD.size()-1] << "\n";
  if(obsnum<2 || long(obsRA.size())!=obsnum || long(obsDec.size())!=obsnum || long(observer_heliostate.size())!=obsnum) {
    cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << obsnum << " " << obsRA.size() << " " << obsDec.size() << " " << observer_heliostate.size() << "\n";
    return(1);
  }
  
  // Make sure that mjdref matches a point in the Everhart-sampled planet files
  refpoint = -99;
  for(j=0;j<long(planet_forward_mjd.size());j++) {
    if(fabs(planet_forward_mjd[j]-mjdref) < STATEMJD_TIMETOL) refpoint = j;
  }
  if(refpoint<0) {
    cerr << "ERROR: reference mjd " << mjdref << " from the Keplerian fit did not match any point in the input Everhart-sampled ephemeris vectors\n";
  }
  cout << "Reference MJD " << mjdref << " corresponds to point " << refpoint << " in the Everhart-sampled ephemeris vectors\n";

  startstate = scale_statevel(starting_statevec,timescale);
  chisq = everchi01(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, starting_statevec, mjdref, obsMJD, observer_statevecs, obsRA, obsDec, sigastrom, fitRA, fitDec, timestep, hnum, hspace, astromrms, verbose);
  if(chisq==LARGERR2) {
    cerr << "ERROR: everchi01 returned error status\n";
    return(5);
  }
  // Initialize eigenvector matrix to the identity matrix
  for(i=0;i<6;i++) {
    for(j=0;j<6;j++) {
      if(i==j) eigenvecs[i][j] = 1.0l;
      else eigenvecs[i][j] = 0.0l;
    }
  }
  cout << "Initial chi-square value is " << chisq << "\n";

  if(outfile.size()>0) outstream1.open(outfile);
  changepar=0;
  for(i=0;i<6;i++) normvec[i] = keepvec[i] = 0.0;
  resetct=covarct=0;
  accepted_states={};
  accepted_all={};
  statechi={};
  for(itct=0;itct<maxnum;itct++) {
    tempstate =  startstate;
    ldval = gaussian_deviate_mt(generator);
    for(i=0;i<6;i++) {
      tempstate[i] += eigenvecs[i][changepar]*stepvec[changepar]*ldval;
    }
    realstate = unscale_statevel(tempstate,timescale);
    newchi = everchi01(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, realstate, mjdref, obsMJD, observer_statevecs, obsRA, obsDec, sigastrom, fitRA, fitDec, timestep, hnum, hspace, astromrms, verbose);
    // Squiggle
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
      if(outfile.size()>0) outstream1 << fixed << setprecision(10) << mjdref << " " << fixed << setprecision(3) << tempstate[0] << " " << tempstate[1] << " " << tempstate[2] << " "  << fixed << setprecision(10) << tempstate[3]/timescale << " " << tempstate[4]/timescale << " " << tempstate[5]/timescale << " " << newchi << "\n";
      if(itct%PRINTNUM==0 || verbose>0) cout << "ACCEPTED" << " itct " << itct << "\n";
      if(newchi<=bestchi) {
	bestchi = newchi;
	beststate = startstate;
	bestRA = fitRA;
	bestDec = fitDec;
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
	if(outfile.size()>0) outstream1 << fixed << setprecision(10) << mjdref << " " << fixed << setprecision(3) << tempstate[0] << " " << tempstate[1] << " " << tempstate[2] << " "  << fixed << setprecision(10) << tempstate[3]/timescale << " " << tempstate[4]/timescale << " " << tempstate[5]/timescale << " " << newchi << "\n";
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
	  covariance_matrix[i][j]/=(double)(accepted_states.size());
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
  if(outfile.size()>0) outstream1.close();

  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromrms << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << mjdref << " : " << fixed << setprecision(3) << beststate[0] << " " << beststate[1] << " " << beststate[2] << " "  << fixed << setprecision(10) << beststate[3]/timescale << " " << beststate[4]/timescale << " " << beststate[5]/timescale << "\n";
  cout << "State vectors correspond to MJD " << mjdref << "\n";

  // Calculate the covariance matrix of the entire accepted sample set
  for(i=0;i<6;i++) {
    for(j=0;j<6;j++) {
      covariance_matrix[i][j] = 0.0l;
      for(k=0;k<long(accepted_all.size());k++) {
	covariance_matrix[i][j] += (accepted_all[k][i]-beststate[i])*(accepted_all[k][j]-beststate[j]);
      }
      covariance_matrix[i][j]/=(double)(accepted_all.size());
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
  double delem = 0.0l;
  long index = 0;
  double_index ldi = double_index(delem,index);
  vector <double_index> ldivec0 ={};
  
  // Take the dot-product of eigenvector 0 with every accepted state,
  // and store in the vector ldivec0.
  for(i=0;i<long(accepted_all.size());i++) {
    // Take the dot product of (accepted_all[i]-beststate) with eigenvector 0
    delem = 0.0l;
    for(j=0;j<6;j++) delem += (accepted_all[i][j]-beststate[j])*eigenvecs[j][0];
    index = i;
    ldi = double_index(delem,index);
    ldivec0.push_back(ldi);
  }

  // Sort ldivec0 by increasing value of the dot-product
  sort(ldivec0.begin(),ldivec0.end(),lower_double_index());
  // Calculate the full range in the dot-product with this eigenvector
  eigenrange0 = ldivec0[ldivec0.size()-1].delem - ldivec0[0].delem;
  main_sample_step = eigenrange0/double(repnum);
  long eigenrepnum0 = repnum;
  long eigenrepct0 = 0;
  if(DEBUG>0) cout << "eigenrange0 = " << eigenrange0 << ", main_sample_step = " << main_sample_step << "\n";
  // Split ldivec0 into eigenrepnum0 equally spaced parts
  long i0=0;
  for(eigenrepct0=1;eigenrepct0<=eigenrepnum0;eigenrepct0++) {
    vector <double_index> ldivec1 ={};
    while(i0<long(ldivec0.size()) && ldivec0[i0].delem < ldivec0[0].delem + eigenrange0*double(eigenrepct0)/double(eigenrepnum0)) {
      // Take the dot-product of this state relative to eigenvector 1
      k = ldivec0[i0].index;
      delem = 0.0l;
      for(j=0;j<6;j++) delem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][1];
      index = k;
      ldi = double_index(delem,index);
      ldivec1.push_back(ldi);
      i0++;
    }
    // We are working with a subset of points, split by the direction defined by eigenvector 0
    if(DEBUG>0) cout << "Loaded ldivec1 with " << ldivec1.size() << " points\n";
    eigenrange1 = 0.0l;
    // Sort ldivec1 by increasing value of the dot-product with eigenvector 1
    if(ldivec1.size()>0) {
      sort(ldivec1.begin(),ldivec1.end(),lower_double_index());
      eigenrange1 = ldivec1[ldivec1.size()-1].delem - ldivec1[0].delem;
      if(DEBUG>0) cout << "eigenrepct0 = " << eigenrepct0 << ", eigenrange1 = " << eigenrange1 << ", main_sample_step = " << main_sample_step << "\n";
    }
    if(ldivec1.size()>0 && eigenrange1>main_sample_step) {
      // Must divide this subset further
      long eigenrepnum1 = eigenrange1/main_sample_step + 1.0;
      long eigenrepct1 = 0;
      long i1=0;
      if(DEBUG>0) cout << "eigenrepnum1 = " << eigenrepnum1 << "\n";
      for(eigenrepct1=1;eigenrepct1<=eigenrepnum1;eigenrepct1++) {
	vector <double_index> ldivec2 ={};
	while(i1<long(ldivec1.size()) && ldivec1[i1].delem < ldivec1[0].delem + eigenrange1*double(eigenrepct1)/double(eigenrepnum1)) {
	  // Take the dot-product of this state relative to eigenvector 2
	  k = ldivec1[i1].index;
	  delem = 0.0l;
	  for(j=0;j<6;j++) delem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][2];
	  index = k;
	  ldi = double_index(delem,index);
	  ldivec2.push_back(ldi);
	  i1++;
	}
	// We are working with a subset of points, split by the direction defined by eigenvector 1
	// Sort ldivec2 by increasing value of the dot-product with eigenvector 2
	eigenrange2 = 0.0l;
	if(ldivec2.size()>0) {
	  sort(ldivec2.begin(),ldivec2.end(),lower_double_index());
	  eigenrange2 = ldivec2[ldivec2.size()-1].delem - ldivec2[0].delem;
	  if(DEBUG>0) cout << "eigenrepct1 = " << eigenrepct1 << ", eigenrange2 = " << eigenrange2 << ", main_sample_step = " << main_sample_step << "\n";
	}
	if(ldivec2.size()>0 && eigenrange2>main_sample_step) {
	  // Must divide this subset further
	  long eigenrepnum2 = eigenrange2/main_sample_step + 1.0;
	  long eigenrepct2 = 0;
	  long i2=0;
	  if(DEBUG>0) cout << "eigenrepnum2 = " << eigenrepnum2 << "\n";
	  for(eigenrepct2=1;eigenrepct2<=eigenrepnum2;eigenrepct2++) {
	    vector <double_index> ldivec3 ={};
	    while(i2<long(ldivec2.size()) && ldivec2[i2].delem < ldivec2[0].delem + eigenrange2*double(eigenrepct2)/double(eigenrepnum2)) {
	      // Take the dot-product of this state relative to eigenvector 3
	      k = ldivec2[i2].index;
	      delem = 0.0l;
	      for(j=0;j<6;j++) delem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][3];
	      index = k;
	      ldi = double_index(delem,index);
	      ldivec3.push_back(ldi);
	      i2++;
	    }
	    // We are working with a subset of points, split by the direction defined by eigenvector 2
	    // Sort ldivec3 by increasing value of the dot-product with eigenvector 3
	    eigenrange3 = 0.0l;
	    if(ldivec3.size()>0) {
	      sort(ldivec3.begin(),ldivec3.end(),lower_double_index());
	      eigenrange3 = ldivec3[ldivec3.size()-1].delem - ldivec3[0].delem;
	      if(DEBUG>0) cout << "eigenrepct2 = " << eigenrepct2 << ", eigenrange3 = " << eigenrange3 << ", main_sample_step = " << main_sample_step << "\n";
	    }
	    if(ldivec3.size()>0 && eigenrange3>main_sample_step) {
	      // Must divide this subset further
	      long eigenrepnum3 = eigenrange3/main_sample_step + 1.0;
	      long eigenrepct3 = 0;
	      long i3=0;
	      if(DEBUG>0) cout << "eigenrepnum3 = " << eigenrepnum3 << "\n";
	      for(eigenrepct3=1;eigenrepct3<=eigenrepnum3;eigenrepct3++) {
		vector <double_index> ldivec4 ={};
		while(i3<long(ldivec3.size()) && ldivec3[i3].delem < ldivec3[0].delem + eigenrange3*double(eigenrepct3)/double(eigenrepnum3)) {
		  // Take the dot-product of this state relative to eigenvector 4
		  k = ldivec2[i3].index;
		  delem = 0.0l;
		  for(j=0;j<6;j++) delem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][4];
		  index = k;
		  ldi = double_index(delem,index);
		  ldivec4.push_back(ldi);
		  i3++;
		}
		// We are working with a subset of points, split by the direction defined by eigenvector 3
		eigenrange4 = 0.0l;
		// Sort ldivec4 by increasing value of the dot-product with eigenvector 4
		if(ldivec4.size()>0) {
		  sort(ldivec4.begin(),ldivec4.end(),lower_double_index());
		  eigenrange4 = ldivec4[ldivec4.size()-1].delem - ldivec4[0].delem;
		  if(DEBUG>0) cout << "eigenrepct3 = " << eigenrepct3 << ", eigenrange4 = " << eigenrange4 << ", main_sample_step = " << main_sample_step << "\n";
		}
		if(ldivec4.size()>0 && eigenrange4>main_sample_step) {
		  // Must divide this subset further
		  long eigenrepnum4 = eigenrange4/main_sample_step + 1.0;
		  long eigenrepct4 = 0;
		  long i4=0;
		  if(DEBUG>0) cout << "eigenrepnum4 = " << eigenrepnum4 << "\n";
		  for(eigenrepct4=1;eigenrepct4<=eigenrepnum4;eigenrepct4++) {
		    vector <double_index> ldivec5 ={};
		    while(i4<long(ldivec4.size()) && ldivec4[i4].delem < ldivec4[0].delem + eigenrange4*double(eigenrepct4)/double(eigenrepnum4)) {
		      // Take the dot-product of this state relative to eigenvector 5
		      k = ldivec4[i4].index;
		      delem = 0.0l;
		      for(j=0;j<6;j++) delem += (accepted_all[k][j]-beststate[j])*eigenvecs[j][5];
		      index = k;
		      ldi = double_index(delem,index);
		      ldivec5.push_back(ldi);
		      i4++;
		    }
		    // We are working with a subset of points, split by the direction defined by eigenvector 4
		    eigenrange5 = 0.0l;
		    // Sort ldivec5 by increasing value of the dot-product with eigenvector 5
		    if(ldivec5.size()>0) {
		      sort(ldivec5.begin(),ldivec5.end(),lower_double_index());
		      eigenrange5 = ldivec5[ldivec5.size()-1].delem - ldivec5[0].delem;
		      if(DEBUG>0) cout << "eigenrepct4 = " << eigenrepct4 << ", eigenrange5 = " << eigenrange5 << ", main_sample_step = " << main_sample_step << "\n";
		    }
		    if(ldivec5.size()>0 && eigenrange5>main_sample_step) {
		      // Must divide this subset further
		      long eigenrepnum5 = eigenrange5/main_sample_step + 1.0;
		      long eigenrepct5 = 0;
		      long i5=0;
		      if(DEBUG>0) cout << "eigenrepnum5 = " << eigenrepnum5 << "\n";
		      for(eigenrepct5=1;eigenrepct5<=eigenrepnum5;eigenrepct5++) {
			while(i5<long(ldivec5.size()) && ldivec5[i5].delem < ldivec5[0].delem + eigenrange5*(double(eigenrepct5)-0.5l)/double(eigenrepnum5)) i5++;
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
    cout << fixed << setprecision(10) << "Representative statevec at MJD: " << mjdref << " : " << fixed << setprecision(3) << representative_states[i][0] << " " << representative_states[i][1] << " " << representative_states[i][2] << " "  << fixed << setprecision(10) << representative_states[i][3]/timescale << " " << representative_states[i][4]/timescale << " " << representative_states[i][5]/timescale << "\n";
    outstream1 << fixed << setprecision(10) << mjdref << " " << fixed << setprecision(3) << representative_states[i][0] << " " << representative_states[i][1] << " " << representative_states[i][2] << " "  << fixed << setprecision(10) << representative_states[i][3]/timescale << " " << representative_states[i][4]/timescale << " " << representative_states[i][5]/timescale << "\n";
  }
  outstream1.close();

  return(0);
}


vector <double> scale_statevel(const vector <double> &statevec, double timescale)
{
  vector <double> outstate = statevec;
  for(int k=3;k<6;k++) outstate[k] *= timescale;
  return(outstate);
}

vector <double> unscale_statevel(const vector <double> &statevec, double timescale)
{
  vector <double> outstate = statevec;
  for(int k=3;k<6;k++) outstate[k] /= timescale;
  return(outstate);
}
