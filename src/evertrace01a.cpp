// October 20, 2025: evertrace01a.cpp
// Given an input file giving observer positions and measured RA and Dec
// for an object at various times, given as MJD, fit the object's orbit.
// Start by identifying a short sequence of observations that can be
// reasonably fit with a Keplerian orbit, and performing this fit using
// the Method of Herget. Then, using the approximate state vectors output
// from the Keplerian fit to perform a full, 6-dimensional orbit fit
// including planetary perturbations. 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define DEBUG 0
#define HNUM 4

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain observations in heliolinc's hldet format.
static void show_usage()
{
  cerr << "Usage: evertrace01a -cfg configfile -observations obsfile -kepspan time_span_for_Keplerian_fit(day) -minchi min_chi_change -rmsthresh astrometric_rms_threshold -obscode obscodefile -maxiter maxiter -outfile outfile -statefile statefile -verbose verbosity\n";
}

int main(int argc, char *argv[])
{
  double kepspan=15.0;
  long double minchichange = 0.001l;
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
  string outfile,statefile;
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
  vector <double> out_statevec;
  vector <double> targMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> fitresid;
  vector <double> obsMJD;
  vector <double> obsTDB;
  vector <double> sigastrom;
  string lnfromfile;
  vector <string> linestringvec;
  ofstream outstream1;
  double timestep = 5.0;
  long maxiter = 10;
  vector <vector <double>> tidemat;
  make_dvec(6, starting_statevec);
  int hnum = HNUM;
  vector <double> hspace;
  int verbose = 0;
  vector <double> obsstate;
  vector <double> sunstate;
  vector <vector <double>> observer_statevecs;
  vector <vector <double>> observer_heliostate;
  long obsnum,obsct,refpoint;
  double obslon,plxcos,plxsin,dist;
  double astromrms=0.0;
  double astromrmsthresh = 0.1;

  vector <vector <double>> Kepobserverpos;
  vector <double> KepMJD;
  vector <double> KepRA;
  vector <double> KepDec;
  vector <double> Kepsig;
  double geodist1 = 1.0;
  double geodist2 = 1.1;
  double ftol = FTOL_HERGET_SIMPLEX;
  double simplex_scale = SIMPLEX_SCALEFAC;
  int bestpoint,point1;
  int simptype,point2,kepnum,kepmax;
  bestpoint=point1=0;
  simptype=point2=kepnum=kepmax=1;
  double kepmetric,kepmetbest;
  kepmetric=kepmetbest=1.0;
  double stateMJD = 0.0;
  long itnum = 0;
  double chisq=0.0l;

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
    if(string(argv[i]) == "-o" || string(argv[i]) == "-obs" || string(argv[i]) == "-obsfile" || string(argv[i]) == "-observations" || string(argv[i]) == "--obsfile" || string(argv[i]) == "--observation" || string(argv[i]) == "--observationfile" || string(argv[i]) == "--observations") {
      if(i+1 < argc) {
	// There is still something to read;
	obsfile=argv[++i];
	i++;
      } else {
	cerr << "Observation file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-kepspan") {
      if(i+1 < argc) {
	//There is still something to read;
	kepspan=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Keplerian time-span keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minchi") {
      if(i+1 < argc) {
	//There is still something to read;
	minchichange=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum chi change keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rmsthresh") {
      if(i+1 < argc) {
	//There is still something to read;
	astromrmsthresh=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Astrometric RMS threshold keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-maxiter") {
      if(i+1 < argc) {
	//There is still something to read;
	maxiter=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "max iterations keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-statefile" || string(argv[i]) == "-state") {
      if(i+1 < argc) {
	//There is still something to read;
	statefile=argv[++i];
	i++;
      }
      else {
	cerr << "State vector file keyword supplied with no corresponding argument\n";
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
  cout << "input configuration file " << configfile << "\n";
  cout << "input observations file " << obsfile << "\n";
  cout << "input observatory code file " << obscodefile << "\n";
  cout << "output file " << outfile << "\n";

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
    sigastrom.push_back(1.0);
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
  
  // SET UP FOR PRELIMINARY KEPLERIAN FIT
  // Find the time period of length kepspan
  // with the largest number of observations
  bestpoint=kepnum=kepmax=0;
  kepmetric=kepmetbest=0.0;
  for(i=0;i<obsnum;i++) {
    j=i;
    kepnum=0;
    kepmetric=0.0;
    while(j<obsnum && obsMJD[j]-obsMJD[i] < kepspan) {
      kepnum++;
      j++;
    }
    kepmetric = double(kepnum)*(obsMJD[j-1]-obsMJD[i]); 
    if(kepmetric>kepmetbest) {
      kepmetbest=kepmetric;
      kepmax=kepnum;
      bestpoint=i;
    }
  }
  kepnum=kepmax;
  cout << "Best sequence contains " << kepnum << " points, and runs from obsMJD[" << bestpoint << "]=" << obsMJD[bestpoint] << " to obsMJD[" << bestpoint+kepnum-1 << "]=" << obsMJD[bestpoint+kepnum-1] << "\n";
  KepMJD = KepRA = KepDec = Kepsig = {};
  Kepobserverpos = {};
  for(i=bestpoint;i<bestpoint+kepnum;i++) {
    KepMJD.push_back(obsMJD[i]);
    Kepobserverpos.push_back(observer_statevecs[i]);
    KepRA.push_back(obsRA[i]);
    KepDec.push_back(obsDec[i]);
    Kepsig.push_back(sigastrom[i]);
  }
  cout << "kepnum = " << kepnum << " " << KepMJD.size() << "\n";
  for(i=0;i<kepnum;i++) {
    cout << Kepobserverpos[i][1] << " " << Kepobserverpos[i][2] << " " << Kepobserverpos[i][3] << " " << KepMJD[i] << " " << KepRA[i] << " " << KepDec[i] << "\n";
  }
  fitDec = fitRA = fitresid = out_statevec = {};
  cout << "Launching Hergetfit_vstarSV\n";
  chisq = Hergetfit_vstarSV(geodist1, geodist2, simplex_scale, simptype, ftol, 1, kepnum, Kepobserverpos, KepMJD, KepRA, KepDec, Kepsig, fitRA, fitDec, fitresid, out_statevec, stateMJD, itnum, verbose);
  cout << "Keplerian fit produced chisq = " << chisq << "\n";

  // Choose a suitable reference MJD near the middle of the Keplerian fit.
  // Require that it also fall exactly on an integer time step for the Everhart
  // sampling read from the configuration file.
  mjdref = 0.5*(KepMJD[0]+KepMJD[kepnum-1]);
  mjdref = round(mjdref/timestep)*timestep;
  // Make sure that mjdref matches a point in the Everhart-sampled planet files
  refpoint = -99;
  for(j=0;j<long(planet_forward_mjd.size());j++) {
    if(fabs(planet_forward_mjd[j]-mjdref) < STATEMJD_TIMETOL) refpoint = j;
  }
  if(refpoint<0) {
    cerr << "ERROR: reference mjd " << mjdref << " from the Keplerian fit did not match any point in the input Everhart-sampled ephemeris vectors\n";
  }
  cout << "Reference MJD " << mjdref << " corresponds to point " << refpoint << " in the Everhart-sampled ephemeris vectors\n";

  // Integrate the Keplerian state vectors from their current reference MJD
  // to the MJD corresponding to planetfile_refpoint
  Kepler_univ_int_SV(GMSUN_KM3_SEC2, stateMJD, out_statevec, mjdref, starting_statevec, verbose);
  if(verbose>0) cout << "Kepler_univ_int_SV output statevec for MJD " << mjdref << ":\n";
  cout << starting_statevec[0] << " " << starting_statevec[1] << " " << starting_statevec[2] << " " << starting_statevec[3] << " " << starting_statevec[4] << " " << starting_statevec[5] << "\n";

  out_statevec={};
  status = evertrace01(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, starting_statevec, mjdref, obsMJD, observer_statevecs, obsRA, obsDec, sigastrom, fitRA, fitDec, out_statevec, timestep, hnum, hspace, minchichange, astromrmsthresh, maxiter, itnum, astromrms, chisq, verbose);
  if(status!=0) {
    cerr << "ERROR: evertrace01 returned error status " << status << "\n";
    return(status);
  }
  
  outstream1.open(outfile);
  // Write final best-fit to output file
  outstream1 << "MJD RA Dec RA_resid Dec_resid total_resid\n";
  cout << "MJD RA Dec RA_resid Dec_resid total_resid\n";
  astromrms = chisq = 0.0;
  for(obsct=0;obsct<obsnum;obsct++) {
    dist = 3600.0*distradec01(obsRA[obsct],obsDec[obsct],fitRA[obsct],fitDec[obsct]);
    astromrms += dist*dist;
    chisq += dist*dist/sigastrom[obsct]/sigastrom[obsct];
    if(verbose>0) cout << fixed << setprecision(10) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << fitRA[obsct] << " " << fitDec[obsct] << " ";
    if(verbose>0) cout << fixed << setprecision(10) << (obsRA[obsct]-fitRA[obsct])*cos(obsDec[obsct]/DEGPRAD)*3600.0 << " " << (obsDec[obsct]-fitDec[obsct])*3600.0 << " " << dist << " " << detvec[obsct].obscode << " " << detvec[obsct].mag << " " << detvec[obsct].band << "\n";
    outstream1 << fixed << setprecision(10) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << fitRA[obsct] << " " << fitDec[obsct] << " ";
    outstream1 << fixed << setprecision(10) << (obsRA[obsct]-fitRA[obsct])*cos(obsDec[obsct]/DEGPRAD)*3600.0 << " " << (obsDec[obsct]-fitDec[obsct])*3600.0 << " " << dist << " " << detvec[obsct].obscode << " " << detvec[obsct].mag << " " << detvec[obsct].band << "\n";
  }
  outstream1.close();
  astromrms = sqrt(astromrms/double(obsnum));  
  cout << "Writing state vector filed called " << statefile << "\n";
  outstream1.open(statefile);
  outstream1 << fixed << setprecision(10) << planet_forward_mjd[refpoint] << " " << fixed << setprecision(3) << out_statevec[0] << " " << out_statevec[1] << " " << out_statevec[2] << " "  << fixed << setprecision(10) << out_statevec[3] << " " << out_statevec[4] << " " << out_statevec[5] << " " << fixed << setprecision(5) << astromrms << "\n";
  outstream1.close();

  double a,e,incl;
  statevec2kep_easy(GMSUN_KM3_SEC2,starting_statevec, a, e, incl);
  cout << "Final state vector corresponds to a Keplerian orbit with a = " << a/AU_KM << " AU, e = " << e << " and incl = " << incl << " degrees\n";
  cout << fixed << setprecision(6) << "Final RMS is " << astromrms << " arcsec, chi square = " << chisq << "\n";
  return(0);
}

