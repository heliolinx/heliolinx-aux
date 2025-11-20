// October 11, 2025: testorb05e.cpp:
// Test different versions of orbit integrators for small
// bodies in the Solar System. This version attempts to calculate
// an observational ephemeris. It als uses custom-interpolated
// plantary ephemeris files at the exact times needed by the Everhart
// integrator, in order to speed up the integration by not
// requiring any interpolation of the ephemeris files as
// part of the integration.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define HNUM 4

// Note: configfile contains the masses and ephemerides for all of the planets.
// configfile also contains the name of an observatory code file. The input
// units of the state vector are AU for positions and km/sec for velocities.
// Internally, the program will use AU for positions and AU/timescale for
// velocities. The default value of timescale is one year, but it can be
// set to other values as needed. The observation file must contain
// MJD, RA, Dec, magnitude, astrometric uncertainty, magnitude uncertainty,
// and obscode. The RA and Dec must be in decimal degrees. The astrometric
// uncertainty must be in arcseconds.
static void show_usage()
{
  cerr << "Usage: testorb06e -cfg configfile -observations obsfile -obscode obscodefile -mjdref mjdref -statevec x y z vx vy vz -testtimes testtimes -outfile outfile -verbose verbosity\n";
}

int main(int argc, char *argv[])
{
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
  string outfile;
  double mjdstart = 0.0L;
  double mjdend = 0.0L;
  double mjdref = 0.0L;
  int configread=0;
  int planetnum=0;
  int planetct=0;
  double ldval=0.0L;
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
  vector <vector <double>> targ_statevecs;
  vector <double> starting_statevec;
  vector <double> targMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> obsMJD;
  vector <double> obsTDB;
  string lnfromfile;
  vector <string> linestringvec;
  ofstream outstream1;
  double timestep = 5.0;
  long testtimes = 1;
  vector <vector <double>> tidemat;
  make_dvec(6, starting_statevec);
  vector <vector <double>> vareq_mat;
  int hnum = HNUM;
  vector <double> hspace;
  int verbose = 0;
  vector <double> obsstate;
  vector <double> sunstate;
  vector <vector <double>> observer_statevecs;
  vector <vector <double>> observer_heliostate;
  long obsnum,obsct,refpoint;
  double obslon,plxcos,plxsin;
  double light_travel_time,RA,Dec;
  vector <double> RA_deriv;
  vector <double> Dec_deriv;
  vector <vector <double>> RA_deriv_mat;
  vector <vector <double>> Dec_deriv_mat;
  vector <vector <double>> Aobs_mat;
  vector <vector <double>> Aobs_transpose;
  vector <vector <double>> Qmat;
  vector <vector <double>> Qinv;
  vector <double> resid_B;
  vector <double> AtransposeB;
  vector <double> Xcor;
  double astromrms=0.0;
  vector <double> rmsvec;
  long testct;
  double a,e,incl;
  a=e=incl=0.0;

  if(argc<18) {
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
    } else if(string(argv[i]) == "-mjdref" || string(argv[i]) == "-mr" || string(argv[i]) == "-MJDref" || string(argv[i]) == "--mjdref" || string(argv[i]) == "--MJDref" || string(argv[i]) == "--ModifiedJulianDayref" || string(argv[i]) == "--modifiedjuliandayref") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stold(argv[++i]);
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
	  starting_statevec[j]=stold(argv[++i]);
	  cout << "read sv entry as " << starting_statevec[j] << "\n";
	} else {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
      }
      i++;
    } else if(string(argv[i]) == "-testtimes") {
      if(i+1 < argc) {
	//There is still something to read;
	testtimes=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "testtimes keyword supplied with no corresponding argument\n";
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
  cout << "input reference MJD for state vectors " << mjdref << "\n";
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
  
  // Load the MJD, RA, and Dec vectors needed by arctrace02.
  obsMJD = obsTDB = obsRA = obsDec = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD.push_back(detvec[obsct].MJD);
    obsTDB.push_back(detvec[obsct].MJD + TTDELTAT/SOLARDAY);
    obsRA.push_back(detvec[obsct].RA);
    obsDec.push_back(detvec[obsct].Dec);
  }
  // Check that MJD values are time-ordered.
  for(obsct=1;obsct<obsnum;obsct++) {
    if(obsMJD[obsct]<obsMJD[obsct-1]) {
      cerr << "ERROR: input observations are not properly sorted by time.\n";
      return(2);
    }
    if(verbose>0) cout << "obsct, MJD, RA, Dec: " << obsct << " " << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << "\n";
  }
  // Make sure that MJDref matches a point in the Everhart-sampled planet files
  long imjdref = round(mjdref);
  if(fabs(mjdref - static_cast<double>(imjdref)) > STATEMJD_TIMETOL) {
    cerr << "ERROR: input reference mjd " << mjdref << " does not fall on an integer timestep\n";
    return(3);
  }
  refpoint = -99;
  for(j=0;j<long(planet_forward_mjd.size());j++) {
    if(fabs(planet_forward_mjd[j]-mjdref) < STATEMJD_TIMETOL) refpoint = j;
  }
  if(refpoint<0) {
    cerr << "ERROR: input reference mjd " << mjdref << " did not match any point in the input Everhart-sampled ephemeris vectors\n";
  }
  cout << "Input reference MJD " << mjdref << " corresponds to point " << refpoint << " in the Everhart-sampled ephemeris vectors\n";
  
  if(obsMJD[0] < planet_forward_mjd[0]) {
    cerr << "ERROR: your first observation, at MJD " << obsMJD[0] << " is before the\n";
    cerr << "begining of the planetary ephemeris files, at MJD " << planet_forward_mjd[0] << "\n";
    return(4);
  } if(obsMJD[obsMJD.size()-1] > planet_forward_mjd[planet_forward_mjd.size()-1]) {
    cerr << "ERROR: your last observation, at MJD " << obsMJD[obsMJD.size()-1] << " is after the\n";
    cerr << "end of the planetary ephemeris files, at MJD " << planet_forward_mjd[planet_forward_mjd.size()-1] << "\n";
    return(4);
  }
  // Find an appropriate value for mjdstart, to make sure the integration
  // will begin before any of the observation times.
  i=1;
  while(refpoint-i*hnum >= 0 && planet_forward_mjd[refpoint-i*hnum] >= obsTDB[0]) i++;
  if(refpoint-i*hnum >= 0) mjdstart = planet_forward_mjd[refpoint-i*hnum];
  else mjdstart = planet_forward_mjd[0];
  // Find an appropriate value for mjdend, to make sure the integration
  // will end after all of the observation times.
  i=1;
  while(refpoint+i*hnum < long(planet_forward_mjd.size()) && planet_forward_mjd[refpoint+i*hnum] <= obsTDB[obsTDB.size()-1]) i++;
  if(refpoint+i*hnum < long(planet_forward_mjd.size())) mjdend = planet_forward_mjd[refpoint+i*hnum];
  else mjdend = planet_forward_mjd[planet_forward_mjd.size()-1];

  cout << "Integration will be performed from MJD_TDB " << mjdstart << " to " << mjdend << "\n";
   
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

  // Open output file and write header
  outstream1.open(outfile);
  outstream1.precision(17);  
  outstream1 << "MJD RA Dec\n";

  cout << "Launching obsint_vareq01()\n";
  rmsvec = {};
  for(testct=0;testct<testtimes;testct++) {
    cout << "Speed-test run " << i <<"\n";
    obsint_vareq01(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, starting_statevec, mjdstart, mjdref, mjdend, obsTDB, targ_statevecs, vareq_mat, timestep, hnum, hspace, verbose);

    fitRA = fitDec = {};
    RA_deriv_mat = Dec_deriv_mat = {};
    for(i=0; i<obsnum; i++) {
      // Initial approximation of the coordinates relative to the observer
      vector <double> relpos;
      make_dvec(3,relpos);
      for(k=0;k<3;k++) relpos[k] = targ_statevecs[i][k] - observer_statevecs[i][k];
      // Initial approximation of the observer-target distance
      ldval = nvecabs(relpos);
      // Convert to meters and divide by the speed of light to get the light travel time.
      light_travel_time = ldval*1000.0/CLIGHT;
      // Light-travel-time corrected version of coordinates relative to the observer
      for(k=0;k<3;k++) relpos[k] = targ_statevecs[i][k] - light_travel_time*targ_statevecs[i][3+k] - observer_statevecs[i][k];
      // Project onto the celestial sphere.
      statevec_to_celederiv(relpos, RA, Dec, RA_deriv, Dec_deriv);
      fitRA.push_back(RA);
      fitDec.push_back(Dec);
      RA_deriv_mat.push_back(RA_deriv);
      Dec_deriv_mat.push_back(Dec_deriv);
    }
    // Construct vector of residuals, O-C, and calculate astrometric RMS
    make_dvec(2*obsnum,resid_B);
    astromrms = 0.0;
    for(obsct=0;obsct<obsnum;obsct++) {
      resid_B[2*obsct] = obsRA[obsct] - fitRA[obsct];
      resid_B[2*obsct+1] = obsDec[obsct] - fitDec[obsct];
      astromrms += intpowD(resid_B[2*obsct]*cos(obsDec[obsct]/DEGPRAD),2) + intpowD(resid_B[2*obsct+1],2);
      cout << scientific << "Residvec " << obsct << " " << resid_B[2*obsct] << " " << resid_B[2*obsct+1] << "\n";
    }
    astromrms/=(double)obsnum;
    astromrms = sqrt(astromrms);
    rmsvec.push_back(astromrms);
    outstream1 << scientific << "Iteration " << testct << " astrometric RMS = " << astromrms*3600.0 << " arcsec\n";

    // Calculate the matrix A. First index is rows (2*obsnum), second is columns (6)
    make_dmat(2*obsnum,6,Aobs_mat);
    for(obsct=0;obsct<obsnum;obsct++) {
      for(i=0;i<6;i++) {
	// Derivative of RA[obsct] w.r.t. each of the initial state vectors
	Aobs_mat[2*obsct][i] = 0.0l;
	for(k=0;k<3;k++) Aobs_mat[2*obsct][i] += RA_deriv_mat[obsct][k] * vareq_mat[obsct][k*6+i];
	// Derivative of Dec[obsct] w.r.t. each of the initial state vectors
	Aobs_mat[2*obsct+1][i] = 0.0l;
	for(k=0;k<3;k++) Aobs_mat[2*obsct+1][i] += Dec_deriv_mat[obsct][k] * vareq_mat[obsct][k*6+i];
      }
    }
    Aobs_transpose = {};
    matrix_transpose(Aobs_mat, Aobs_transpose);
    Qmat = {};
    matXmat(Aobs_transpose, Aobs_mat, Qmat);
    Qinv = {};
    status = invertmatrix01(Qmat, 6, Qinv, verbose);
    AtransposeB = {};
    matXvec(Aobs_transpose, resid_B, AtransposeB);
    Xcor = {};
    matXvec(Qinv, AtransposeB, Xcor);
    cout << "Xcor";
    for(k=0;k<3;k++) cout << " " << Xcor[k];
    for(k=3;k<6;k++) cout << " " << Xcor[k];
    cout << "\n";

    // Apply new correction
    for(k=0;k<6;k++) starting_statevec[k] += Xcor[k];
  }

  statevec2kep_easy(GMSUN_KM3_SEC2,starting_statevec, a, e, incl);
  cout << "Final state vector corresponds to a Keplerian orbit with a = " << a/AU_KM << " AU, e = " << e << " and incl = " << incl << " degrees\n";
  cout << fixed << setprecision(6) << "Final RMS is " << astromrms*3600.0 << " arcsec\n";
  
  // Write final best-fit to output file
  outstream1 << "MJD RA Dec RA_resid Dec_resid total_resid\n";
  cout << "MJD RA Dec RA_resid Dec_resid total_resid\n";
  for(obsct=0;obsct<obsnum;obsct++) {
    cout << fixed << setprecision(10) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " ";
    cout << fixed << setprecision(10) << resid_B[2*obsct]*cos(obsDec[obsct]/DEGPRAD)*3600.0 << " " << resid_B[2*obsct+1]*3600.0 << " " << sqrt(intpowD(resid_B[2*obsct]*cos(obsDec[obsct]/DEGPRAD),2) + intpowD(resid_B[2*obsct+1],2))*3600.0 << "\n";
    outstream1 << fixed << setprecision(10) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " ";
    outstream1 << fixed << setprecision(10) << resid_B[2*obsct]*cos(obsDec[obsct]/DEGPRAD)*3600.0 << " " << resid_B[2*obsct+1]*3600.0 << " " << sqrt(intpowD(resid_B[2*obsct]*cos(obsDec[obsct]/DEGPRAD),2) + intpowD(resid_B[2*obsct+1],2))*3600.0 << "\n";
    outstream1 << "Final state vector corresponds to a Keplerian orbit with a = " << a/AU_KM << " AU, e = " << e << " and incl = " << incl << " degrees\n";
    outstream1 << fixed << setprecision(6) << "Final RMS is " << astromrms*3600.0 << " arcsec\n";
  }
  outstream1.close();
  
  return(0);
}
