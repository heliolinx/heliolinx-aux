// April 20, 2026
// Given a representative MCMC file from MCMC_everhart01a.cpp,
// integrate all of the orbits to a specific time and output an
// observational uncertainty cloud.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define MINCOLS 9
#define DEBUG 0
#define TIMERANGEFAC 2000000000
#define TWINFIT_DEFAULT 10
#define HNUM 4
#define MINTRACKSEP 1.0 // Minimum separation for a valid pair of tracklets, in days.
                        // If they are closer than this, they are not sufficiently independent,
                        // and do not validate the reality of the linkage.

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain
// MJD, observer X, Y, Z, VX, VY, VZ, RA, and Dec. The observer coordinates and velocities must
// be in km and km/sec, relative to the Sun. The RA and Dec must be in decimal degrees.
static void show_usage()
{
  cerr << "Usage: everephem01a -cfg configfile -infile infile -obscode obscodefile -ephemcode ephemcode -mjd mjd -outfile outfile -verbose verbosity\n";
}


int main(int argc, char *argv[])
{
  ifstream instream1;
  string infile;
  string configfile;
  string outfile;
  long i, j, k;
  int status, planetnum, planetct, configread;
  long badread,reachedeof,obsnum,obsct,tct;
  i=j=status=configread=badread=reachedeof=obsnum=obsct=tct=0;
  vector <observatory> observatory_list = {};
  string obscodefile;
  string ephemeris_code;
  char ephemcode[MINSTRINGLEN];
  int hnum = HNUM;
  int verbose = 0;
  int verbose2 = 0;
  vector <double> hspace;
  string stest;
  double ldval=0.0L;
  double timestep = 5.0;
  vector <double> planetmasses;
  vector <double> mjdtest;
  vector <vector <double>> statevecs;
  vector <vector <double>> Earth_statevecs;
  vector <vector <double>> Sun_statevecs;
  vector <double> Earth_mjd;
  vector <double> Sun_mjd;
  vector <vector <double>> targ_statevecs;
  double mjd_statevecs;
  vector <double> planet_backward_mjd;
  vector <vector <double>> planet_backward_statevecs;
  vector <double> planet_forward_mjd;
  vector <vector <double>> planet_forward_statevecs;
  vector <double> obsstate;
  vector <double> sunstate;
  vector <vector <double>> eph_obsstate;
  vector <double> one_statevec;
  vector <vector <double>> starting_statevecs;
  vector <double> out_statevec;
  vector <double> master_statevec;
  vector <double> best_statevec;
  vector <double> obsMJD;
  vector <double> fitRA;
  vector <double> fitDec;
  string lnfromfile;
  vector <string> linestringvec;
  vector <double> mjd_sv_vec;
  double mjdstart = 0;
  double mjdend = 0;
  double light_travel_time = 0.0l;
  double timediff = 0.0l;
  string Earthfile,Sunfile,planetfile;
  double dist,pa;
  long startpoint,endpoint;
  long repnum,repct;

  double ephmjd;
  double mjdref;
  double RAvel,Decvel;
  ofstream outstream1;
  double obslon = 0.0l;
  double plxcos = 0.0l;
  double plxsin = 0.0l;
  
 
  make_dvec(6,one_statevec);

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
	instream1.open(configfile);
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
	Sun_mjd={};
	Sun_statevecs={};
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
	Earth_mjd={};
	Earth_statevecs={};
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
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-in" || string(argv[i]) == "-inputfile" || string(argv[i]) == "-input_file" || string(argv[i]) == "--infile" || string(argv[i]) == "--input_file" || string(argv[i]) == "-inf" || string(argv[i]) == "--inputfile") {
      if(i+1 < argc) {
	// There is still something to read;
	infile=argv[++i];
	i++;
      } else {
	cerr << "Input file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-mjd") {
      if(i+1 < argc) {
	//There is still something to read;
	ephmjd=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Ephemeris MJD keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-ephemcode") {
      if(i+1 < argc) {
	//There is still something to read;
	ephemeris_code = argv[++i];
	// Copy input value for obscode, making sure it's not too long
	strncpy(ephemcode, ephemeris_code.c_str(), sizeof(ephemcode));
	ephemcode[sizeof(ephemcode)-1] = 0;
	i++;
      }
      else {
	cerr << "Ephemeris obscode keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--outfile") {
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
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose2=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
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
  cout << "input MCMC file " << infile << "\n";
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

  if(verbose2<=0) verbose=0;
  else verbose = verbose2-1;
  
  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, verbose);
  if(status!=0) {
    cerr << "ERROR reading observatory code file " << obscodefile << "\n";
    return(1);
  }
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";

  // Read the ibout representative state vector file
  starting_statevecs = {};
  instream1.open(infile);
  if(!instream1) {
    cerr << "can't open input file " << infile << "\n";
    return(1);
  }
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Read on.
      // Read the starting MJD
      startpoint=0;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { mjd_statevecs = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read MJD " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the six components of the state vector
      for(i=0;i<6;i++) {
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { one_statevec[i] = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read state vector x component " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	} else badread=1;
      }
      if(badread==0) {
	mjd_sv_vec.push_back(mjd_statevecs);
	starting_statevecs.push_back(one_statevec);
      }
    } else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(badread!=0) {
      cerr << "ERROR reading representative state vector file " << infile << "\n";
      cerr << "Last point was " << starting_statevecs.size() << "; last file line was " << lnfromfile << "\n";
      return(badread);
    }
  }
  instream1.close();

  if(badread!=0) {
    cerr << "ERROR reading representative state vector file " << infile << "\n";
    return(badread);
  }
  if(reachedeof==1) {
    if(verbose>=1) cout << "Input file " << infile << " read successfully to the end.\n";
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << infile << " before the end\n";
  } else if(reachedeof==-1) {
    cerr << "ERROR: file read failed\n";
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
  }
  repnum = starting_statevecs.size();
  if(repnum!=long(mjd_sv_vec.size())) {
    cerr << "ERROR: vector mismatch from representative state vector file\n";
    cerr << repnum << " != " << mjd_sv_vec.size() << "\n";
    return(2);
  }
  cout << "Read " << repnum << " lines from representative state vector file " << infile << "\n";
  // Calculate observer's barycentric position at ephmjd
  status = obscode_lookup(observatory_list,ephemcode,obslon,plxcos,plxsin);
  observer_barystate01(ephmjd, 5, obslon, plxcos, plxsin, Earth_mjd, Earth_statevecs, obsstate, verbose);
  cout << "Observer state vector at time " << ephmjd << " is " << obsstate[0] << " " << obsstate[1] << " " << obsstate[2] << " " << obsstate[3] << " " << obsstate[4] << " " << obsstate[5] << "\n";
  // Calculate RA and Dec for all of the orbit clones at time ephmjd
  outstream1.open(outfile);
  for(repct=0;repct<repnum;repct++) {
    one_statevec = starting_statevecs[repct];
    mjdref = mjd_sv_vec[repct];
    long refpoint = -99;
    for(j=0;j<long(planet_forward_mjd.size());j++) {
      if(fabs(planet_forward_mjd[j]-mjdref) < STATEMJD_TIMETOL) refpoint = j;
    }
    if(refpoint<0) {
      cerr << "ERROR: Input reference mjd " << mjdref << " did not match any point in the input Everhart-sampled ephemeris vectors\n";
      return(1);
    }
    if(verbose>0) cout << "Input reference MJD " << mjdref << " corresponds to point " << refpoint << " in the Everhart-sampled ephemeris vectors\n";
    i = refpoint;
    mjdstart = planet_forward_mjd[i];
    while(i>=hnum && (mjdstart>=ephmjd-timestep || mjdstart>=mjdref-timestep)) {
      i-=hnum;
      mjdstart = planet_forward_mjd[i];
    }
    i=refpoint;
    mjdend = planet_forward_mjd[i];
    while(i<long(planet_forward_mjd.size()-hnum) && (mjdend<=ephmjd+timestep || mjdend<=mjdref+timestep)) {
      i+=hnum;
      mjdend = planet_forward_mjd[i];
    }
    obsMJD = {};
    obsMJD.push_back(ephmjd);
    obsMJD.push_back(ephmjd+0.1);
   targ_statevecs = {};
    cout << "Integrating representative state vector " << one_statevec[0] << " "<< one_statevec[1] << " "<< one_statevec[2] << " "<< one_statevec[3] << " "<< one_statevec[4] << " "<< one_statevec[5] << "\n";
    cout << "from MJD " << mjdref << " to " << obsMJD[0] << "\n";
    obsint_everuse01(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, one_statevec, mjdstart, mjdref, mjdend, obsMJD,  targ_statevecs, timestep, hnum, hspace, verbose);
    cout << "Output state vectors " << targ_statevecs[0][0] << " " << targ_statevecs[0][1] << " " << targ_statevecs[0][2] << " " << targ_statevecs[0][3] << " " << targ_statevecs[0][4] << " " << targ_statevecs[0][5] << "\n";
    // Initial approximation of the coordinates relative to the observer
    vector <double> relpos;
    double RA1,RA2,Dec1,Dec2;
    make_dvec(3,relpos);
    for(k=0;k<3;k++) relpos[k] = targ_statevecs[0][k] - obsstate[k];
    // Initial approximation of the observer-target distance
    ldval = nvecabs(relpos);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    for(k=0;k<3;k++) relpos[k] = targ_statevecs[0][k] - light_travel_time*targ_statevecs[0][3+k] - obsstate[k];
    // Project onto the celestial sphere.
    statevec_to_celestial(relpos, RA1, Dec1);
    // Calculate the position a little bit later
    for(k=0;k<3;k++) relpos[k] += TTDELTAT*(targ_statevecs[0][3+k] - obsstate[3+k]);
    // Project onto the celestial sphere.
    statevec_to_celestial(relpos, RA2, Dec2);
    timediff = TTDELTAT/SOLARDAY; // TTDELTAT is in seconds, timediff is in days
    distradec02(RA1, Dec1, RA2, Dec2, &dist, &pa);
    RAvel = dist*sin(pa/DEGPRAD)/timediff; // Degrees per day
    Decvel = dist*cos(pa/DEGPRAD)/timediff; // Degrees per day
    outstream1 << fixed << setprecision(6) << RA1 << " " << Dec1 << " " << RAvel << " " << Decvel << "\n";
  }
  outstream1.close();
  return(0);
}

