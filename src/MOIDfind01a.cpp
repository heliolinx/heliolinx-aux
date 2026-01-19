// January 06, 2026: MOIDfind01a.cpp
// Given input state vectors at a moment in time, integrate the object's
// orbit, and calculate its closest approach to Earth during
// a specified period.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define DEBUG 0
#define HNUM 4
#define FINESTEP 0.001

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain observations in heliolinc's hldet format.
static void show_usage()
{
  cerr << "Usage: MOIDfind01a -cfg configfile -mjdstart mjdstart -mjdref mjd_epoch -mjdend mjdend -statevec x y z vx vy vz -outfile outfile -verbose verbosity\n";
}

int main(int argc, char *argv[])
{
  long i,j,k;
  i=j=k=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string Earthfile,Sunfile,planetfile;
  string spaceobsfile;
  string outfile,statefile;
  double mjdref = 0.0;
  double mjdstart = 0.0;
  double mjdend = 0.0;
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
  vector <double> targMJD;
  vector <vector <double>> targ_statevecs;
  string lnfromfile;
  vector <string> linestringvec;
  ofstream outstream1;
  double timestep = 5.0;
  vector <vector <double>> tidemat;
  make_dvec(6, starting_statevec);
  int hnum = HNUM;
  vector <double> hspace;
  int verbose = 0;
  long refpoint;
  double mjdnow=0.0;
  double mjdmax=0.0;
  vector <double> fineMJD;
  vector <vector <double>> fine_statevecs;
  vector <vector <double>> fine_Earth_statevecs;
  double dist,mindist,appvel,mintime;
  vector <double> targnow;
  vector <double> Earthnow;
  
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
    if(string(argv[i]) == "-mjdstart" || string(argv[i]) == "-ms" || string(argv[i]) == "-MJDstart" || string(argv[i]) == "--mjdstart" || string(argv[i]) == "--MJDstart" || string(argv[i]) == "--ModifiedJulianDaystart" || string(argv[i]) == "--modifiedjuliandaystart") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdstart=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD start keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-mjdend" || string(argv[i]) == "-me" || string(argv[i]) == "-MJDend" || string(argv[i]) == "--mjdend" || string(argv[i]) == "--MJDend" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdend=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
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
  


  // Make sure that mjdref matches a point in the Everhart-sampled planet files
  refpoint = -99;
  for(j=0;j<long(planet_forward_mjd.size());j++) {
    if(fabs(planet_forward_mjd[j]-mjdref) < STATEMJD_TIMETOL) refpoint = j;
  }
  if(refpoint<0) {
    cerr << "ERROR: MJD at the epoch " << mjdref << " did not match any point in the input Everhart-sampled ephemeris vectors\n";
  }
  cout << "MJD at the epoch " << mjdref << " corresponds to point " << refpoint << " in the Everhart-sampled ephemeris vectors\n";

  targMJD = {};
  targ_statevecs={};
  status = integrate_statevec04(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, starting_statevec, mjdstart, mjdref, mjdend, targMJD,  targ_statevecs, timestep, hnum, hspace);
  if(status!=0) {
    cerr << "ERROR: integrate_statevec04 returned error status " << status << "\n";
    return(status);
  }
  
  outstream1.open(outfile);
  // Write integrated state vectors to output file
  outstream1 << "MJD X Y Z VX VY VZ\n";
  for(i=0;i<long(targMJD.size());i++) {
    outstream1 << fixed << setprecision(10) << targMJD[i] << " ";
    outstream1 << fixed << setprecision(3) << targ_statevecs[i][0] << " " << targ_statevecs[i][1] << " " << targ_statevecs[i][2] << " ";
    outstream1 << fixed << setprecision(9) << targ_statevecs[i][3] << " " << targ_statevecs[i][4] << " " << targ_statevecs[i][5] << "\n";
  }
  outstream1.close();

  
  mindist = appvel = dist = mintime = LARGERR;
  long closept=0;
  long earthoff=-10;
  for(i=0;i<long(Earth_mjd.size());i++) {
    if(Earth_mjd[i] == targMJD[0]) earthoff=i;
  }
  if(earthoff<0) {
    cerr << "ERROR: Earth MJD vector does not match target vector\n";
    return(2);
  }
  for(i=0;i<long(targMJD.size());i++) {
    dist = sqrt(DSQUARE(targ_statevecs[i][0] - Earth_statevecs[i+earthoff][0]) + DSQUARE(targ_statevecs[i][1] - Earth_statevecs[i+earthoff][1]) + DSQUARE(targ_statevecs[i][2] - Earth_statevecs[i+earthoff][2]));
    if(dist<mindist) {
      mindist = dist;
      mintime = targMJD[i];
      appvel = sqrt(DSQUARE(targ_statevecs[i][3] - Earth_statevecs[i+earthoff][3]) + DSQUARE(targ_statevecs[i][4] - Earth_statevecs[i+earthoff][4]) + DSQUARE(targ_statevecs[i][5] - Earth_statevecs[i+earthoff][5]));
      closept=i;
    }
  }
  cout << "Closest approach is " << mindist << " km at MJD " << mintime << " at relative velocity " << appvel << " km/sec\n";
  mjdnow = mintime;
  Earthnow = {};
  planetpos02(mjdnow, 3, Earth_mjd, Earth_statevecs, Earthnow);
  // Interpolate target position
  targnow ={};
  planetpos02(mjdnow, 3, targMJD, targ_statevecs, targnow);
  dist = sqrt(DSQUARE(targnow[0] - Earthnow[0]) + DSQUARE(targnow[1] - Earthnow[1]) + DSQUARE(targnow[2] - Earthnow[2]));
  cout << "mjd, dist, mindist = " << mjdnow << " " << dist << " " << mindist << "\n\n";

  if(closept>2) mjdnow = targMJD[closept-2];
  else mjdnow = targMJD[0];
  if(closept<long(targMJD.size()-2)) mjdmax = targMJD[closept+2];
  else mjdmax = targMJD[targMJD.size()-1];
  while(mjdnow<mjdmax) {
    // Interpolate Earth position
    Earthnow = {};
    planetpos02(mjdnow, 3, Earth_mjd, Earth_statevecs, Earthnow);
    // Interpolate target position
    targnow ={};
    planetpos02(mjdnow, 3, targMJD, targ_statevecs, targnow);
    dist = sqrt(DSQUARE(targnow[0] - Earthnow[0]) + DSQUARE(targnow[1] - Earthnow[1]) + DSQUARE(targnow[2] - Earthnow[2]));
    //cout << "mjd, dist, mindist = " << mjdnow << " " << dist << " " << mindist << "\n";
    if(dist<mindist) {
      mindist = dist;
      mintime = mjdnow;
      appvel = sqrt(DSQUARE(targnow[3] - Earthnow[3]) + DSQUARE(targnow[4] - Earthnow[4]) + DSQUARE(targnow[5] - Earthnow[5]));
    }
    mjdnow += FINESTEP;
  }
  cout << "Closest approach is " << mindist << " km at MJD " << mintime << " at relative velocity " << appvel << " km/sec\n";

  return(0);
}

