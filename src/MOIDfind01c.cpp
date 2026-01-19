// January 08, 2026: MOIDfind01c.cpp
// Given input Keplerian orbital elements in the format of MPCORB.DAT,
// integrate the object's orbit, and calculate its closest approach to Earth during
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
  cerr << "Usage: MOIDfind01c -cfg configfile -mjdstart mjdstart -mjdend mjdend -mpcdat MPCORB_file -outfile outfile -verbose verbosity\n";
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
  string mpcorbfile;
  string outfile,statefile;
  string epochstring;
  double mjdref = 0.0;
  double mjdstart = 0.0;
  double mjdend = 0.0;
  int configread=0;
  int planetnum=0;
  int planetct=0;
  double ldval=0.0;
  long double M, Peri, Node, Incl, e, n, a;
  M = Peri = Node = Incl = e = n = a = 0.0;
  long double mjd_epoch = mjdref;
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
  vector <double> old_statevec;
  vector <double> targMJD;
  vector <vector <double>> targ_statevecs;
  string lnfromfile;
  vector <string> linestringvec;
  ofstream outstream1;
  ofstream outstream2;
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
  double dist,mindist,dminvel,dmintime;
  double vdist,minvdist,vmindist,vmintime;
  vector <double> targnow;
  vector <double> Sunnow;
  vector <double> Earthnow;
  string desig,arc,fullname;
  double H,G;
  long nobs,nopp;
  int reachedeof=0;
  
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
    } else if(string(argv[i]) == "-mpcorb" || string(argv[i]) == "-MPCORB" || string(argv[i]) == "-mpcdat" || string(argv[i]) == "-MPCDAT" || string(argv[i]) == "--mpcorb" || string(argv[i]) == "--MPCORB" || string(argv[i]) == "--MPCORBDAT") {
      if(i+1 < argc) {
	//There is still something to read;
	mpcorbfile=argv[++i];
	i++;
      }
      else {
	cerr << "MPCORB file keyword supplied with no corresponding argument\n";
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

  
  outstream1.open(outfile);
  outstream1 << "aproach_time approach_dist approach_dvel mindvel_time mindvel_dist mindvel desig H G epochstring mjd_epoch M Peri Node Incl e n a nobs nopp arc fullname\n";
  // Read the MPCORB.DAT file
  instream1.open(mpcorbfile);
  if(!instream1) {
    cerr << "can't open input file " << mpcorbfile << "\n";
    return(1);
  }
  reachedeof = 0;
  long linect=0;
  while(reachedeof==0) {
    desig = arc = fullname = epochstring = "";
    H = G = 0.0;
    nobs = nopp = 0;
    a = e = Incl = Node = Peri = M = mjd_epoch = n = 0.0;
    getline(instream1,lnfromfile);
    linect++;
    int badread=0;
    int badorb=0;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad() && lnfromfile.size()>=194) {
      // Read on.
      // Read the designation
      desig = lnfromfile.substr(0,8);
      cout << "desig: " << desig << "\n";
      // Read the absolute magnitude H
      stest = lnfromfile.substr(8,5);
      cout << "H: " << stest << "\n";
      try { H = stod(stest); }
      catch(...) { cerr << "Warning: unable to read absolute magnitude H from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read the slope parameter G
      stest = lnfromfile.substr(15,4);
      cout << "G: " << stest << "\n";
      try { G = stod(stest); }
      catch(...) { cerr << "Warning: unable to read slope parameter G from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read the epoch of the orbit
      epochstring = lnfromfile.substr(20,5);
      cout << "epochstring: " << epochstring << "\n";
      status = packepoch2MJD(epochstring, mjdref);
      mjd_epoch = mjdref;
      if(status!=0) {
	cerr << "Warning: unable to parse epoch string " << epochstring << " on line " << linect << "\n";
	badread=1;
      }
      // Read the mean anomaly at the epoch
      stest = lnfromfile.substr(26,9);
      cout << "M: " << stest << "\n";
      try { M = stold(stest); }
      catch(...) { cerr << "Warning: unable to read mean anomaly from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read the argument of perihelion
      stest = lnfromfile.substr(37,9);
      cout << "Peri: " << stest << "\n";
      try { Peri = stold(stest); }
      catch(...) { cerr << "Warning: unable to read argument of perihelion from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read the longitude of the ascending node
      stest = lnfromfile.substr(48,9);
      cout << "Node: " << stest << "\n";
      try { Node = stold(stest); }
      catch(...) { cerr << "Warning: unable to read longitude of ascending node from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read the orbital inclination
      stest = lnfromfile.substr(59,9);
      cout << "Incl: " << stest << "\n";
      try { Incl = stold(stest); }
      catch(...) { cerr << "Warning: unable to read orbital inclination from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read the orbital eccentricity
      stest = lnfromfile.substr(70,9);
      cout << "e: " << stest << "\n";
      try { e = stold(stest); }
      catch(...) { cerr << "Warning: unable to read orbital eccentricity from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read mean daily motion
      stest = lnfromfile.substr(81,10);
      cout << "n: " << stest << "\n";
      try { n = stold(stest); }
      catch(...) { cerr << "Warning: unable to read mean daily motion from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read orbital semimajor axis
      stest = lnfromfile.substr(94,9);
      cout << "a: " << stest << "\n";
      try { a = stold(stest); }
      catch(...) { cerr << "Warning: unable to read orbital semimajor axis from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read number of observations
      stest = lnfromfile.substr(117,5);
      cout << "nobs: " << stest << "\n";
      try { nobs = stol(stest); }
      catch(...) { cerr << "Warning: unable to read number of observations from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read number of oppositions
      stest = lnfromfile.substr(123,3);
      try { nopp = stol(stest); }
      catch(...) { cerr << "Warning: unable to read number of oppositions from string " << stest << " on line " << linect << "\n";
	badread=1;
      }
      // Read arc
      arc = lnfromfile.substr(127,9);
      // Read full name
      fullname = lnfromfile.substr(167,26);
      if(badread==0) {
	// Proceed with analysis.
	cout  << fixed << setprecision(8) << "MJD at Epoch is " << mjd_epoch << "\n";
	keplerian_orbit keporb=keplerian_orbit(a,e,Incl,Node,Peri,M,mjd_epoch,n);
	point3LD outpos = point3LD(0.0L,0.0L,0.0L);
	point3LD outvel = point3LD(0.0L,0.0L,0.0L);

	// keporb.semimaj_axis        in AU
	// keporb.eccentricity        unitless
	// keporb.inclination         in degrees
	// keporb.long_ascend_node    Longitude of the ascending node, in degrees
	// keporb.arg_perihelion      Argument of perihelion, in degrees
	// keporb.mean_anom           Mean anomaly at the epoch, in degrees
	// keporb.mjd_epoch           Epoch for the orbit in MJD
	// keporb.mean_daily_motion   in degrees/day

	Kepler2dyn(mjd_epoch, keporb, outpos, outvel);
	starting_statevec[0] = outpos.x*AU_KM;
	starting_statevec[1] = outpos.y*AU_KM;
	starting_statevec[2] = outpos.z*AU_KM;
	starting_statevec[3] = outvel.x*AU_KM/SOLARDAY;
	starting_statevec[4] = outvel.y*AU_KM/SOLARDAY;
	starting_statevec[5] = outvel.z*AU_KM/SOLARDAY;

	// Make sure that mjdref matches a point in the Everhart-sampled planet files
	refpoint = -99;
	for(j=0;j<long(planet_forward_mjd.size());j++) {
	  if(fabs(planet_forward_mjd[j]-mjdref) < STATEMJD_TIMETOL) refpoint = j;
	}
	if(refpoint<0) {
	  cerr << "WARNING: MJD at the epoch " << mjdref << " did not match any point in the input Everhart-sampled ephemeris vectors\n";
	  // Select a nearby reference time that DOES match, if possible.
	  j = long(double(planet_forward_mjd.size())/4.0)*4;
	  while(j>=long(planet_forward_mjd.size())) j-=4;
	  while(j>=0 && planet_forward_mjd[j]>mjdref) j-=4;
	  if(j<0 || planet_forward_mjd[j]>mjdref) {
	    cerr << "WARNING: failed to find a point near mjdref = " << mjdref << "\n";
	    badorb=1;
	  } else refpoint=j;
	  if(badorb==0) {
	    // Integrate the Keplerian state vectors from their current reference MJD
	    // to the MJD corresponding to planetfile_refpoint
	    cout << "A short Keplerian integration will convert the input orbit from epoch MJD " << mjdref << "\n";
	    cout << "to MJD " << planet_forward_mjd[refpoint] << ", which matches point " << refpoint << " in the input Everhart-sampled ephemeris vectors\n";
	    old_statevec = starting_statevec;
	    Kepler_univ_int_SV(GMSUN_KM3_SEC2, mjdref, old_statevec, planet_forward_mjd[refpoint], starting_statevec, verbose);
	    mjdref = planet_forward_mjd[refpoint];
	  }
	} else cout << "MJD at the epoch " << mjdref << " corresponds to point " << refpoint << " in the Everhart-sampled ephemeris vectors\n";
	if(badorb==0) {
	  // Convert from heliolcentric to barycentric coordinates by adding in the position of the Sun.
	  planetpos02(mjdref, 3, Sun_mjd, Sun_statevecs, Sunnow);
	  for(i=0;i<6;i++) starting_statevec[i] += Sunnow[i];
	  targMJD = {};
	  targ_statevecs={};
	  status = integrate_statevec04(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, starting_statevec, mjdstart, mjdref, mjdend, targMJD,  targ_statevecs, timestep, hnum, hspace);
	  if(status!=0) {
	    cerr << "WARNING: integrate_statevec04 returned error status " << status << "\n";
	    badorb=1;
	  } else {
	    outstream2.open("outtest.txt");
	    // Write integrated state vectors to output file
	    outstream2 << "MJD X Y Z VX VY VZ\n";
	    for(i=0;i<long(targMJD.size());i++) {
	      outstream2 << fixed << setprecision(10) << targMJD[i] << " ";
	      outstream2 << fixed << setprecision(3) << targ_statevecs[i][0] << " " << targ_statevecs[i][1] << " " << targ_statevecs[i][2] << " ";
	      outstream2 << fixed << setprecision(9) << targ_statevecs[i][3] << " " << targ_statevecs[i][4] << " " << targ_statevecs[i][5] << "\n";
	    }
	    outstream2.close();
  
	    mindist = dminvel = dist = dmintime = LARGERR;
	    vmindist = minvdist = vdist = vmintime = LARGERR;
	    long closept=0;
	    long minvpt=0;
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
	      vdist = sqrt(DSQUARE(targ_statevecs[i][3] - Earth_statevecs[i+earthoff][3]) + DSQUARE(targ_statevecs[i][4] - Earth_statevecs[i+earthoff][4]) + DSQUARE(targ_statevecs[i][5] - Earth_statevecs[i+earthoff][5]));
	      if(dist<mindist) {
		mindist = dist;
		dmintime = targMJD[i];
		dminvel = vdist;
		closept=i;
	      }
	      if(vdist<minvdist) {
		minvdist = vdist;
		vmintime = targMJD[i];
		vmindist = dist;
		minvpt=i;
	      }
	    }
	    cout << "Closest approach is " << mindist << " km at MJD " << dmintime << " at relative velocity " << dminvel << " km/sec\n";
	    cout << "Best velocity match is " << minvdist << " km/sec at MJD " << vmintime << " at a distance of " << vmindist << " km\n";
	    mjdnow = dmintime;
	    Earthnow = {};
	    //cout << "Earthstart: " << mjdnow << " " << Earth_mjd[0] << " " << Earth_mjd[Earth_mjd.size()-1] << "\n";
	    planetpos02(mjdnow, 3, Earth_mjd, Earth_statevecs, Earthnow);
	    // Interpolate target position
	    targnow ={};
	    //cout << "targstart: " << mjdnow << " " << targMJD[0] << " " << targMJD[targMJD.size()-1] << "\n";
	    planetpos02(mjdnow, 3, targMJD, targ_statevecs, targnow);
	    dist = sqrt(DSQUARE(targnow[0] - Earthnow[0]) + DSQUARE(targnow[1] - Earthnow[1]) + DSQUARE(targnow[2] - Earthnow[2]));
	    //cout << "mjd, dist, mindist = " << mjdnow << " " << dist << " " << mindist << "\n\n";
	    
	    if(closept>4) mjdnow = targMJD[closept-2];
	    else mjdnow = targMJD[2];
	    if(closept<long(targMJD.size()-4)) mjdmax = targMJD[closept+2];
	    else mjdmax = targMJD[targMJD.size()-3];
	    while(mjdnow<mjdmax) {
	      // Interpolate Earth position
	      Earthnow = {};
	      //cout << "Earthint: " << mjdnow << " " << Earth_mjd[0] << " " << Earth_mjd[Earth_mjd.size()-1] << "\n";
	      planetpos02(mjdnow, 3, Earth_mjd, Earth_statevecs, Earthnow);
	      // Interpolate target position
	      targnow ={};
	      //cout << "targint: " << mjdnow << " " << targMJD[0] << " " << targMJD[targMJD.size()-1] << "\n";
	      planetpos02(mjdnow, 3, targMJD, targ_statevecs, targnow);
	      dist = sqrt(DSQUARE(targnow[0] - Earthnow[0]) + DSQUARE(targnow[1] - Earthnow[1]) + DSQUARE(targnow[2] - Earthnow[2]));
	      if(dist<mindist) {
		mindist = dist;
		dmintime = mjdnow;
		dminvel = sqrt(DSQUARE(targnow[3] - Earthnow[3]) + DSQUARE(targnow[4] - Earthnow[4]) + DSQUARE(targnow[5] - Earthnow[5]));
	      }
	      mjdnow += FINESTEP;
	    }

	    if(minvpt>4) mjdnow = targMJD[minvpt-2];
	    else mjdnow = targMJD[2];
	    if(minvpt<long(targMJD.size()-4)) mjdmax = targMJD[minvpt+2];
	    else mjdmax = targMJD[targMJD.size()-3];
	    while(mjdnow<mjdmax) {
	      // Interpolate Earth position
	      Earthnow = {};
	      //cout << "Earthint: " << mjdnow << " " << Earth_mjd[0] << " " << Earth_mjd[Earth_mjd.size()-1] << "\n";
	      planetpos02(mjdnow, 3, Earth_mjd, Earth_statevecs, Earthnow);
	      // Interpolate target position
	      targnow ={};
	      //cout << "targint: " << mjdnow << " " << targMJD[0] << " " << targMJD[targMJD.size()-1] << "\n";
	      planetpos02(mjdnow, 3, targMJD, targ_statevecs, targnow);
	      vdist = sqrt(DSQUARE(targnow[3] - Earthnow[3]) + DSQUARE(targnow[4] - Earthnow[4]) + DSQUARE(targnow[5] - Earthnow[5]));
	      if(vdist<minvdist) {
		minvdist = vdist;
		vmintime = mjdnow;
		vmindist = sqrt(DSQUARE(targnow[0] - Earthnow[0]) + DSQUARE(targnow[1] - Earthnow[1]) + DSQUARE(targnow[2] - Earthnow[2]));
	      }
	      mjdnow += FINESTEP;
	    }
	    cout << "Closest approach is " << mindist << " km at MJD " << dmintime << " at relative velocity " << dminvel << " km/sec\n";
	    cout << "Best velocity match is " << minvdist << " km/sec at MJD " << vmintime << " at a distance of " << vmindist << " km\n";
	    outstream1 << fixed << setprecision(6) << dmintime << " ";
	    outstream1 << fixed << setprecision(2) << mindist << " ";
	    outstream1 << fixed << setprecision(6) << dminvel << " ";
	    outstream1 << fixed << setprecision(2) << vmintime << " ";
	    outstream1 << fixed << setprecision(2) << vmindist << " ";
	    outstream1 << fixed << setprecision(6) << minvdist << " ";
	    outstream1 << fixed << setprecision(2) << desig << " " << H << " " << G << " " << epochstring << " ";
	    outstream1 << fixed << setprecision(6) << mjd_epoch << " " << M << " " << Peri << " " << Node << " " << Incl << " ";
	    outstream1 << fixed << setprecision(8) << e << " " << n << " " << a << " " << nobs << " " << nopp << " " << arc << " " << fullname << "\n";
	  }
	}
      }
    } else {
      if(instream1.eof()) reachedeof=1;
      if(instream1.fail()) {
	reachedeof=2;
	//cerr << "ERROR: failure reading input MPCORB file\n";
	return(0);
      }
      if(instream1.bad()) {
	reachedeof=3;
	cerr << "ERROR stream reading input MPCORB file in bad state\n";
	return(3);
      }
    }
  }
  return(0);
}

