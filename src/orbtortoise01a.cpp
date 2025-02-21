// April 06, 2022: orbtortoise01a.cpp
// Fit an input observations file using the downhill simplex method.
// The observation file must contain MJD, RA, Dec,
// astrometric uncertainty, obscode.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3
#define MJDCOL 3
#define RACOL 6
#define DECCOL 8
#define MINOBSINTERVAL 1.0 // Minimum time-between-images in seconds
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME 1.5 // Default max inter-image time interval
                    // for tracklets, in hours (will be converted
                    // to days before use).
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define DEBUG 0
#define TOLSCALE 50.0 // if 

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
  cerr << "Usage: orbtortoise01a -cfg configfile -observations obsfile -mjd mjdstart -statevec x y z vx vy vz -statestep starting step for state vectors -timescale timescale to convert velocity to distance -ftol ftol -bumpnum 1000 -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string obsfile;
  string observatory_code_file;
  string planetfile;
  string outfile;
  long double x,y,z,vx,vy,vz;
  x=y=z=vx=vy=vz=0.0L;
  long double simplex_side = 0.0L;
  long double timescale = 0.0L;
  long double mjdstart = 0.0L;
  int configread=0;
  int polyorder=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
  long double ldval=0.0L;
  double dval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> planetpos_reverse;
  vector <point3LD> Earthpos;
  vector <point3LD> Sunpos;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  vector <point3LD> observerpos;
  vector <long double> obsMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> magnitude;
  vector <double> sigastrom;
  vector <double> sigmag;
  vector <string> obscodevec;
  int obsnum=0;
  int obsct=0;
  string lnfromfile;
  vector <string> linestringvec;
  int MJDcol,RAcol,Deccol,magcol,sigastromcol,sigmagcol,obscodecol;
  MJDcol=RAcol=Deccol=magcol=sigastromcol=sigmagcol=obscodecol=0;
  int badread=0;
  int reachedeof=0;
  vector <observatory> observatory_list = {};
  char obscode[MINSTRINGLEN];
  double obslon = 289.26345L;
  double plxcos = 0.865020L;
  double plxsin = -0.500901L;
  point3LD outpos = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  long double outRA,outDec;
  int fieldnum=0;
  vector <long double> scalestate;
  make_LDvec(60,scalestate);
  vector <long double> newstate;
  make_LDvec(60,newstate);
  vector <long double> chivec;
  make_LDvec(70,chivec);
  vector <long double> beststatevec;
  long double chisq=0L;
  long double newchi=0L;
  long double bestchi=0L;
  long double worstchi=0L;
  int worstpoint=0;
  int bestpoint=0;
  vector <vector <long double>> simplex;
  make_LDmat(70,60,simplex);
  long double simprange=0L;
  long double ftol=1e-5;
  int countworse=0;
  int bumpnum=0;
  int bumpct=0;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> fitresid;
  vector <double> bestRA;
  vector <double> bestDec;
  vector <double> bestresid;
  long double globalchimin=LARGERR;

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
	// Read the name of the observatory code file
	status=readconfigstring(instream1,observatory_code_file);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigstring(instream1,observatory_code_file);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observatory code file read as " << observatory_code_file << "\n";      
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
	  if(planetct == pctEarth) Earthpos = temppos;
	  if(planetct == pctSun) Sunpos = temppos;
	  cout << "Finished reading ephemeris file " << planetfile << "\n";
	}
	// Read default state vector x coordinate from config file
	status=readconfigLD(instream1,&x);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&x);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default state vector coordinate x read as " << x << " AU\n";
	// Read default state vector y coordinate from config file
	status=readconfigLD(instream1,&y);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&y);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default state vector coordinate y read as " << y << "AU\n";
	// Read default state vector z coordinate from config file
	status=readconfigLD(instream1,&z);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&z);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default state vector coordinate z read as " << z << "AU\n";
	// Read default state vector vx coordinate from config file
	status=readconfigLD(instream1,&vx);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&vx);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default state vector coordinate vx read as " << vx << " km/sec\n";
	// Read default state vector vy coordinate from config file
	status=readconfigLD(instream1,&vy);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&vy);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default state vector coordinate vy read as " << vy << " km/sec\n";
	// Read default state vector vz coordinate from config file
	status=readconfigLD(instream1,&vz);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&vz);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default state vector coordinate vz read as " << vz << " km/sec\n";
	// Read initial side length for the simplex
	status=readconfigLD(instream1,&simplex_side);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&simplex_side);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Initial side length for simplex read as " << simplex_side << "\n";
	// Read timescale for converting velocities to distance units
	status=readconfigLD(instream1,&timescale);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&timescale);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Timescale for converting velocities to\ndistance units read as " << timescale << " days\n";
	// Read observation file column holding the MJD;
	status=readconfigint(instream1,&MJDcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&MJDcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding MJD is " << MJDcol << " \n";
	// Read observation file column holding the RA;
	status=readconfigint(instream1,&RAcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&RAcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding RA is " << RAcol << " \n";
	// Read observation file column holding the Dec;
	status=readconfigint(instream1,&Deccol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&Deccol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding Dec is " << Deccol << " \n";
	// Read observation file column holding the magnitude;
	status=readconfigint(instream1,&magcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&magcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding magnitude is " << magcol << " \n";
	// Read observation file column holding the astrometric uncertainty;
	status=readconfigint(instream1,&sigastromcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&sigastromcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding astrometric uncertainty is " << sigastromcol << " \n";
	// Read observation file column holding the photometric uncertainty;
	status=readconfigint(instream1,&sigmagcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&sigmagcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding photometric uncertainty is " << sigmagcol << " \n";
	// Read observation file column holding the observatory code;
	status=readconfigint(instream1,&obscodecol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&obscodecol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding observatory code is " << obscodecol << " \n";
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
    } else if(string(argv[i]) == "-m" || string(argv[i]) == "-mjd" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjdstart" || string(argv[i]) == "--MJDstart" || string(argv[i]) == "--ModifiedJulianDay" || string(argv[i]) == "--modifiedjulianday") {
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
    }  else if(string(argv[i]) == "-statevec" || string(argv[i]) == "-sv" || string(argv[i]) == "-statevecs" || string(argv[i]) == "--startpos" || string(argv[i]) == "--statevector" || string(argv[i]) == "--statevectors" || string(argv[i]) == "--startingstatevectors") {
      if(i+6 < argc) {
	//There is still something to read;
	x=stold(argv[++i]);
	y=stold(argv[++i]);
	z=stold(argv[++i]);
	vx=stold(argv[++i]);
	vy=stold(argv[++i]);
	vz=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Starting state vector keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-statestep" || string(argv[i]) == "-simpstep" || string(argv[i]) == "-ss" || string(argv[i]) == "-simplexstep" || string(argv[i]) == "-simpside" || string(argv[i]) == "-simplexside" || string(argv[i]) == "--initialsimplexstep" || string(argv[i]) == "--startingsimplexside") {
      if(i+1 < argc) {
	//There is still something to read;
	simplex_side=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Initial simplex side length keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timescale" || string(argv[i]) == "-ts" || string(argv[i]) == "-tscale" || string(argv[i]) == "-chartime" || string(argv[i]) == "-tchar" || string(argv[i]) == "--timescale" || string(argv[i]) == "--characteristictime" || string(argv[i]) == "--chartimescale") {
      if(i+1 < argc) {
	//There is still something to read;
	timescale=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Timescale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-ftol") {
      if(i+1 < argc) {
	//There is still something to read;
	ftol=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "ftol keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-bumpnum") {
      if(i+1 < argc) {
	//There is still something to read;
	bumpnum=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "bumpnum keyword supplied with no corresponding argument\n";
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
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  cout.precision(17);  
  cout << "input configuration file " << configfile << "\n";
  cout << "input observation file " << obsfile << "\n";
  cout << "input reference MJD " << mjdstart << "\n";
  cout << "input state vector position " << x << " " << y << " " << z << " AU\n";
  cout << "input state vector velocity " << vx << " " << vy << " " << vz << " km/sec\n";
  cout << "side length for initial simplex " << simplex_side << "\n";
  cout << "time scale for converting velocities to distance units " << timescale << " days\n";
  cout << "MJD, RA, and Dec are in columns " << MJDcol << ", " << RAcol << ", " << Deccol << ", respectively\n";
  cout << "magnitude and astrometric uncertainty (arcsec) are in columns " << magcol << " and " << sigastromcol << "\n";
  cout << "magnitude uncertainty and observatory code are in columns " << sigmagcol << " and " << obscodecol << "\n";
  cout << "ftol is " << ftol << "\n";
  cout << "bumpnum is " << bumpnum << "\n";
  cout << "output file " << outfile << "\n";

  // Read input observation file.
  instream1.open(obsfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << obsfile << "\n";
    return(1);
  }
  badread=0;
  while(reachedeof==0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(lnfromfile.size()>20 && (reachedeof==0 || reachedeof==1)) {
      fieldnum = get_col_vector01(lnfromfile, linestringvec);
      if(fieldnum>=MJDcol) {
	ldval = stold(linestringvec[MJDcol-1]);
	obsMJD.push_back(ldval);
      } else badread=1;
      if(fieldnum>=RAcol) {
	dval = stod(linestringvec[RAcol-1]);
	obsRA.push_back(dval);
      } else badread=1;
      if(fieldnum>=Deccol) {
	dval = stod(linestringvec[Deccol-1]);
	obsDec.push_back(dval);
      } else badread=1;
      if(fieldnum>=magcol) {
	dval = stod(linestringvec[magcol-1]);
	magnitude.push_back(dval);
      } else badread=1;
      if(fieldnum>=sigastromcol) {
	dval = stod(linestringvec[sigastromcol-1]);
	sigastrom.push_back(dval);
      } else badread=1;
      if(fieldnum>=sigmagcol) {
	dval = stod(linestringvec[sigmagcol-1]);
	sigmag.push_back(dval);
      } else badread=1;
       if(fieldnum>=obscodecol) {
	 obscodevec.push_back(linestringvec[obscodecol-1]);
      } else badread=1;
    }
    cout << obsMJD.size() << " " << obsMJD[obsMJD.size()-1] << " " << obsRA[obsMJD.size()-1] << " " << obsDec[obsMJD.size()-1] << "\n";
  }
  obsnum=obsMJD.size();
  if(obsnum<2 || long(obsRA.size())!=obsnum || long(obsDec.size())!=obsnum || long(magnitude.size())!=obsnum || long(sigastrom.size())!=obsnum || long(sigmag.size())!=obsnum || long(obscodevec.size())!=obsnum) {
    cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << obsnum << " " << obsRA.size() << " " << obsDec.size() << " " << magnitude.size() << " " << sigastrom.size() << " " << sigmag.size() << " " << obscodevec.size() << "\n";
    return(1);
  }
  if(badread>=1){
    cerr << "Error reading " << obsfile << ": apparent short line\n";
    return(1);
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();

  // Read observatory code file
  read_obscode_file2(observatory_code_file, observatory_list,0);

  cout << "Read " << observatory_list.size() << " lines from observatory code file " << observatory_code_file << ":\n";
  
  if(DEBUG>=2) {
    for(i=0;i<long(observatory_list.size());i++) {
      cout << observatory_list[i].obscode << " " << observatory_list[i].obslon << " " << observatory_list[i].plxcos << " " << observatory_list[i].plxsin << "\n";
    }
  }
  
  // Calculate the exact position of the observer at the time of each observation.
  observerpos={};
  for(obsct=0;obsct<obsnum;obsct++) {
    if(obsct==0 || (obsct>0 && obscodevec[obsct]!=obscodevec[obsct-1])) {
      // Observatory has changed: get observatory coordinates for this image.
      stringncopy01(obscode, obscodevec[obsct],MINSTRINGLEN);
      status = obscode_lookup(observatory_list,obscode,obslon,plxcos,plxsin);
      if(status>0) {
	cerr << "ERROR: obscode_lookup failed for observatory code " << obscode << "\n";
	return(3);
      }
    }
    observer_barycoords01LD(obsMJD[obsct], 5, obslon, plxcos, plxsin, planetmjd, Earthpos, outpos);
    observerpos.push_back(outpos);
  }

  // Integrator wants units of km and km/sec.
  startpos = point3LD(x*AU_KM,y*AU_KM,z*AU_KM);
  startvel = point3LD(vx,vy,vz);
  obspos=obsvel={};
  integrate_orbit03LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, obsMJD, startpos, startvel, mjdstart, obspos, obsvel);
  
  for(obsct=0;obsct<obsnum;obsct++) {
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    stateunitLD_to_celestial(outpos, outRA, outDec);
    //cout  << fixed << setprecision(6) << "Input MJD " << mjdstart << ": " << x << " "  << y << " "  << z << " " << " Output: " << obsMJD[obsct] << ": " << obspos[obsct].x/AU_KM << " " << obspos[obsct].y/AU_KM << " " << obspos[obsct].z/AU_KM << "\n";
    cout  << fixed << setprecision(6) << "Input MJD " << obsMJD[obsct] << ": " << obsRA[obsct] << " "  << obsDec[obsct] << " "  << " Output: " << outRA << ": " << outDec <<  "\n";
  }
  
  scalestate[0] = x;
  scalestate[1] = y;
  scalestate[2] = z;
  scalestate[3] = vx*SOLARDAY*timescale/AU_KM;
  scalestate[4] = vy*SOLARDAY*timescale/AU_KM;
  scalestate[5] = vz*SOLARDAY*timescale/AU_KM;
  
  chisq = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,scalestate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
  cout << "Chi-square value for input state vector is " << chisq << "\n";

  // Create input simplex.
  simplex[0] = scalestate;
  for(i=1;i<=6;i++) {
    simplex[i] = scalestate;
    simplex[i][i-1] += simplex_side;
  }

  cout << "bumpnum = " << bumpnum << "\n";

  // Calculate chi-squared for all the points
  chivec[0] = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[0],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
  for(i=1;i<=6;i++) {
    chivec[i] = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[i],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
  }
  // Launch downhill simplex search
  simprange = LARGERR;
  globalchimin=LARGERR;
  while(simprange>ftol) {
    // Find best and worst points
    worstpoint=bestpoint=0;
    bestchi = worstchi = chivec[0];
    for(i=1;i<=6;i++) {
      if(chivec[i]<bestchi) {
      bestchi = chivec[i];
      bestpoint=i;
      }
      if(chivec[i]>worstchi) {
	worstchi = chivec[i];
	worstpoint=i;
      }
    }
    simprange = (worstchi-bestchi)/bestchi;
    cout << fixed << setprecision(6) << "Best chi-square value is " << bestchi << ", range is " << simprange << ", vector is " << simplex[bestpoint][0] << " "  << simplex[bestpoint][1] << " " << simplex[bestpoint][2] << " "   << simplex[bestpoint][3]*AU_KM/timescale/SOLARDAY << " "  << simplex[bestpoint][4]*AU_KM/timescale/SOLARDAY << " "  << simplex[bestpoint][5]*AU_KM/timescale/SOLARDAY << " bumpct " << bumpct << "\n";
      
    // Try to reflect away from worst point
    // Find mean over all the points except the worst one
    for(j=0;j<6;j++) scalestate[j]=0L;
    for(i=0;i<=6;i++) {
      if(i!=worstpoint) {
	for(j=0;j<6;j++) scalestate[j] += simplex[i][j];
      }
    }
    for(j=0;j<6;j++) scalestate[j]/=6.0L;
    // Calculate new trial point
    for(j=0;j<6;j++) newstate[j] = scalestate[j] - (simplex[worstpoint][j] - scalestate[j]);
    // Calculate chi-square value at this new point
    chisq = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,newstate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
    if(chisq<globalchimin) {
      globalchimin=chisq;
      bestRA=fitRA;
      bestDec=fitDec;
      bestresid=fitresid;
      beststatevec=newstate;
    }

    bumpct++;
    if(chisq<bestchi) {
      // Very good result. Let this point replace worstpoint in the simplex
      	for(j=0;j<6;j++) simplex[worstpoint][j] = newstate[j];
	chivec[worstpoint]=chisq;
      // Extrapolate further in this direction: maybe we can do even better
      for(j=0;j<6;j++) newstate[j] = scalestate[j] - 2.0L*(simplex[worstpoint][j] - scalestate[j]);
      newchi = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,newstate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
      if(newchi<globalchimin) {
	globalchimin=newchi;
	bestRA=fitRA;
	bestDec=fitDec;
	bestresid=fitresid;
	beststatevec=newstate;
      }
      bumpct++;
      if(newchi<chisq) {
	// Let this even better point replace worstpoint in the simplex
	for(j=0;j<6;j++) simplex[worstpoint][j] = newstate[j];
	chivec[worstpoint]=newchi;
      }
    } else {
      // How many points are worse than the new one?
      countworse=0;
      for(i=0;i<=6;i++) {
	if(chivec[i]>chisq) countworse++;
      }
      if(countworse>=2) { 
	// The new point is better than the previous worst and second worst
	// Add it to the simplex in place of the worst point
      	for(j=0;j<6;j++) simplex[worstpoint][j] = newstate[j];
	chivec[worstpoint]=chisq;
      } else {
	// This new point is not much good. Try going less far in the same direction.
	for(j=0;j<6;j++) newstate[j] = scalestate[j] - 0.5L*(simplex[worstpoint][j] - scalestate[j]);
	// Calculate chi-square value at this new point
	chisq = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,newstate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
	if(chisq<globalchimin) {
	  globalchimin=chisq;
	  bestRA=fitRA;
	  bestDec=fitDec;
	  bestresid=fitresid;
	  beststatevec=newstate;
	}
 	bumpct++;
	// How many points are worse than the new one?
	countworse=0;
	for(i=0;i<=6;i++) {
	  if(chivec[i]>chisq) countworse++;
	}
	if(countworse>=2) { 
	  // The new point is better than the previous worst and second worst
	  // Add it to the simplex in place of the worst point
	  for(j=0;j<6;j++) simplex[worstpoint][j] = newstate[j];
	  chivec[worstpoint]=chisq;
	} else {
	  // Still no good. Try contracting around the best point
	  for(i=0;i<=6;i++) {
	    if(i!=bestpoint) {
	      for(j=0;j<=6;j++) simplex[i][j] -= 0.5L*(simplex[i][j] - simplex[bestpoint][j]);
	      chivec[i] = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[i],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
	      if(chivec[i]<globalchimin) {
		globalchimin=chivec[i];
		bestRA=fitRA;
		bestDec=fitDec;
		bestresid=fitresid;
		beststatevec=simplex[i];
	      }
	      bumpct++;
	    }
	  }
	}
      }
    }
    if(bumpnum>0 && bumpct>=bumpnum && simprange<TOLSCALE*ftol) {
      // Bump up the size of the simplex, to help us explore
      // the parameter space faster and not get stuck.
      // Find the centroid of the current simplex
      for(j=0;j<6;j++) scalestate[j]=0L;
      for(i=0;i<=6;i++) {
	for(j=0;j<6;j++) scalestate[j] += simplex[i][j];
      }
      for(j=0;j<6;j++) scalestate[j]/=7.0L;
      // Pick three random dimensions in which to tweak the mean.
      j = 6*rand()/RAND_MAX;
      if(j>=0 && j<6) scalestate[j] += simplex_side*((0.5*rand()-RAND_MAX)/RAND_MAX);
      j = 6*rand()/RAND_MAX;
      if(j>=0 && j<6) scalestate[j] += simplex_side*((0.5*rand()-RAND_MAX)/RAND_MAX);
      j = 6*rand()/RAND_MAX;
      if(j>=0 && j<6) scalestate[j] += simplex_side*((0.5*rand()-RAND_MAX)/RAND_MAX);
      // Reconstruct simplex      
      simplex[0] = scalestate;
      simplex[0][0] -= simplex_side;
      for(i=1;i<=6;i++) {
	simplex[i] = scalestate;
	simplex[i][i-1] += simplex_side;
      }
      // Re-calculate all the chi-squared values
      for(i=0;i<=6;i++) {
	chivec[i] = tortoisechi01(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[i],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
	if(chivec[i]<globalchimin) {
	  globalchimin=chivec[i];
	  bestRA=fitRA;
	  bestDec=fitDec;
	  bestresid=fitresid;
	  beststatevec=simplex[i];
	}
      }
      bumpct=0; // Reset the count
    }
  }
  for(obsct=0;obsct<obsnum;obsct++)
    {
      cout << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << bestRA[obsct] << " " << bestDec[obsct] << " " << bestresid[obsct] << "\n";
    }
  cout << fixed << setprecision(6) << "Best chi-square value is " << globalchimin << ", range is " << simprange << ", vector is " << beststatevec[0] << " "  << beststatevec[1] << " " << beststatevec[2] << " "   << beststatevec[3]*AU_KM/timescale/SOLARDAY << " "  << beststatevec[4]*AU_KM/timescale/SOLARDAY << " "  << beststatevec[5]*AU_KM/timescale/SOLARDAY << "\n";

  return(0);
}


