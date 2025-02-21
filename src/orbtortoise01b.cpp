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

int Keplerint02(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, const vector <long double> &SunMJD, const vector <point3LD> &Sunpos, const vector <point3LD> &Sunvel, vector <point3LD> &obspos, vector <point3LD> &obsvel);

int integrate_orbit03LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <long double> &obsMJD, point3LD startpos, point3LD startvel, long double mjdstart, vector <point3LD> &obspos, vector <point3LD> &obsvel);

long double tortoisechi01(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <long double> &scalestate, long double timescale, long double mjdstart, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid);

long double tortoisechi02(const vector <long double> &SunMJD, const vector <point3LD> &Sunpos, const vector <point3LD> &Sunvel, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <long double> &scalestate, long double timescale, long double mjdstart, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid);

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
  vector <point3LD> Sunvel;
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
  observatory obs1 = observatory("I11",0l,0l,0l);
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
  make_LDvec(6,scalestate);
  vector <long double> newstate;
  make_LDvec(6,newstate);
  vector <long double> chivec;
  make_LDvec(7,chivec);
  vector <long double> beststatevec;
  long double chisq=0L;
  long double newchi=0L;
  long double bestchi=0L;
  long double worstchi=0L;
  int worstpoint=0;
  int bestpoint=0;
  vector <vector <long double>> simplex;
  make_LDmat(7,6,simplex);
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
	  status=status=readconfigstring(instream1,observatory_code_file);
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
	    for(j=0;j<planetmjd.size();j++) {
	      if(mjdtest[j]!=planetmjd[j]) {
		cout << "ERROR: time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(1);
	      }
	    }
	  }
	  for(j=0;j<temppos.size();j++) {
	    planetpos.push_back(temppos[j]);
	  }
	  if(planetct == pctEarth) Earthpos = temppos;
	  if(planetct == pctSun) Sunpos = temppos;
	  if(planetct == pctSun) Sunvel = tempvel;
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
  if(obsnum<2 || obsRA.size()!=obsnum || obsDec.size()!=obsnum || magnitude.size()!=obsnum || sigastrom.size()!=obsnum || sigmag.size()!=obsnum || obscodevec.size()!=obsnum) {
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
  instream1.open(observatory_code_file);
  if(!instream1) {
    cerr << "can't open input file " << observatory_code_file << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  while (!instream1.eof() && !instream1.fail() && !instream1.bad())
    {
      instream1 >> stest;
      stringncopy01(obscode,stest,MINSTRINGLEN);
      instream1 >> obslon;
      instream1 >> plxcos;
      instream1 >> plxsin;
      obs1 = observatory(obscode,obslon,plxcos,plxsin);
      observatory_list.push_back(obs1);
      // Skip the rest of the line
      getline(instream1,lnfromfile);
    }
  instream1.close();
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << observatory_code_file << ":\n";
  
  if(DEBUG>=2) {
    for(i=0;i<observatory_list.size();i++) {
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
  
  chisq = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,scalestate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
  
  cout << "Chi-square value for input state vector is " << chisq << "\n";

  // Create input simplex.
  simplex[0] = scalestate;
  for(i=1;i<=6;i++) {
    simplex[i] = scalestate;
    simplex[i][i-1] += simplex_side;
  }

  cout << "bumpnum = " << bumpnum << "\n";

  // Calculate chi-squared for all the points
  chivec[0] = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[0],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
  for(i=1;i<=6;i++) {
    chivec[i] = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[i],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
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
    chisq = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,newstate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
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
      newchi = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,newstate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
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
	chisq = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,newstate,timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
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
	      chivec[i] = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[i],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
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
	chivec[i] = tortoisechi02(planetmjd,Sunpos,Sunvel,observerpos,obsMJD,obsRA,obsDec,sigastrom,simplex[i],timescale,mjdstart,fitRA,fitDec,fitresid)/double(obsMJD.size());
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


// integrate_orbit03LD: April 08, 2022
// Uses modeling of the acceleration as a polynomial of order n>1
// to integrate the orbit of a massless test particle (e.g. asteroid)
// under the gravity of multiple 'planets'. It is assumed that
// in general these 'planets' will consist of the Sun, the
// eight major planets, and the Moon (possibly needed for
// cases of NEOs closely approaching the Earth). However,
// more or fewer planets may be used as desired.
// Note that the vector obsMJD is assumed to be time-sorted, and
// serious failures will result if it is not.
int integrate_orbit03LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <long double> &obsMJD, point3LD startpos, point3LD startvel, long double mjdstart, vector <point3LD> &obspos, vector <point3LD> &obsvel)
{
  vector <point3LD> planetsalltimes;
  vector <point3LD> planetsonce;
  vector <point3LD> targaccel;
  vector <point3LD> accelfit;
  vector <point3LD> targvel;
  vector <point3LD> targpos;
  vector <point3LD> accelmod;
  vector <long double> temptime;
  vector <long double> ppxvec;
  vector <long double> ppyvec;
  vector <long double> ppfitvec;
  vector <long double> obsMJD2 = obsMJD; // Mutable copy of immutable input vector.
  point3LD singleaccel = point3LD(0L,0L,0L);
  point3LD singlevel = point3LD(0L,0L,0L);
  point3LD singlepos = point3LD(0L,0L,0L);
  int i=0;
  int j=0;
  int k=0;
  int obsnum=obsMJD.size();
  int obsct=0;
  int endhere=-1;
  int planetpointnum = planetmjd.size();
  vector <long double> forwardmjd;
  vector <long double> backwardmjd;
   int planetpointct = 0;
  int pointafter=0;
  int pointbefore=0;
  int latestpoint=0;
  int stepsin=0;
  long double dt0=0L;
  long double dt1=0L;
  long double dt2=0L;
  long double timemult=0L;
  point3LD accelslope = point3LD(0L,0L,0L);
  int obsbefore=0;
  int obsafter=0;
  long double lastobstime=0L;
  long double firstobstime=0L;
  int obspoint=0;

  // Correct the input observed times, assumed to be UT1, by the
  // corretion to TT (terrestrial time) used, e.g. by JPL Horizons.
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD2[obsct] += TTDELTAT/SOLARDAY; // UT is always behind TT, hence UT must be adjusted forward.
  }
    
  firstobstime = obsMJD2[0];
  lastobstime = obsMJD2[obsnum-1];
  if(DEBUG>1) cout << "Integration will span times from MJD " << firstobstime << " to " << lastobstime << "\n";

  if(polyorder<2) {
    cerr << "ERROR: integrate_orbit03LD called with polyorder = " << polyorder << "\n";
    cerr << "polyorder must be at least 2!\n";
    return(1);
  }
  
  if(lastobstime<firstobstime) {
    cerr << "ERROR: integrate_orbit03LD called with end time (" << lastobstime << ") before start time (" << firstobstime << ")\n";
    return(1);
  } else if(mjdstart<=planetmjd[1] || firstobstime<=planetmjd[1] || lastobstime>=planetmjd[planetpointnum-1]) {
    cerr << "ERROR: integrate_orbit03LD called with start time " << mjdstart << " or time range )" << firstobstime << "-" << lastobstime << ") outside range of planet vectors (" << planetmjd[1] << "-" << planetmjd[planetpointnum-1] << ")\n";
    return(1);
  }
  // Make sure that relevant vectors for the polynomial fitting
  // are all large enough.
  for(i=0;i<=polyorder+1;i++) {
    targaccel.push_back(singleaccel);
    accelfit.push_back(singleaccel);
    targvel.push_back(singlevel);
    targpos.push_back(singlepos);
    accelmod.push_back(singleaccel);
    temptime.push_back(0L);
    ppfitvec.push_back(0L);
  }

  // Make sure the output vectors are large enough
  obspos={};
  obsvel={};
  for(obsct=0;obsct<=obsnum;obsct++) {
    obsvel.push_back(singlevel);
    obspos.push_back(singlepos);
  }
  // Are the observations before or after the starting time, or both?
  obsbefore=obsafter=0;
  for(obsct=0;obsct<=obsnum;obsct++) {
    if(obsMJD2[obsct]<mjdstart) obsbefore=1;
    if(obsMJD2[obsct]>=mjdstart) obsafter=1;
  }

  if(DEBUG>1) cout << "Checking for forward integration\n";
  if(obsafter==1) {
  if(DEBUG>1) cout << "Forward integration will be performed\n";
    // Integrate forward in time to find the position of
    // the target asteroid at all positions after starttime.

    // Load the initial time and planet position vectors
    planetsonce={};
    forwardmjd={};
    planetsalltimes={};
    temptime[0] = mjdstart;
    forwardmjd.push_back(mjdstart);
    nplanetpos01LD(temptime[0]-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
    for(i=0;i<planetnum;i++) planetsalltimes.push_back(planetsonce[i]);
    if(DEBUG>1) cout  << fixed << setprecision(6) << "mjdstart = " << mjdstart << " loaded: " << forwardmjd[0] << "\n";
 
    // Find the first observation time simultaneous with or after starttime
    obsct=0;
    while(obsMJD2[obsct]<mjdstart) obsct++;
    obspoint=obsct;
    
    // Find the first point in the planet files that is after mjdstart
    for(i=0;i<planetpointnum;i++) {
      if(planetmjd[i]>mjdstart) break;
    }
    pointafter = i; // first point after mjdstart

    if(DEBUG>1) cout  << fixed << setprecision(6) << "Starting points for forward integration:\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Observations, point " << obspoint << ", time " << obsMJD2[obspoint] << "\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Planet file, point " << pointafter << ", time " << planetmjd[pointafter] << "\n";
    
    // Load a vector with times and planet positions for all the planet file
    // points spanning the times of the observations.
    
    dt0 = (planetmjd[pointafter+1] - planetmjd[pointafter])*SOLARDAY/sqrt(M_PI);
    
    j=1;
    i=0;
    obsct=obspoint;
    // j counts elements actually loaded into forwardmjd and planetsalltimes
    // i counts steps in the planetfile past pointafter
    // obsct indexes the observation vector.
    for(obsct=obspoint;obsct<obsnum;obsct++) {
      while(planetmjd[pointafter+i]<obsMJD2[obsct]) {
	// Load any planet file points until we get past obsMJD2[obsct]
	forwardmjd.push_back(planetmjd[pointafter+i]);
	if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " > planetmjd[" << pointafter+i << "] = " << planetmjd[pointafter+i] << ": loaded forwardmjd[" << j << "] = " << forwardmjd[j] << "\n";
	planetsonce={};
	nplanetgrab01LD(pointafter+i, planetnum, planetmjd, planetpos, planetsonce);
	for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
	i++;
	j++;
      }
      // Load the next observation point.
      forwardmjd.push_back(obsMJD2[obsct]);
      planetsonce={};
      nplanetpos01LD(obsMJD2[obsct]-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      j++;
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " < planetmjd[" << pointafter+i << "] = " << planetmjd[pointafter+i] << ": loaded forwardmjd[" << j-1 << "] = " << forwardmjd[j-1] << "\n";
    }
    // Add extra planet points until we have polyorder+2
    while(j<polyorder+2) {
      forwardmjd.push_back(planetmjd[pointafter+i]);
      planetsonce={};
      nplanetgrab01LD(pointafter+i, planetnum, planetmjd, planetpos, planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " < planetmjd[" << pointafter+i << "] = " << planetmjd[pointafter+i] << ": loaded forwardmjd[" << j-1 << "] = " << forwardmjd[j-1] << "\n";
      i++;
      j++;
    }
    if(DEBUG>1) {
      cout << fixed << setprecision(6)  << "Loaded " << forwardmjd.size() << " points in forwardmjd\n";
      for(i=0;i<forwardmjd.size();i++) {
	cout  << fixed << setprecision(6) << "Forward MJD = " << forwardmjd[i] << "\n";
      }
    }
    // Now we've loaded all the times we want for the forward integration
    // Load temporary time vector
    for(j=0;j<=polyorder+1;j++) temptime[j] = forwardmjd[j];

    //cout << fixed << setprecision(6)  << "Loaded " << forwardmjd.size() << " points in forwardmjd\n";
    //cout << fixed << setprecision(6)  << "Loaded " << planetsalltimes.size() << " points in planetsalltimes\n";
    
    // Load starting position and velocity
    targvel[0] = startvel;
    targpos[0] = startpos;
    j=0;
    //cout << temptime[j] << " " << targpos[j].x << " " << targpos[j].y << " " << targpos[j].z << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
    // Bootstrap up to a fit of order polyorder.
    // Calculate acceleration at starting point, loading planet positions from big vector.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[0], targaccel[0]);
    dt1 = (temptime[1]-temptime[0])*SOLARDAY;
    // First Approx: estimate next position, assuming constant acceleration.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1;
    // Calculate acceleration at this new position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*1 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);
  
    // Second approx: linearly varying acceleration.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;

    // Improved position for next time step.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Re-calculate acceleration at this improved position.
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);

    // Re-calculate improved acceleration slope.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;
  
    // Improved velocity for next time step
    targvel[1].x = targvel[0].x + targaccel[0].x*dt1 + accelslope.x*0.5L*dt1*dt1;
    targvel[1].y = targvel[0].y + targaccel[0].y*dt1 + accelslope.y*0.5L*dt1*dt1;
    targvel[1].z = targvel[0].z + targaccel[0].z*dt1 + accelslope.z*0.5L*dt1*dt1;

    // Use linearly extrapolated acceleration to estimate position for
    // the next time step.
    dt1 = (temptime[2]-temptime[1])*SOLARDAY;
    targpos[2].x = targpos[1].x + targvel[1].x*dt1 + targaccel[1].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[2].y = targpos[1].y + targvel[1].y*dt1 + targaccel[1].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[2].z = targpos[1].z + targvel[1].z*dt1 + targaccel[1].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Calculate acceleration for this extrapolated position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*2 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[2], targaccel[2]);

    // Now we have three acceleration points: can load for a full polynomial fit.
    for(stepsin=3;stepsin<=polyorder+1;stepsin++) {
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      cout.precision(17);
  
      // Perform new fits to revised accelerations
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      for(i=0;i<stepsin;i++) {
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
    }

    // We are now set up for a full-order polynomial integration.
    // We have valid positions in targpos, targvel, and temptime
    // for indices from 0 to polyorder+1.
    // See if any of the target points have already been calculated
    j=0;
    obsct=obspoint;
    while(j<=polyorder+1 && obsct<obsnum) {
      if(obsMJD2[obsct]==temptime[j]) {
	obspos[obsct]=targpos[j];
	obsvel[obsct]=targvel[j];
	j++;
	obsct++;
      } else if(obsMJD2[obsct]<temptime[j]) obsct++;
      else if(temptime[j]<obsMJD2[obsct]) j++;
      else {
	cerr << "Impossible time comparison case: " << temptime[j] << " " << obsMJD2[obsct] << "\n";
	return(4);
      }
    }
    latestpoint=polyorder+1;
    // Proceed with the full polynomial integration.
    while(latestpoint<forwardmjd.size()-1) {
      latestpoint++;
      // Cycle the dynamical vectors
      for(i=0;i<polyorder+1;i++) {
	temptime[i] = temptime[i+1];
	targaccel[i] = targaccel[i+1];
	targvel[i] = targvel[i+1];
	targpos[i] = targpos[i+1];
      }
      // Load a new point into temptime
      temptime[polyorder+1] = forwardmjd[latestpoint];
      // Fit for acceleration
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	//cout << "Using planet point " << latestpoint+j-polyorder-1 << " Earth at x = " << planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + 3].x << "\n";
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }

      // Use these revised accelerations to re-do the fits
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
	//cout << "Using planet point " << latestpoint+j-polyorder-1 << " Earth at x = " << planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + 3].x << "\n";
      }
      j=polyorder+1;
      //cout << temptime[j] << " " << targpos[j].x << " " << targpos[j].y << " " << targpos[j].z << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
      // Load any target points that have newly been calculated,
      // or re-calculated
     j=0;
      obsct=obspoint;
       while(j<=polyorder+1 && obsct<obsnum) {
	if(obsMJD2[obsct]==temptime[j]) {
	  obspos[obsct]=targpos[j];
	  obsvel[obsct]=targvel[j];
	  j++;
	  obsct++;
	} else if(obsMJD2[obsct]<temptime[j]) obsct++;
	else if(temptime[j]<obsMJD2[obsct]) j++;
	else {
	  cerr << "Impossible time comparison case: " << temptime[j] << " " << obsMJD2[obsct] << "\n";
	  return(4);
	}
      }
      // We have now gone through two iterations of extrapolation
      // to predict the next acceleration point as accurately as possible.
      // The next step of the loop will move the extrapolated point back
      // by one step, and use it to start extrapolating a new point,
      // at the same time refining the former extrapolated points.
    }
  } // END OF FORWARD INTEGRATION
  if(DEBUG>1) {
    cout << "Results of forward integration:\n";
    for(j=0;j<temptime.size();j++) {
      cout  << fixed << setprecision(6) << temptime[j] << " " << targpos[j].x/AU_KM << " " << targpos[j].y/AU_KM << " " << targpos[j].z/AU_KM << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
    }
  }

  // NOW PEFORM BACKWARD INTEGRATION, IF NECESSARY
  if(DEBUG>1) cout << "Checking for backward integration\n";
  if(obsbefore==1) {
  if(DEBUG>1) cout << "Backward integration will be performed\n";
    // Integrate backward in time to find the position of
    // the target asteroid at all positions before starttime.

    // Load the initial time and planet position vectors
    planetsonce={};
    backwardmjd={};
    planetsalltimes={};
    temptime[0] = -mjdstart;
    backwardmjd.push_back(-mjdstart);
    nplanetpos01LD(mjdstart-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
    for(i=0;i<planetnum;i++) planetsalltimes.push_back(planetsonce[i]);
    if(DEBUG>1) cout  << fixed << setprecision(6) << "mjdstart = " << mjdstart << " loaded: " << backwardmjd[0] << "\n";

    // Find the first observation time before starttime
    obsct=obsnum-1;
    while(obsMJD2[obsct]>=mjdstart) obsct--;
    obspoint=obsct;
    
    // Find the first point in the planet files that is before mjdstart
    for(i=planetpointnum-1;i>0;i--) {
      if(planetmjd[i]<mjdstart) break;
    }
    pointbefore = i; // first point before mjdstart

    if(DEBUG>1) cout  << fixed << setprecision(6) << "Starting points for backward integration:\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Observations, point " << obspoint << ", time " << obsMJD2[obspoint] << "\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Planet file, point " << pointbefore << ", time " << planetmjd[pointbefore] << "\n";
    
    // Load a vector with times and planet positions for all the planet file
    // points spanning the times of the observations.
    
    // Time scaling factor, designed to avoid integers or near-zero values
    dt0 = (planetmjd[pointbefore+1] - planetmjd[pointbefore])*SOLARDAY/sqrt(M_PI);
    
    j=1;
    i=0;
    obsct=obspoint;
    // j counts elements actually loaded into forwardmjd and planetsalltimes
    // i counts steps in the planetfile before pointbefore
    // obsct indexes the observation vector.
    for(obsct=obspoint; obsct>=0; obsct--) {
      while(planetmjd[pointbefore-i]>obsMJD2[obsct]) {
	// Load any planet file points until we get to one before obsMJD2[obsct]
      	backwardmjd.push_back(-planetmjd[pointbefore-i]);
	if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " < planetmjd[" << pointbefore-i << "] = " << planetmjd[pointbefore-i] << ": loaded backwardmjd[" << j << "] = " << backwardmjd[j] << "\n";
	planetsonce={};
	nplanetgrab01LD(pointbefore-i, planetnum, planetmjd, planetpos, planetsonce);
	for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
	i++;
	j++;
      }
      // Load the next observation point.
      backwardmjd.push_back(-obsMJD2[obsct]);
      planetsonce={};
      nplanetpos01LD(obsMJD2[obsct]-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      j++;
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " > planetmjd[" << pointbefore-i << "] = " << planetmjd[pointbefore-i] << ": loaded backwardmjd[" << j-1 << "] = " << backwardmjd[j-1] << "\n";
    }
    // Add extra planet points until we have polyorder+2
    while(j<polyorder+2) {
      backwardmjd.push_back(-planetmjd[pointbefore-i]);
      planetsonce={};
      nplanetgrab01LD(pointbefore-i, planetnum, planetmjd, planetpos, planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " > planetmjd[" << pointbefore-i << "] = " << planetmjd[pointbefore-i] << ": loaded backwardmjd[" << j-1 << "] = " << backwardmjd[j-1] << "\n";
      i++;
      j++;
    }
    if(DEBUG>1) {
      cout << fixed << setprecision(6)  << "Loaded " << backwardmjd.size() << " points in backwardmjd\n";
      for(i=0;i<backwardmjd.size();i++) {
	cout  << fixed << setprecision(6) << "Backward MJD = " << backwardmjd[i] << "\n";
      }
    }

    // Now we've loaded all the times we want for the backward integration
    // Load temporary time vector
    for(j=0;j<=polyorder+1;j++) temptime[j] = backwardmjd[j];

    // Load starting position and velocity, negative velocity for backward integration
    targpos[0] = startpos;
    targvel[0].x = -startvel.x;
    targvel[0].y = -startvel.y;
    targvel[0].z = -startvel.z;
    // Bootstrap up to a fit of order polyorder.
    // Calculate acceleration at starting point, loading planet positions from big vector.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[0], targaccel[0]);
    dt1 = (temptime[1]-temptime[0])*SOLARDAY;
    // First Approx: estimate next position, assuming constant acceleration.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1;
    // Calculate acceleration at this new position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*1 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);
  
    // Second approx: linearly varying acceleration.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;

    // Improved position for next time step.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Re-calculate acceleration at this improved position.
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);

    // Re-calculate improved acceleration slope.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;
  
    // Improved velocity for next time step
    targvel[1].x = targvel[0].x + targaccel[0].x*dt1 + accelslope.x*0.5L*dt1*dt1;
    targvel[1].y = targvel[0].y + targaccel[0].y*dt1 + accelslope.y*0.5L*dt1*dt1;
    targvel[1].z = targvel[0].z + targaccel[0].z*dt1 + accelslope.z*0.5L*dt1*dt1;

    // Use linearly extrapolated acceleration to estimate position for
    // the next time step.
    dt1 = (temptime[2]-temptime[1])*SOLARDAY;
    targpos[2].x = targpos[1].x + targvel[1].x*dt1 + targaccel[1].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[2].y = targpos[1].y + targvel[1].y*dt1 + targaccel[1].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[2].z = targpos[1].z + targvel[1].z*dt1 + targaccel[1].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Calculate acceleration for this extrapolated position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*2 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[2], targaccel[2]);

    // Now we have three acceleration points: can load for a full polynomial fit.
    for(stepsin=3;stepsin<=polyorder+1;stepsin++) {
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      cout.precision(17);
  
      // Perform new fits to revised accelerations
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      for(i=0;i<stepsin;i++) {
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      
    }

    // We are now set up for a full-order polynomial integration.
    // We have valid positions in targpos, targvel, and temptime
    // for indices from 0 to polyorder+1.
    // See if any of the target points have already been calculated
    j=0;
    obsct=obspoint;
    while(j<=polyorder+1 && obsct>=0) {
      if(obsMJD2[obsct]==-temptime[j]) {
	obspos[obsct]=targpos[j];
	obsvel[obsct].x=-targvel[j].x;
	obsvel[obsct].y=-targvel[j].y;
	obsvel[obsct].z=-targvel[j].z;
	j++;
	obsct--;
      } else if(obsMJD2[obsct]>-temptime[j]) obsct--;
      else if(-temptime[j]>obsMJD2[obsct]) j++;
      else {
	cerr << "Impossible time comparison case: " << -temptime[j] << " " << obsMJD2[obsct] << "\n";
	return(4);
      }
    }
    latestpoint=polyorder+1;
    // Proceed with the full polynomial integration.
    while(latestpoint<backwardmjd.size()-1) {
      latestpoint++;
      // Cycle the dynamical vectors
      for(i=0;i<polyorder+1;i++) {
	temptime[i] = temptime[i+1];
	targaccel[i] = targaccel[i+1];
	targvel[i] = targvel[i+1];
	targpos[i] = targpos[i+1];
      }
      // Load a new point into temptime
      temptime[polyorder+1] = backwardmjd[latestpoint];
      // Fit for acceleration
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }

      // Use these revised accelerations to re-do the fits
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      j=polyorder+1;
      //cout << temptime[j] << " " << targpos[j].x << " " << targpos[j].y << " " << targpos[j].z << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
      // Load any target points that have newly been calculated,
      // or re-calculated
      j=0;
      obsct=obspoint;
      while(j<=polyorder+1 && obsct>=0) {
	if(obsMJD2[obsct]==-temptime[j]) {
	  obspos[obsct] = targpos[j];
	  obsvel[obsct].x = -targvel[j].x;
	  obsvel[obsct].y = -targvel[j].y;
	  obsvel[obsct].z = -targvel[j].z;
	  j++;
	  obsct--;
	} else if(obsMJD2[obsct]>-temptime[j]) obsct--;
	else if(-temptime[j]>obsMJD2[obsct]) j++;
	else {
	  cerr << "Impossible time comparison case: " << -temptime[j] << " " << obsMJD2[obsct] << "\n";
	  return(4);
	}
      }
      // We have now gone through two iterations of extrapolation
      // to predict the next acceleration point as accurately as possible.
      // The next step of the loop will move the extrapolated point back
      // by one step, and use it to start extrapolating a new point,
      // at the same time refining the former extrapolated points.
    }
  } // Closes the statement doing backward integration.
 
  return(0);
}

// tortoisechi01: April 11, 2022:
// Get chi-square value based on an input simplex point.
long double tortoisechi01(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <long double> &scalestate, long double timescale, long double mjdstart, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid)
{
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  int obsct;
  int obsnum = obsMJD.size();
  long double light_travel_time;
  point3LD outpos = point3LD(0,0,0);
  long double outRA=0L;
  long double outDec=0L;
  long double ldval=0L;
  long double chisq=0L;
  double dval;
  resid = fitRA = fitDec = {};
  
  // Input point scalestate is supposed to hold a 6-D state vector
  // in units of AU and AU/timescale, where timescale is a value in
  // days used to convert velocity to distance units

  // Convert to km and km/sec
  point3LD startpos = point3LD(scalestate[0]*AU_KM,scalestate[1]*AU_KM,scalestate[2]*AU_KM);
  point3LD startvel = point3LD(scalestate[3]*AU_KM/SOLARDAY/timescale,scalestate[4]*AU_KM/SOLARDAY/timescale,scalestate[5]*AU_KM/SOLARDAY/timescale);

  // Integrate orbit.
  integrate_orbit03LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, obsMJD, startpos, startvel, mjdstart, obspos, obsvel);
		  
  for(obsct=0;obsct<obsnum;obsct++) {
    // Initial approximation of the coordinates relative to the observer
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    // Initial approximation of the observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = obspos[obsct].x - light_travel_time*obsvel[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - light_travel_time*obsvel[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - light_travel_time*obsvel[obsct].z - observerpos[obsct].z;
    // Light-travel-time corrected observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, outRA, outDec);
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Input MJD " << obsMJD[obsct] << ": " << obsRA[obsct] << " "  << obsDec[obsct] << " "  << " Output: " << outRA << ": " << outDec <<  "\n";
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    fitRA.push_back(outRA);
    fitDec.push_back(outDec);
    resid.push_back(dval);
  }
  chisq=0.0L;
  for(obsct=0;obsct<obsnum;obsct++) {
    chisq += LDSQUARE(resid[obsct]/sigastrom[obsct]);
    if(DEBUG>0) cout << "Residual for point " << obsct << " is " << resid[obsct] << "\n";
  }
  return(chisq);
}

  
// Keplerint02: April 13, 2022: Like Keplerint, but does the
// calculation for a bunch of points simultaneously. Note that
// we assume the observation times are in UT1 but mjdstart is
// in TT.
// Description of ancestor program Keplerint:
// Integrate an orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Keplerint02(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, const vector <long double> &SunMJD, const vector <point3LD> &Sunpos, const vector <point3LD> &Sunvel, vector <point3LD> &obspos, vector <point3LD> &obsvel)
{
  long double e,E,a,lscalar,r0,v0,r1,v1;
  point3LD startpos2 = startpos; // Mutable copy of immutable input
  point3LD startvel2 = startvel; // Mutable copy of immutable input
  point3LD lvec = point3LD(0L,0L,0L);
  point3LD lunit = point3LD(0L,0L,0L);
  point3LD r0unit = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  point3LD SunposNow = point3LD(0L,0L,0L);
  point3LD SunvelNow = point3LD(0L,0L,0L);
  point3LD targpos = point3LD(0L,0L,0L);
  point3LD targvel = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  int status = 0;
  long double junkra,junkdec,sinev,thetav,v1ra,v1dec;
  junkra = junkdec = sinev = thetav = v1ra = v1dec = 0L;
  int obsct=0;
  int obsnum = obsMJD.size();

  // Get sun position at mjdstart, assuming it is in TT.
  planetpos01LD(mjdstart-TTDELTAT/SOLARDAY,3,SunMJD,Sunpos,SunposNow);
  planetpos01LD(mjdstart-TTDELTAT/SOLARDAY,3,SunMJD,Sunvel,SunvelNow);
  startpos2.x -= SunposNow.x;
  startpos2.y -= SunposNow.y;
  startpos2.z -= SunposNow.z;
  startvel2.x -= SunvelNow.x;
  startvel2.y -= SunvelNow.y;
  startvel2.z -= SunvelNow.z;
 
  // Calculate scalar input position
  r0 = sqrt(dotprod3LD(startpos2,startpos2));
  v0 = sqrt(dotprod3LD(startvel2,startvel2));
  
  // Calculate specific energy and angular momentum
  E = 0.5L*v0*v0 - MGsun/r0;
  lvec = crossprod3LD(startpos2,startvel2);
  lscalar = sqrt(dotprod3LD(lvec,lvec));
  if(E>=0L) {
    //cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5L/E;
  e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
  if(e<0L || e>=1.0L) {
    cerr << "ERROR: Keplerint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  //cout << "semimajor axis = " << a/AU_KM << " and eccentricity = " << e << "\n";
	       
  // Calculate angle theta0 from perihelion using simple ellipse geometry.
  if(e>0L) {
    costheta = ((a-a*e*e)/r0 - 1.0L)/e;
    if(costheta>=-1.0L && costheta<=1.0L) theta0 = acos(costheta);
    else {
      cerr << "ERROR: Keplerint finds costheta = " << costheta << "\n";
      return(1);
    }
  }
  radvel = dotprod3LD(startpos2,startvel2)/r0;
  //cout << "Radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: theta will be correct.
    //cout << "Moving outward from perihelion.\n";
    ;
  } else {
    // We are moving inward towards perihelion: theta needs adjustment.
    theta0 = 2.0L*M_PI - theta0;
    //cout << "Moving inward towards perihelion.\n";
  }
  //cout << "theta0 = " << theta0*DEGPRAD << "\n";
  
  // Calculate Goldstein's psi variable from theta.
  cospsi = (costheta + e)/(1.0L + costheta*e);
  if(cospsi>=-1.0L && cospsi<=1.0L) psi = acos(cospsi);
  else {
    cerr << "ERROR: Keplerint finds cospsi = " << cospsi << "\n";
    return(1);
  }
  if(radvel<0) {
    // We are moving inward towards perihelion: psi needs adjustment.
    psi = 2.0L*M_PI - psi;
  }
  //cout << "psi = " << psi*DEGPRAD << "\n";
 
  // Calculate time since perihelion using psi.
  omega = sqrt(MGsun/(a*a*a));
  //cout << "Period = " << 2.0L*M_PI/omega/SOLARDAY/365.25 << " years\n";
  t0omega = psi - e*sin(psi);
  //cout << "t0omega = " << t0omega;

  // Loop on all times-of-observation, and calculate the target position at those times
  obspos = obsvel = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    // The new time t1 for which we want to re-evaluate psi is
    // given by t0 + obsMJD[obsct]-mjdstart, with the TTDELTAT correction
    // applied to obsMJD assuming it is UT1.
    t1omega = t0omega + (obsMJD[obsct]+TTDELTAT/SOLARDAY-mjdstart)*SOLARDAY*omega;
    while(t1omega > 2.0L*M_PI) t1omega -= 2.0L*M_PI;
    while(t1omega < 0.0L) t1omega += 2.0L*M_PI;
    // Solve Kepler's equation for psi(t1)
    psi = kep_transcendental(t1omega,e,KEPTRANSTOL);
    cospsi = cos(psi);
    // Calculate theta(t1) from psi(t1)
    if(1.0L - e*cospsi != 0.0L) {
      costheta = (cospsi - e)/(1.0L - e*cospsi);
      if(costheta >= -1.0L && costheta <= 1.0L) theta1 = acos(costheta);
      else if(costheta < -1.0L) {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = M_PI;
      } else {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = 0.0L;
      }
      if(psi>M_PI && theta1<=M_PI) theta1 = 2.0L*M_PI - theta1;
    } else {
      cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0L - e*cospsi << "\n";
      theta1 = 0.0L;
    }
    while(theta1<0.0L) theta1 += 2.0L*M_PI;
    while(theta1>=2.0L*M_PI) theta1 -= 2.0L*M_PI;

    // Calculate r(t1) from psi(t1)
    r1 = a*(1.0L - e*cospsi);
    // Calculate v1 from r1 and the known energy
    v1 = sqrt((E +  MGsun/r1)*2.0L);
  
    // Use vector algebra to find the full vector r(t1).
    // This vector is perpendicular to lvec, and is angled by theta1-theta0
    // relative to startpos2.
    // Convert angular momentum vector to spherical coordinates
    celedeproj01LD(lvec,&lra,&ldec); // Note that output is in degrees.
    celedeproj01LD(startpos2,&r0ra,&r0dec); // Note that output is in degrees.
    // Transform the starting unit vector into a coordinate system with
    // the angular momentum vector at the pole, and the old pole at RA=0
    poleswitch01LD(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
    // Rotate starting unit vector around the angular momentum axis by
    // the calculated angle.
    newra += theta1-theta0;
    // The unit vector for the new position r1 is on the equator at this RA,
    // in the coordinate system that has the angular momentum vector at the pole.
    // Convert back to the original coordinate system.
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
    // Now for the velocity. If the velocity is at right angle to the vector r1,
    // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
    // of the angle between v1 and r1.

    sinev = lscalar/v1/r1;
    if(sinev>=1.0L) thetav = 0.5L*M_PI;
    else if(sinev<0.0L) {
      cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
      thetav = 0.0L;
    }
    else thetav = asin(sinev);
    if(theta1<=M_PI) {
      // Outward bound from perihelion.
      newra += thetav;
    } else {
      // Inward bound to perihelion
      newra += (M_PI - thetav);
    }
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

    r1unit = celeproj01LD(r1ra*DEGPRAD,r1dec*DEGPRAD);
    v1unit =celeproj01LD(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
    targpos.x = r1unit.x*r1;
    targpos.y = r1unit.y*r1;
    targpos.z = r1unit.z*r1;
    targvel.x = v1unit.x*v1;
    targvel.y = v1unit.y*v1;
    targvel.z = v1unit.z*v1;
    // Add back in the Sun's position at time obsMJD[obsct]
    planetpos01LD(obsMJD[obsct],3,SunMJD,Sunpos,SunposNow);
    planetpos01LD(obsMJD[obsct],3,SunMJD,Sunvel,SunvelNow);
    targpos.x += SunposNow.x;
    targpos.y += SunposNow.y;
    targpos.z += SunposNow.z;
    targvel.x += SunvelNow.x;
    targvel.y += SunvelNow.y;
    targvel.z += SunvelNow.z;
    obspos.push_back(targpos);
    obsvel.push_back(targvel);
  }
    
  return(0);
}


	


// tortoisechi02: April 13, 2022:
// Get chi-square value based on an input simplex point.
// Like tortoisechi01, but uses Keplerian integration
// rather than n-body.
long double tortoisechi02(const vector <long double> &SunMJD, const vector <point3LD> &Sunpos, const vector <point3LD> &Sunvel, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <long double> &scalestate, long double timescale, long double mjdstart, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid)
{
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  int obsct;
  int obsnum = obsMJD.size();
  long double light_travel_time;
  point3LD outpos = point3LD(0,0,0);
  long double outRA=0L;
  long double outDec=0L;
  long double ldval=0L;
  long double chisq=0L;
  double dval;
  resid = fitRA = fitDec = {};
  int status=0;
  
  // Input point scalestate is supposed to hold a 6-D state vector
  // in units of AU and AU/timescale, where timescale is a value in
  // days used to convert velocity to distance units

  // Convert to km and km/sec
  point3LD startpos = point3LD(scalestate[0]*AU_KM,scalestate[1]*AU_KM,scalestate[2]*AU_KM);
  point3LD startvel = point3LD(scalestate[3]*AU_KM/SOLARDAY/timescale,scalestate[4]*AU_KM/SOLARDAY/timescale,scalestate[5]*AU_KM/SOLARDAY/timescale);

  // Integrate orbit.
  status=0;
  status = Keplerint02(GMSUN_KM3_SEC2,mjdstart,obsMJD,startpos,startvel,SunMJD,Sunpos,Sunvel,obspos,obsvel);
  if(status!=0) {
    // Keplerint02 failed, likely because input state vectors lead
    // to an unbound orbit.
    return(LARGERR);
  }
  for(obsct=0;obsct<obsnum;obsct++) {
    // Initial approximation of the coordinates relative to the observer
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    // Initial approximation of the observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = obspos[obsct].x - light_travel_time*obsvel[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - light_travel_time*obsvel[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - light_travel_time*obsvel[obsct].z - observerpos[obsct].z;
    // Light-travel-time corrected observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, outRA, outDec);
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Input MJD " << obsMJD[obsct] << ": " << obsRA[obsct] << " "  << obsDec[obsct] << " "  << " Output: " << outRA << ": " << outDec <<  "\n";
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    fitRA.push_back(outRA);
    fitDec.push_back(outDec);
    resid.push_back(dval);
  }
  chisq=0.0L;
  for(obsct=0;obsct<obsnum;obsct++) {
    chisq += LDSQUARE(resid[obsct]/sigastrom[obsct]);
    if(DEBUG>0) cout << "Residual for point " << obsct << " is " << resid[obsct] << "\n";
  }
  return(chisq);
}
