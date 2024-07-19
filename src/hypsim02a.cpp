/// June 10, 2024: hypsim02a.cpp
// Given a set of heliolinc hypotheses, generate simulated asteroids
// exactly matching the hypotheses, with angular positions and velocities
// randomly distributed on the celestial sphere.
// This version differs from hypsim01a.cpp in that it calculates magnitudes,
// including both phase and distance effects. This requires an input power
// law for the H magnitudes. The objective is to produce a more realistic
// distribution of solar elongations.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define DEBUG 0
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define IDNUMLEN 7
#define SHORTVECTOR 0.3 // Threshold for random vector being 'short'
#define PHASE_G 0.15

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
  cerr << "Usage: hypsim02a -cfg configfile -ranseed random_number_seed -times timefile -hyp hypothesis_file -unbound 0 -mjdref mjdref -mjdstart mjdstart -mjdend mjdend -astromerr 1-D astrometric error (arcsec) -maglim limiting_magnitude -Hmin min_Hmag -Hmax max_Hmag -Hslope Hmag_powerlaw_slope -simnum simnum -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile, timefile, hypfile, planetfile, outfile;
  long double mjdref = 0.0L;
  long double mjdstart = 0.0L;
  long double mjdend = 0.0L;
  long double GMsun = 0.0L;
  long simnum=0;
  long simct=0;
  int configread=0;
  int polyorder=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
  int unbound = 0;
  long double ldval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> mjdobs;
  vector <long double> hypdist;
  vector <long double> hypvel;
  vector <long double> hypacc;
  vector <long double> zenithRA;
  vector <long double> zenithDec;
  
  long double hdist,hvel,hacc;
  long double localg,tanvel,physacc,vesc,n;
  localg=tanvel=physacc=vesc=n=0.0L;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> targpos;
  vector <point3LD> targvel;
  point3LD targ_to_sun = point3LD(0,0,0);
  point3LD targvel_to_sunvel = point3LD(0,0,0);
  point3LD targ_to_obs = point3LD(0,0,0);
  point3LD obs_to_targ = point3LD(0,0,0);
  point3LD obs_to_sun = point3LD(0,0,0);
  point3LD Sunposnow = point3LD(0,0,0);
  point3LD Sunvelnow = point3LD(0,0,0);
  long double obsdist=1.0L;
  long double sundist=1.0L;
  long double obsfromsun=1.0L;
  long double light_travel_time;
  vector <point3LD> observerpos;
  vector <long double> targMJD;
  long timenum, timect, hypnum, hypct;
  timenum = timect = hypnum = hypct = 0;
  string lnfromfile;
  vector <string> linestringvec;
  //  double obslon = 289.26345L;
  //  double plxcos = 0.865020L;
  //  double plxsin = -0.500901L;
  double obslon = 289.25058l; // Changed to match X05 obscode
  double plxcos = 0.864981l;
  double plxsin = -0.500958l;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  point3LD lttpos = point3LD(0,0,0);
  point3LD temp3D = point3LD(0,0,0);
  point3LD unitZ = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  long double outRA,outDec;
  long double a,e,incl,phaseang,cosphase,sunelong;
  ofstream outstream1;
  long double astromsigma=0.1L;
  string seedstring;
  int goodsimct=0;
  int idnl = IDNUMLEN;
  int compfac = 1;
  string idnumstring;
  double minalt = 20.0l;
  double Hmin = 15.0;
  double Hmax = 20.0;
  double Hslope = 0.4;
  double limiting_mag = 0.0l;
  double x,xmin,xmax,acoef,absmag,obsmag,phi1,phi2;
  long double phaseslope = PHASE_G;

  if(argc<15) {
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
	  } else {
	    cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
	    if(planetct==pctSun) GMsun = planetmasses[planetct];
	  }
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
	  if(planetct == pctSun) {
	    Sunpos = temppos;
	    Sunvel = tempvel;
	  }
	  cout << "Finished reading ephemeris file " << planetfile << "\n";
	}
	// Read default reference MJD
	status=readconfigLD(instream1,&mjdref);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdref);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default reference MJD read as " << mjdref << "\n";
	// Read default starting MJD
	status=readconfigLD(instream1,&mjdstart);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdstart);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default starting MJD read as " << mjdstart << "\n";
	// Read default ending MJD
	status=readconfigLD(instream1,&mjdend);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdend);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default ending MJD read as " << mjdend << "\n";
	// Read astrometric error in arcseconds
	status=readconfigLD(instream1,&astromsigma);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&astromsigma);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Astrometric error read as " << astromsigma << " arcsec\n";
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
    } else if(string(argv[i]) == "-t" || string(argv[i]) == "-time" || string(argv[i]) == "-times" || string(argv[i]) == "-timefile" || string(argv[i]) == "--timefile" || string(argv[i]) == "--times" || string(argv[i]) == "--time" || string(argv[i]) == "--timesfile") {
      if(i+1 < argc) {
	// There is still something to read;
	timefile=argv[++i];
	i++;
      } else {
	cerr << "Time file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-h" || string(argv[i]) == "-hyp" || string(argv[i]) == "-hypfile" || string(argv[i]) == "--hypothesisfile" || string(argv[i]) == "--hypfile") {
      if(i+1 < argc) {
	// There is still something to read;
	hypfile=argv[++i];
	i++;
      } else {
	cerr << "Hypothesis file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-unbound" || string(argv[i]) == "--unbound") {
      if(i+1 < argc) {
	// There is still something to read;
	unbound=stoi(argv[++i]);
	i++;
      } else {
	cerr << "Unbound keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdref") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Reference MJD keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdstart" || string(argv[i]) == "-ms" || string(argv[i]) == "-MJDstart" || string(argv[i]) == "--mjdstart" || string(argv[i]) == "--MJDstart" || string(argv[i]) == "--ModifiedJulianDaystart" || string(argv[i]) == "--modifiedjuliandaystart") {
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
    } else if(string(argv[i]) == "-astromerr" || string(argv[i]) == "-astromsig" || string(argv[i]) == "-as"  || string(argv[i]) == "-sigast"  || string(argv[i]) == "-sigastrom" || string(argv[i]) == "--astromerror" || string(argv[i]) == "--astromsigma") {
      if(i+1 < argc) {
	//There is still something to read;
	astromsigma=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Astrometric error keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minalt") {
      if(i+1 < argc) {
	//There is still something to read;
	minalt=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum altitude keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maglim") {
      if(i+1 < argc) {
	//There is still something to read;
	limiting_mag=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Limiting magnitude keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-Hmin") {
      if(i+1 < argc) {
	//There is still something to read;
	Hmin=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Min Hmag keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-Hmax") {
      if(i+1 < argc) {
	//There is still something to read;
	Hmax=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Max Hmag keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-Hslope") {
      if(i+1 < argc) {
	//There is still something to read;
	Hslope=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Absolute magnitude slope keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-simnum" || string(argv[i]) == "-sn" || string(argv[i]) == "-snum" || string(argv[i]) == "--simnum") {
      if(i+1 < argc) {
        //There is still something to read;
        simnum=stol(argv[++i]);
        i++;
      }
      else {
        cerr << "Simulation number keyword supplied with too few corresponding arguments\n";
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
  cout << "randum number seed string " << seedstring << "\n";
  cout << "input configuration file " << configfile << "\n";
  cout << "input time file " << timefile << "\n";
  cout << "input hypothesis file " << hypfile << "\n";
  cout << "input starting MJD " << mjdstart << "\n";
  cout << "input ending MJD " << mjdstart << "\n";
  cout << "1-D Gaussian error added to output astrometry: " << astromsigma << " arcsec\n";
  cout << "number of encounters to simulate per hypothesis: " << simnum << "\n";
  cout << "output file " << outfile << "\n";

  // Read input time file
  mjdobs = {};
  instream1.open(timefile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << timefile << "\n";
    return(1);
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    ldval=0.0L;
    instream1 >> ldval;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad() && ldval>mjdstart && ldval<mjdend) mjdobs.push_back(ldval);
  }
  instream1.close();
  timenum = mjdobs.size();
  if(timenum<=0) {
    cout << "ERROR: no times found in input time file\n";
    return(1);
  }
  cout << "Read " << timenum << " valid times from input file " << timefile << "\n";
  
  // Read input hypothesis file
  hypdist = hypvel = hypacc = {};
  instream1.open(hypfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << hypfile << "\n";
    return(1);
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    instream1 >> hdist;
    instream1 >> hvel;
    instream1 >> hacc;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      localg = GMsun/LDSQUARE(hdist*AU_KM);
      physacc = -localg*hacc;
      vesc = 2.0*GMsun/hdist/AU_KM; // this is the square of the escape velocity
      tanvel = hdist*AU_KM*(physacc+localg); // this is the square of the tangential velocity
      if(tanvel<0.0l) {
	cerr << "ERROR: hypothesis point " << hdist << ", " << hvel << ", " << hacc << " is not possible for any trajectory\n";
	return(1);
      }
      if(vesc < tanvel + LDSQUARE(hvel*AU_KM/SOLARDAY) && unbound<=0) {
	cerr << "ERROR: hypothesis point " << hdist << ", " << hvel << ", " << hacc << " is not possible for a bound orbit\n";
	cerr << "vesc = " << sqrt(vesc) << ", tanvel = " << sqrt(tanvel) << ", hvel = " << hvel*AU_KM/SOLARDAY << ", totvel = " << sqrt(tanvel + LDSQUARE(hvel*AU_KM/SOLARDAY)) << "\n";
	return(1);
      }
      hypdist.push_back(hdist);
      hypvel.push_back(hvel);
      hypacc.push_back(hacc);    
    }
  }
  instream1.close();
  hypnum = hypdist.size();
  if(hypnum<=0) {
    cerr << "ERROR: no valid lines read from hypothesis file " << hypfile << "\n";
    return(3);
  }
  cout << "Read " << hypnum << " valid lines from input file " << hypfile << "\n";
 
  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine


  // Calculate the exact position of the observer at the time of each image.
  observerpos = {};
  zenithRA = zenithDec = {};
  for(timect=0; timect<timenum; timect++) {
    observer_barycoords01LD(mjdobs[timect], 5, obslon, plxcos, plxsin, planetmjd, Earthpos, outpos);
    observerpos.push_back(outpos);
    long double RA, Dec;
    zenith_radecLD(mjdobs[timect], obslon, plxcos, plxsin, RA, Dec);
    zenithRA.push_back(RA);
    zenithDec.push_back(Dec);
  }

  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "stringID a e incl q Q mjdref n absmag sundist obsdist sunelong phaseang mjdobs outRA outDec mag alt Ast-Sun(J2000x)(km) Ast-Sun(J2000y)(km) Ast-Sun(J2000z)(km) Ast-Sun(J2000vx)(km/s) Ast-Sun(J2000vy)(km/s) Ast-Sun(J2000vz)(km/s) perfectRA perfectDec\n";

  // Keep track of leading zeros for string IDs
  idnl = IDNUMLEN-1;
  compfac = 10;
  idnumstring={};
  for(i=0; i<idnl; i++) idnumstring.push_back('0');

  // Loop on hypotheses
  for(hypct=0; hypct<hypnum; hypct++) {
    // Loop on simulated objects.
    simct=0;
    while(simct<simnum) {
      // Randomly assign an absolute magnitude H
      acoef = 1.0l/(Hslope*log(10.0));
      xmin = exp(Hmin/acoef);
      xmax = exp(Hmax/acoef);
      x = xmin + (xmax-xmin)*unitvar(generator);
      absmag = acoef*log(x);
      // Obtain a random vector on the unit sphere
      int inunitsphere=0;
      while(!inunitsphere) {
	startpos.x = 2.0L*unitvar(generator) - 1.0L;
	startpos.y = 2.0L*unitvar(generator) - 1.0L;
	startpos.z = 2.0L*unitvar(generator) - 1.0L;
	ldval = vecabs3LD(startpos);
	if(ldval<=1.0L && ldval>SHORTVECTOR) inunitsphere=1;
      }
      vecnorm3LD(startpos);
      // Obtain a random perpendicular vector, for the tangential velocity
      int perpendicularcomp=0;
      while(!perpendicularcomp) {
	temp3D.x = 2.0L*unitvar(generator) - 1.0L;
	temp3D.y = 2.0L*unitvar(generator) - 1.0L;
	temp3D.z = 2.0L*unitvar(generator) - 1.0L;
	ldval = vecabs3LD(temp3D);
	if(ldval<=1.0L && ldval>SHORTVECTOR) {
	  // It's inside the unit sphere, and not too short
	  vecnorm3LD(temp3D);
	  ldval = dotprod3LD(startpos,temp3D);
	  if(ldval<SHORTVECTOR) perpendicularcomp=1; // It's reasonably close to perpendicular
	                                             // to the position vector startpos.
	}
      }
      outvel = crossprod3LD(temp3D,startpos); // Randomly oriented, perpendicular to startpos
      vecnorm3LD(outvel); // Unit vector for the tangential velocity
      // Calculate the tangential velocity in km/sec
      localg = GMsun/LDSQUARE(hypdist[hypct]*AU_KM); // Local gravitational acceleration in km/sec^2
      physacc = -localg*hypacc[hypct]; // r-double-dot as specified by hypothesis, in km/sec^2
      tanvel = hypdist[hypct]*AU_KM*(physacc+localg); // this is the square of the tangential velocity
      if(tanvel<0.0L) tanvel=0.0L; //Earlier check is supposed to have caught this: we are only hedging against roundoff here.
      else tanvel = sqrt(tanvel);
      // Convert radial velocity to km/sec, and distance to km
      hvel = hypvel[hypct]*AU_KM/SOLARDAY;
      hdist = hypdist[hypct]*AU_KM;
      // Sum the tangential and radial components to get the full velocity.
      startvel.x = outvel.x*tanvel + startpos.x*hvel;
      startvel.y = outvel.y*tanvel + startpos.y*hvel;
      startvel.z = outvel.z*tanvel + startpos.z*hvel;
      // Above 3-D velocity is guaranteed to produce the specified instantaneous values of rdot and r-double-dot
      // Scale the distance
      startpos.x *= hdist;
      startpos.y *= hdist;
      startpos.z *= hdist;

      // Calculate Keplerian semimajor axis and eccentricity
      long double r0 = vecabs3LD(startpos);
      long double v0 = vecabs3LD(startvel);
      long double u = dotprod3LD(startvel,startpos);
      a = r0*GMsun/(2.0l*GMsun-v0*v0*r0);
      long double alpha = GMsun/a;
      long double EC, ES, CH, SH;
      EC = ES = CH = SH = e = 0.0l;
  
      if(alpha>0.0l) {
	// Bound, elliptical orbit
	n = sqrt(GMsun/a/a/a);
	EC = 1.0l - r0/a;
	ES = u/n/a/a;
	e = sqrt(EC*EC + ES*ES);
      } else if (alpha<0.0l) {
	// Unbound, hyperbolic orbit
	CH = 1.0l - r0/a;
	SH = u/sqrt(-GMsun*a);
	e = sqrt(CH*CH - SH*SH);
	n=0.0L;
      }
      // Calculate Keplerian inclination
      temp3D = crossprod3LD(startpos,startvel);
      vecnorm3LD(temp3D);
      unitZ = point3LD(0.0L,0.0L,1.0L);
      ldval = dotprod3LD(unitZ,temp3D);
      if(ldval<=-1.0L) incl = 180.0L;
      else if(ldval>=1.0) incl = 0.0L;
      else incl = DEGPRAD*acos(ldval);
      
      // Add back in the Solar coordinates at mjdref
      planetposvel01LD(mjdref,polyorder,planetmjd,Sunpos,Sunvel,Sunposnow,Sunvelnow);
      startpos.x += Sunposnow.x;
      startpos.y += Sunposnow.y;
      startpos.z += Sunposnow.z;
      startvel.x += Sunvelnow.x;
      startvel.y += Sunvelnow.y;
      startvel.z += Sunvelnow.z;
      
      // Perform n-body integration starting from the reference time
      // with the state vectors produced above. Note that integrate_orbit3D
      // can handle integration backward in time as well as forward.
      integrate_orbit03LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, mjdobs, startpos, startvel, mjdref, targpos, targvel);
      int goodobject=1;
      // Initial test integration
      for(timect=0;timect<timenum;timect++) {
	planetposvel01LD(mjdobs[timect],polyorder,planetmjd,Sunpos,Sunvel,Sunposnow,Sunvelnow);
	// Calculate approximate target-to-object distance
	targ_to_obs.x = observerpos[timect].x - targpos[timect].x;
	targ_to_obs.y = observerpos[timect].y - targpos[timect].y;
	targ_to_obs.z = observerpos[timect].z - targpos[timect].z;
	obsdist = vecabs3LD(targ_to_obs);
	// Calculate and apply a light-travel-time correction
	light_travel_time = obsdist*1000.0/CLIGHT; // Factor of 1000 converts obsdist to meters
	lttpos.x = targpos[timect].x - light_travel_time*targvel[timect].x;
	lttpos.y = targpos[timect].y - light_travel_time*targvel[timect].y;
	lttpos.z = targpos[timect].z - light_travel_time*targvel[timect].z;
	// lttpos is now the target position with the light-travel-time correction applied
	targ_to_obs.x = observerpos[timect].x - lttpos.x;
	targ_to_obs.y = observerpos[timect].y - lttpos.y;
	targ_to_obs.z = observerpos[timect].z - lttpos.z;
	targ_to_sun.x = Sunposnow.x - lttpos.x;
	targ_to_sun.y = Sunposnow.y - lttpos.y;
	targ_to_sun.z = Sunposnow.z - lttpos.z;
	targvel_to_sunvel.x = Sunvelnow.x - targvel[timect].x;
	targvel_to_sunvel.y = Sunvelnow.y - targvel[timect].y;
	targvel_to_sunvel.z = Sunvelnow.z - targvel[timect].z;
	obs_to_targ.x = -targ_to_obs.x;
	obs_to_targ.y = -targ_to_obs.y;
	obs_to_targ.z = -targ_to_obs.z;
	
	obsdist = vecabs3LD(targ_to_obs);
	sundist = vecabs3LD(targ_to_sun);

	// Calculate observed magnitude not accounting for phase
	obsmag = absmag + 2.5*log10(obsdist*obsdist*sundist*sundist/AU_KM/AU_KM/AU_KM/AU_KM);
      
	// Calculate obs_to_targ unit vector
	vecnorm3LD(obs_to_targ);
	// Project onto the celestial sphere.
	stateunitLD_to_celestial(obs_to_targ, outRA, outDec);
	// Find out if object is above the horizon.
	double dist = distradec01(outRA, outDec, zenithRA[timect], zenithDec[timect]);
	if(obsmag>limiting_mag || dist > 90.0l-minalt) {
	  // Object is too faint for detection, or below the functional horizon
	  goodobject=0;
	} else {
	  // Correct the magnitude calculated above for phase effects
	  // Calculate phase angle.
	  cosphase = dotprod3LD(targ_to_obs,targ_to_sun)/obsdist/sundist;
	  if(cosphase>1.0L) {
	    cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
	    phaseang=0.0L;
	  } else if(cosphase<-1.0L) {
	    cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
	    phaseang=M_PI;
	  } else phaseang = acos(cosphase);
	  phi1 = exp(-PHASECONST_A1*pow(tan(phaseang/2.0),PHASECONST_B1));
	  phi2 = exp(-PHASECONST_A2*pow(tan(phaseang/2.0),PHASECONST_B2));
	  //cout << "obsmag = " << obsmag << " ";
	  obsmag -= 2.5*log10((1.0-phaseslope)*phi1 + phaseslope*phi2);
	  //cout << "becomes obsmag = " << obsmag << "\n";
	  if(obsmag>limiting_mag) goodobject=0; // Too faint, due to phase effects
	}
      }
      if(goodobject==1) {
	// Integration where output will occur
	for(timect=0;timect<timenum;timect++) {
	  planetposvel01LD(mjdobs[timect],polyorder,planetmjd,Sunpos,Sunvel,Sunposnow,Sunvelnow);
	  // Calculate approximate target-to-object distance
	  targ_to_obs.x = observerpos[timect].x - targpos[timect].x;
	  targ_to_obs.y = observerpos[timect].y - targpos[timect].y;
	  targ_to_obs.z = observerpos[timect].z - targpos[timect].z;
	  obsdist = vecabs3LD(targ_to_obs);
	  // Calculate and apply a light-travel-time correction
	  light_travel_time = obsdist*1000.0/CLIGHT; // Factor of 1000 converts obsdist to meters
	  lttpos.x = targpos[timect].x - light_travel_time*targvel[timect].x;
	  lttpos.y = targpos[timect].y - light_travel_time*targvel[timect].y;
	  lttpos.z = targpos[timect].z - light_travel_time*targvel[timect].z;
	  // lttpos is now the target position with the light-travel-time correction applied
	  targ_to_obs.x = observerpos[timect].x - lttpos.x;
	  targ_to_obs.y = observerpos[timect].y - lttpos.y;
	  targ_to_obs.z = observerpos[timect].z - lttpos.z;
	  targ_to_sun.x = Sunposnow.x - lttpos.x;
	  targ_to_sun.y = Sunposnow.y - lttpos.y;
	  targ_to_sun.z = Sunposnow.z - lttpos.z;
	  targvel_to_sunvel.x = Sunvelnow.x - targvel[timect].x;
	  targvel_to_sunvel.y = Sunvelnow.y - targvel[timect].y;
	  targvel_to_sunvel.z = Sunvelnow.z - targvel[timect].z;
	  obs_to_targ.x = -targ_to_obs.x;
	  obs_to_targ.y = -targ_to_obs.y;
	  obs_to_targ.z = -targ_to_obs.z;
	
	  obsdist = vecabs3LD(targ_to_obs);
	  sundist = vecabs3LD(targ_to_sun);
	  
	  // Calculate observed magnitude not accounting for phase
	  obsmag = absmag + 2.5*log10(obsdist*obsdist*sundist*sundist/AU_KM/AU_KM/AU_KM/AU_KM);
	
	  // Calculate obs_to_targ unit vector
	  vecnorm3LD(obs_to_targ);
	  // Project onto the celestial sphere.
	  stateunitLD_to_celestial(obs_to_targ, outRA, outDec);
	  // Calculate the object's altitude
	  double alt = 90.0l - distradec01(outRA, outDec, zenithRA[timect], zenithDec[timect]);
	  // Calculate phase angle.
	  cosphase = dotprod3LD(targ_to_obs,targ_to_sun)/obsdist/sundist;
	  if(cosphase>1.0L) {
	    cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
	    phaseang=0.0L;
	  } else if(cosphase<-1.0L) {
	    cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
	    phaseang=M_PI;
	  } else phaseang = acos(cosphase);
	  phi1 = exp(-PHASECONST_A1*pow(tan(phaseang/2.0),PHASECONST_B1));
	  phi2 = exp(-PHASECONST_A2*pow(tan(phaseang/2.0),PHASECONST_B2));
	  cout << "obsmag = " << obsmag << " ";
	  obsmag -= 2.5*log10((1.0-phaseslope)*phi1 + phaseslope*phi2);
	  cout << "becomes obsmag = " << obsmag << "\n";
	  // Also calculate solar elongation, while we have the values handy
	  obs_to_sun.x = Sunposnow.x - observerpos[timect].x;
	  obs_to_sun.y = Sunposnow.y - observerpos[timect].y;
	  obs_to_sun.z = Sunposnow.z - observerpos[timect].z;
	  obsfromsun = vecabs3LD(obs_to_sun);
	  double q = a*(1.0L-e);
	  double Q = a*(1.0L+e);
	  cosphase = dotprod3LD(obs_to_targ,obs_to_sun)/obsfromsun; // obs_to_targ was normalized already.
	  if(cosphase>1.0L) {
	    cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
	    sunelong=0.0L;
	  } else if(cosphase<-1.0L) {
	    cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
	    sunelong=M_PI;
	  } else sunelong = acos(cosphase);
	  // Write a bunch of data to the output file.
	  outstream1 << "ieo" << idnumstring << goodsimct << " " << " ";
	  outstream1 << fixed << setprecision(6) << a/AU_KM << " " << e << " " << incl << " " << q/AU_KM << " " << Q/AU_KM << " " << mjdref << " " << n*SOLARDAY*DEGPRAD << " " << absmag << " ";
	  outstream1 << fixed << setprecision(6) << sundist/AU_KM << " " << obsdist/AU_KM << " " << sunelong*DEGPRAD << " " << phaseang*DEGPRAD << " ";
	  outstream1 << fixed << setprecision(6) << mjdobs[timect] << " " << outRA + astromsigma*gaussian_deviate_mt(generator)/3600.0L/cos(outDec/DEGPRAD) << " " << outDec + astromsigma*gaussian_deviate_mt(generator)/3600.0L << " ";
	  // Note that we reverse the sign here because the previously calculated
	  // quantities are the negatives of the true state vectors: that is, their
	  // vector origin is at the object, pointing toward the sun. The negatives
	  // we apply here solve this problem and produce the correct sun-to-object state vectors.
	  outstream1 << fixed << setprecision(3) << obsmag << " " << alt << " "; 
	  outstream1 << -targ_to_sun.x << " " << -targ_to_sun.y << " " << -targ_to_sun.z << " ";
	  outstream1 << -targvel_to_sunvel.x << " " << -targvel_to_sunvel.y << " " << -targvel_to_sunvel.z << " ";
	  outstream1 << fixed << setprecision(6) << outRA << " " << outDec << "\n";
	  // Close loop over all images
	}

	goodsimct++;
	cout << "Wrote successfully simulated object number " << goodsimct << " to output file\n";
	simct++;
	if(goodsimct>=compfac) {
	  compfac *= 10;
	  if(idnl>0) idnl--;
	  idnumstring={};
	  for(i=0; i<idnl; i++) idnumstring.push_back('0');
	}
      }// else cout << "Simulated object turned out not to be detectable\n";
      // Close loop on count of simulated objects
    }
  }
  outstream1.close();
  
  return(0);
}
