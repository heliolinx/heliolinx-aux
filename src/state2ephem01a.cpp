// September 29, 2024: state2ephem01a.cpp: given input in terms of
// barycentric state vectors, calculate an ephemeris.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define DEBUG 0
#define PHASE_G 0.15

// Note: configfile contains the masses and ephemerides for all of the planets.
// configfile also contains the name of an observatory code file. The input
// units of the state vector are AU for positions and km/sec for velocities.
// Internally, the program will use AU for positions and AU/timescale for
// velocities. The default value of timescale is one year, but it can be
// set to other values as needed. The input MJD file is just a single column
// of times in Modified Juliand Days.

static void show_usage()
{
  cerr << "Usage: state2ephem01a -cfg configfile -mjdfile mjdfile -mjdstart mjdstart -statevec x y z vx vy vz -obslon obslon -plxcos plxcos -plxsin plxsin -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string mjdfile;
  string planetfile;
  string outfile;
  long double mjdstart = 0.0L;
  int configread=0;
  int polyorder=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
  long double ldval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> targpos;
  vector <point3LD> targvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  point3LD targ_to_sun = point3LD(0,0,0);
  point3LD targvel_to_sunvel = point3LD(0,0,0);
  point3LD targ_to_obs = point3LD(0,0,0);
  point3LD obs_to_sun = point3LD(0,0,0);
  point3LD Sunposnow = point3LD(0,0,0);
  point3LD Sunvelnow = point3LD(0,0,0);
  long double absmag = 0.0l;
  long double cosphase=0L;
  long double phaseang=0L;
  long double sunelong=0L;
  long double phi1,phi2;
  long double obsdist=1.0L;
  long double sundist=1.0L;
  long double obsfromsun=1.0L;
  long double phaseslope = PHASE_G;
  long double light_travel_time;
  double obsmag;
  vector <point3LD> observerpos;
  vector <long double> obsMJD;
  int obsnum=0;
  int obsct=0;
  string lnfromfile;
  vector <string> linestringvec;
  int reachedeof=0;
  //  double obslon = 289.26345L;
  //  double plxcos = 0.865020L;
  //  double plxsin = -0.500901L;
  double obslon = 289.25058l; // Changed to match X05 obscode
  double plxcos = 0.864981l;
  double plxsin = -0.500958l;
  point3LD outpos = point3LD(0,0,0);
  long double outRA,outDec;
  ofstream outstream1;
  

  if(argc<16) {
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
    if(string(argv[i]) == "-mjdfile") {
      if(i+1 < argc) {
	// There is still something to read;
	mjdfile=argv[++i];
	i++;
      } else {
	cerr << "MJD file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-state" || string(argv[i]) == "-statevec") {
      if(i+6 < argc) {
	//There is still something to read;
	startpos.x = stold(argv[++i]);
	startpos.y = stold(argv[++i]);
	startpos.z = stold(argv[++i]);
	startvel.x = stold(argv[++i]);
	startvel.y = stold(argv[++i]);
	startvel.z = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "state vector keyword supplied with insufficient corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obslon") {
      if(i+1 < argc) {
	//There is still something to read;
	obslon=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Observatory longitude keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-plxcos") {
      if(i+1 < argc) {
	//There is still something to read;
	plxcos=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Observatory plxcos keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-plxsin") {
      if(i+1 < argc) {
	//There is still something to read;
	plxsin=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Observatory plxsin keyword supplied with too few corresponding arguments\n";
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
  cout << "input MJD file " << mjdfile << "\n";
  cout << "input starting MJD " << mjdstart << "\n";
  cout << "Observatory parameters: " << obslon << " " << plxcos << " " << plxsin << "\n";
  cout << "Input state vectors: " << startpos.x << " " << startpos.y << " " << startpos.z << " " << startvel.x << " " << startvel.y << " " << startvel.z << "\n";
  cout << "output file " << outfile << "\n";


  // Read input MJD file.
  instream1.open(mjdfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << mjdfile << "\n";
    return(1);
  }
  reachedeof=0;
  while(reachedeof==0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
    instream1 >> ldval;
    obsMJD.push_back(ldval);
    if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
  }
  obsnum=obsMJD.size();
  if(reachedeof==1) { 
    cout << obsnum << "lines read from file " << mjdfile << "\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();
  
  // Calculate the exact position of the observer at the time of each observation.
  observerpos={};
  for(obsct=0;obsct<obsnum;obsct++) {
    observer_barycoords01LD(obsMJD[obsct], 5, obslon, plxcos, plxsin, planetmjd, Earthpos, outpos);
    observerpos.push_back(outpos);
  }

  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "sundist obsdist sunelong phaseang obsMJD RA Dec mag Ast-Sun(J2000x)(km) Ast-Sun(J2000y)(km) Ast-Sun(J2000z)(km) Ast-Sun(J2000vx)(km/s) Ast-Sun(J2000vy)(km/s) Ast-Sun(J2000vz)(km/s) Ast-Bary(J2000x)(km) Ast-Bary(J2000y)(km) Ast-Bary(J2000z)(km) Ast-Bary(J2000vx)(km/s) Ast-Bary(J2000vy)(km/s) Ast-Bary(J2000vz)(km/s)\n";

  // Perform n-body integration beginning with the state vectors.
  integrate_orbit03LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, obsMJD, startpos, startvel, mjdstart, targpos, targvel);
  for(obsct=0;obsct<obsnum;obsct++) {
    planetposvel01LD(obsMJD[obsct],polyorder,planetmjd,Sunpos,Sunvel,Sunposnow,Sunvelnow);
 
    targ_to_obs.x = observerpos[obsct].x - targpos[obsct].x;
    targ_to_obs.y = observerpos[obsct].y - targpos[obsct].y;
    targ_to_obs.z = observerpos[obsct].z - targpos[obsct].z;
    targ_to_sun.x = Sunposnow.x - targpos[obsct].x;
    targ_to_sun.y = Sunposnow.y - targpos[obsct].y;
    targ_to_sun.z = Sunposnow.z - targpos[obsct].z;
    targvel_to_sunvel.x = Sunvelnow.x - targvel[obsct].x;
    targvel_to_sunvel.y = Sunvelnow.y - targvel[obsct].y;
    targvel_to_sunvel.z = Sunvelnow.z - targvel[obsct].z;

    obsdist = sqrt(targ_to_obs.x*targ_to_obs.x + targ_to_obs.y*targ_to_obs.y + targ_to_obs.z*targ_to_obs.z);
    sundist = sqrt(targ_to_sun.x*targ_to_sun.x + targ_to_sun.y*targ_to_sun.y + targ_to_sun.z*targ_to_sun.z);
    
    // Calculate observed magnitude not accounting for phase
    obsmag = absmag + 2.5*log10(obsdist*obsdist*sundist*sundist/AU_KM/AU_KM/AU_KM/AU_KM);	  
    // Calculate phase angle.
    cosphase = dotprod3LD(targ_to_obs,targ_to_sun)/obsdist/sundist;
    if(cosphase>1.0L) {
      cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
      phaseang=0.0L;
    } else if(cosphase<-1.0L) {
      cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
      phaseang=M_PI;
    } else phaseang = acos(cosphase);
    // Correct the magnitude calculated above for phase effects
    phi1 = exp(-PHASECONST_A1*pow(tan(phaseang/2.0),PHASECONST_B1));
    phi2 = exp(-PHASECONST_A2*pow(tan(phaseang/2.0),PHASECONST_B2));
    obsmag -= 2.5*log10((1.0-phaseslope)*phi1 + phaseslope*phi2);
    // Also calculate solar elongation, while we have the values handy
    obs_to_sun.x = Sunposnow.x - observerpos[obsct].x;
    obs_to_sun.y = Sunposnow.y - observerpos[obsct].y;
    obs_to_sun.z = Sunposnow.z - observerpos[obsct].z;
    targ_to_obs.x = targpos[obsct].x - observerpos[obsct].x; // Now it is really obs_to_targ
    targ_to_obs.y = targpos[obsct].y - observerpos[obsct].y;
    targ_to_obs.z = targpos[obsct].z - observerpos[obsct].z;
    obsfromsun = sqrt(obs_to_sun.x*obs_to_sun.x + obs_to_sun.y*obs_to_sun.y + obs_to_sun.z*obs_to_sun.z);
    cosphase = dotprod3LD(targ_to_obs,obs_to_sun)/obsfromsun/obsdist;
    if(cosphase>1.0L) {
      cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
      sunelong=0.0L;
    } else if(cosphase<-1.0L) {
      cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
      sunelong=M_PI;
    } else sunelong = acos(cosphase);
    // Initial approximation of the coordinates relative to the observer
    outpos.x = targpos[obsct].x - observerpos[obsct].x;
    outpos.y = targpos[obsct].y - observerpos[obsct].y;
    outpos.z = targpos[obsct].z - observerpos[obsct].z;
    light_travel_time = obsdist*1000.0/CLIGHT; // Factor of 1000 converts obsdist to meters
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x -= light_travel_time*targvel[obsct].x;
    outpos.y -= light_travel_time*targvel[obsct].y;
    outpos.z -= light_travel_time*targvel[obsct].z;
    // Light-travel-time corrected observer-target distance
    obsdist = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= obsdist;
    outpos.y /= obsdist;
    outpos.z /= obsdist;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, outRA, outDec);
    // Write a bunch of data to the output file.
    outstream1 << fixed << setprecision(7) << sundist/AU_KM << " " << obsdist/AU_KM << " ";
    outstream1 << fixed << setprecision(3) << sunelong*DEGPRAD << " " << phaseang*DEGPRAD << " ";
    outstream1 << fixed << setprecision(7) << obsMJD[obsct] << " " << outRA << " " << outDec << " ";
    // Note that we reverse the sign here because the previously calculated
    // quantities are the negatives of the true state vectors: that is, their
    // vector origin is at the object, pointing toward the sun. The negatives
    // we apply here solve this problem and produce the correct sun-to-object state vectors.
    outstream1 << fixed << setprecision(3) << obsmag << " "<< -targ_to_sun.x << " " << -targ_to_sun.y << " " << -targ_to_sun.z << " ";
    outstream1 << fixed << setprecision(6) << -targvel_to_sunvel.x << " " << -targvel_to_sunvel.y << " " << -targvel_to_sunvel.z << " ";
    outstream1 << fixed << setprecision(3) << targpos[obsct].x << " " << targpos[obsct].y << " " << targpos[obsct].z << " ";
    outstream1 << fixed << setprecision(6) << targvel[obsct].x << " " << targvel[obsct].y << " " << targvel[obsct].z << "\n";
  }
  outstream1.close();
  
  return(0);
}
