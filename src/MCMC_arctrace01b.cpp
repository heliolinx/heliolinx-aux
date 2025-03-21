// February 24, 2025: MCMC_arctrace01b.cpp: Use an MCMC to fit
// an orbit to input data in hldet format. Output RA, Dec, and
// angular velocities rather than state vectors.
//
// Uses fitting based on arctrace01_6Dh.cpp, which is
// like arctrace01_6D.cpp, but takes its input file in heliolinc hldet
// format, and requires an observatory code file.
//
// Description of ancestor program arctrace01_6D.cpp
// Given an input file giving observer positions and measured RA and Dec
// for an object at various times, given as MJD, fit the object's orbit
// based on a starting approximate state vector, and the MJD at which this
// matched the object's position. The input MJD must fall exactly on
// on step of the planet files read from the configuration file.
// Perform a full, 6-dimensional orbit fit
// including planetary perturbations. 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define MINCOLS 9
#define DEBUG 0
#define RESETNUM 2000
#define MINACCEPT 0.2
#define MAXACCEPT 0.5
#define POSSTEP 1000.0
#define VELSTEP 0.001
#define PRINTNUM 1000

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain
// MJD, observer X, Y, Z, VX, VY, VZ, RA, and Dec. The observer coordinates and velocities must
// be in km and km/sec, relative to the Sun. The RA and Dec must be in decimal degrees.
static void show_usage()
{
  cerr << "Usage: MCMC_arctrace01a -cfg configfile -ranseed random_number_seed -observations obsfile -mjd starting_mjd -statevec x y z vx vy vz -maxnum maxnum -obscode obscodefile -mjdend mjd_end -outfile outfile \n";
}

int main(int argc, char *argv[])
{
  ifstream instream1;
  int status, planetnum, planetct, pctSun, pctEarth, configread;
  long maxnum,itct,resetct,i, j;
  maxnum = 10000;
  itct = resetct = 0;
  int badread,reachedeof,obsnum,obsct;
  i=j=status=configread=badread=reachedeof=obsnum=obsct=0;
  int polyorder = 3;
  vector <hldet> detvec = {};
  vector <observatory> observatory_list = {};
  string obscodefile;
  double obslon = 0.0l;
  double plxcos = 0.0l;
  double plxsin = 0.0l;
  string stest;
  string configfile;
  string obsfile;
  string planetfile;
  string outfile;
  double kepspan=15.0;
  long double ldval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> observer_barypos;
  vector <long double> obsMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> sigastrom;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> resid;
  string lnfromfile;
  vector <string> linestringvec;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  point3LD obspos = point3LD(0,0,0);
  point3LD obsvel = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  point3LD testpos = point3LD(0,0,0);
  point3LD testvel = point3LD(0,0,0);
  point3LD lttpos = point3LD(0,0,0);
  vector <long double> beststatevec;
  double chisq = 0.0;
  double newchi = 0.0;
  double bestchi = LARGERR;
  vector <double> bestRA;
  vector <double> bestDec;
  vector <double> bestresid;
  long double astromrms = 1.0;
  long double mjdstart = 0.0l;
  long double mjdend = 0.0l;
  int planetfile_refpoint=0;
  int planetfile_endpoint=0;
  int verbose = 0;
  int changepar = 0;
  string seedstring;
  vector <double> stepvec;
  vector <double> keepvec;
  vector <double> normvec;
  make_dvec(6, stepvec);
  make_dvec(6, keepvec);
  make_dvec(6, normvec);
  stepvec[0] = stepvec[1] = stepvec[2] = POSSTEP;
  stepvec[3] = stepvec[4] = stepvec[5] = VELSTEP;
  vector <point3LD> orbit05pos;
  vector <point3LD> orbit05vel;
  vector <long double> orbit05MJD;
  long double light_travel_time,RA1,RA2,Dec1,Dec2,RAvel,Decvel;
  double bestephRA,bestephDec,bestephRAvel,bestephDecvel;
  bestephRA = bestephDec = bestephRAvel = bestephDecvel = 0;
  double timediff,dist,pa;
  
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
	  if(planetct == pctEarth) {
	    Earthpos = temppos;
	    Earthvel = tempvel;
	  }
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
    } else if(string(argv[i]) == "-mjd") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdstart=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Starting MJD keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-statevec" || string(argv[i]) == "-sv" || string(argv[i]) == "-svec" || string(argv[i]) == "--statevector" || string(argv[i]) == "--statevec" || string(argv[i]) == "-statevector" || string(argv[i]) == "--sv") {
      if(i+1 < argc) {
	//There is still something to read;
	startpos.x=stold(argv[++i]);
	cout << "read sv entry as " << startpos.x << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startpos.y=stold(argv[++i]);
	cout << "read sv entry as " << startpos.y << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startpos.z=stold(argv[++i]);
	cout << "read sv entry as " << startpos.z << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
 	startvel.x=stold(argv[++i]);
	cout << "read sv entry as " << startvel.x << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startvel.y=stold(argv[++i]);
	cout << "read sv entry as " << startvel.y << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startvel.z=stold(argv[++i]);
	cout << "read sv entry as " << startvel.z << "\n";
	i++;
      }
      else {
	cerr << "State vector keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
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
    } else if(string(argv[i]) == "-mjdend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdend=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Ending MJD keyword supplied with no corresponding argument\n";
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
  cout << "input observation file " << obsfile << "\n";
  cout << "input timespan for Keplerian fit " << kepspan << "\n";
  cout << "output file " << outfile << "\n";

  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, verbose);
  if(status!=0) {
    cerr << "ERROR reading observatory code file " << obscodefile << "\n";
    return(1);
  }
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";
  if(DEBUG>=2) {
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
  obsMJD = {};
  obsRA = obsDec = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD.push_back(detvec[obsct].MJD);
    obsRA.push_back(detvec[obsct].RA);
    obsDec.push_back(detvec[obsct].Dec);
    sigastrom.push_back(detvec[obsct].sig_across);
  }
  // Calculate the observer's heliocentric position at the time of each observation.
  observer_barypos = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    status = obscode_lookup(observatory_list,detvec[obsct].obscode,obslon,plxcos,plxsin);
    if(status>0) {
      cerr << "ERROR: obscode_lookup failed for observatory code " << detvec[obsct].obscode << "\n";
      return(3);
    }
    // Calculate observer's exact barycentric position and velocity.
    observer_barycoords01LD(obsMJD[obsct], 5, (long double)obslon, (long double)plxcos, (long double)plxsin, planetmjd, Earthpos, obspos);
    observer_barypos.push_back(obspos);
  }
  cout << obsMJD.size() << " " << obsMJD[obsMJD.size()-1] << " " << obsRA[obsMJD.size()-1] << " " << obsDec[obsMJD.size()-1] << "\n";
  if(obsnum<2 || long(obsRA.size())!=obsnum || long(obsDec.size())!=obsnum || long(observer_barypos.size())!=obsnum) {
    cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << obsnum << " " << obsRA.size() << " " << obsDec.size() << " " << observer_barypos.size() << "\n";
    return(1);
  }

  planetfile_refpoint=-99;
  j=0;
  while(j<long(planetmjd.size()) && planetfile_refpoint==-99) {
    j++;
    if(fabs(planetmjd[j]-mjdstart) <= IMAGETIMETOL/SOLARDAY) planetfile_refpoint=j;
  }
  if(planetfile_refpoint>=0) {
    cout << "Input mjd " << mjdstart << " corresponds to planet file entry number " << planetfile_refpoint << "\n";
  } else {
    cerr << "ERROR: input mjd " << mjdstart << " does not match any entry in the planet file\n";
    return(1);
  }
  planetfile_endpoint = -99;
  while(j<long(planetmjd.size()) && planetmjd[j]<mjdend) j++;
  planetfile_endpoint = j; // This is the first point in planetmjd that is after the end of the observations
  planetfile_endpoint += polyorder+1; // Pad endpoint so interpolation works.

  chisq = tortoisechi03(polyorder,planetnum,planetmjd,planetmasses,planetpos,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,startpos,startvel,planetfile_refpoint,fitRA,fitDec,resid);
  cout << "Initial chi-square value is " << chisq << "\n";
  ofstream outstream1 {outfile};
  changepar=0;
  for(i=0;i<6;i++) normvec[i] = keepvec[i] = 0.0;
  resetct=0;
  for(itct=0;itct<maxnum;itct++) {
    testpos = startpos;
    testvel = startvel;
    if(changepar==0) {
      testpos.x += 2.0*stepvec[changepar]*unitvar(generator) - stepvec[changepar];
    } else if(changepar==1) {
      testpos.y += 2.0*stepvec[changepar]*unitvar(generator) - stepvec[changepar];
    } else if(changepar==2) {
      testpos.z += 2.0*stepvec[changepar]*unitvar(generator) - stepvec[changepar];
    } else if(changepar==3) {
      testvel.x += 2.0*stepvec[changepar]*unitvar(generator) - stepvec[changepar];
    } else if(changepar==4) {
      testvel.y += 2.0*stepvec[changepar]*unitvar(generator) - stepvec[changepar];
    } else if(changepar==5) {
      testvel.z += 2.0*stepvec[changepar]*unitvar(generator) - stepvec[changepar];
    } 
    newchi = tortoisechi03(polyorder,planetnum,planetmjd,planetmasses,planetpos,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,testpos,testvel,planetfile_refpoint,fitRA,fitDec,resid);
    if(itct%PRINTNUM==0) cout << fixed << setprecision(10) << planetmjd[planetfile_refpoint] << " " << fixed << setprecision(3) << testpos.x << " " << testpos.y << " " << testpos.z << " "  << fixed << setprecision(10) << testvel.x << " " << testvel.y << " " << testvel.z << " " << newchi << " ";
    if(newchi<=chisq) {
      // The fit has become better: accept the point
      chisq = newchi;
      startpos = testpos;
      startvel = testvel;
      keepvec[changepar]+=1.0;
      normvec[changepar]+=1.0;
      // Predict the position at mjdend
      // Perform new integration.
      integrate_orbit05LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, startpos, startvel, planetfile_refpoint, planetfile_refpoint, planetfile_endpoint, orbit05MJD, orbit05pos, orbit05vel);
      // Interpolate to get object position at mjdend
      planetpos01LD(mjdend, polyorder, orbit05MJD, orbit05pos, testpos);
      planetpos01LD(mjdend, polyorder, orbit05MJD, orbit05vel, testvel);
      // Interpolate to get Earth's position and velocity at mjdend
      planetposvel01LD(mjdend, polyorder, planetmjd, Earthpos, Earthvel, obspos, obsvel);
      // Initial approximation of the coordinates relative to the observer
      lttpos.x = testpos.x - obspos.x;
      lttpos.y = testpos.y - obspos.y;
      lttpos.z = testpos.z - obspos.z;
      // Initial approximation of the observer-target distance
      ldval = vecabs3LD(lttpos);
      // Convert to meters and divide by the speed of light to get the light travel time.
      light_travel_time = ldval*1000.0/CLIGHT;
      // Light-travel-time corrected version of coordinates relative to the observer
      lttpos.x = testpos.x - light_travel_time*testvel.x - obspos.x;
      lttpos.y = testpos.y - light_travel_time*testvel.y - obspos.y;
      lttpos.z = testpos.z - light_travel_time*testvel.z - obspos.z;
       // Calculate unit vector
      vecnorm3LD(lttpos);
      // Project onto the celestial sphere.
      stateunitLD_to_celestial(lttpos, RA1, Dec1);
      // Calculate the position a little bit later
      lttpos.x = testpos.x - light_travel_time*testvel.x - obspos.x + TTDELTAT*(testvel.x - obsvel.x);
      lttpos.y = testpos.y - light_travel_time*testvel.y - obspos.y + TTDELTAT*(testvel.y - obsvel.y);
      lttpos.z = testpos.z - light_travel_time*testvel.z - obspos.z + TTDELTAT*(testvel.z - obsvel.z);
      // Calculate unit vector
      vecnorm3LD(lttpos);
      // Project onto the celestial sphere.
      stateunitLD_to_celestial(lttpos, RA2, Dec2);
      timediff = TTDELTAT/SOLARDAY; // TTDELTAT is in seconds, timediff is in days
      distradec02(RA1, Dec1, RA2, Dec2, &dist, &pa);
      RAvel = dist*sin(pa/DEGPRAD)/timediff; // Degrees per day
      Decvel = dist*cos(pa/DEGPRAD)/timediff; // Degrees per day
      // Write the point to the output file
      outstream1 << fixed << setprecision(6) << mjdend << " " << RA1 << " " << Dec1 << " " << RAvel << " "  << Decvel << " " << newchi << "\n";
      if(itct%PRINTNUM==0) cout << "ACCEPTED" << " itct " << itct << "\n";
      if(newchi<=bestchi) {
	bestchi = newchi;
	outpos = startpos;
	outvel = startvel;
	bestRA = fitRA;
	bestDec = fitDec;
	bestresid = resid;
	bestephRA = RA1;
	bestephDec = Dec1;
	bestephRAvel = RAvel;
	bestephDecvel = Decvel;
	ldval = 0.0;
	for(obsct=0;obsct<obsnum;obsct++) ldval += bestresid[obsct]*bestresid[obsct];
	astromrms = sqrt(ldval/double(obsnum));
      }
    } else {
      double testpar = unitvar(generator);
      if(testpar<exp(0.5*(chisq-newchi))) {
	// The fit has become only slightly worse: accept the point
	chisq = newchi;
	startpos = testpos;
	startvel = testvel;
	keepvec[changepar]+=1.0;
	normvec[changepar]+=1.0;
	// Predict the position at mjdend
	// Perform new integration.
	integrate_orbit05LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, startpos, startvel, planetfile_refpoint, planetfile_refpoint, planetfile_endpoint, orbit05MJD, orbit05pos, orbit05vel);
	// Interpolate to get object position at mjdend
	planetpos01LD(mjdend, polyorder, orbit05MJD, orbit05pos, testpos);
	planetpos01LD(mjdend, polyorder, orbit05MJD, orbit05vel, testvel);
	// Interpolate to get Earth's position and velocity at mjdend
	planetposvel01LD(mjdend, polyorder, planetmjd, Earthpos, Earthvel, obspos, obsvel);
	// Initial approximation of the coordinates relative to the observer
	lttpos.x = testpos.x - obspos.x;
	lttpos.y = testpos.y - obspos.y;
	lttpos.z = testpos.z - obspos.z;
	// Initial approximation of the observer-target distance
	ldval = vecabs3LD(lttpos);
	// Convert to meters and divide by the speed of light to get the light travel time.
	light_travel_time = ldval*1000.0/CLIGHT;
	// Light-travel-time corrected version of coordinates relative to the observer
	lttpos.x = testpos.x - light_travel_time*testvel.x - obspos.x;
	lttpos.y = testpos.y - light_travel_time*testvel.y - obspos.y;
	lttpos.z = testpos.z - light_travel_time*testvel.z - obspos.z;
	// Calculate unit vector
	vecnorm3LD(lttpos);
	// Project onto the celestial sphere.
	stateunitLD_to_celestial(lttpos, RA1, Dec1);
	// Calculate the position a little bit later
	lttpos.x = testpos.x - light_travel_time*testvel.x - obspos.x + TTDELTAT*(testvel.x - obsvel.x);
	lttpos.y = testpos.y - light_travel_time*testvel.y - obspos.y + TTDELTAT*(testvel.y - obsvel.y);
	lttpos.z = testpos.z - light_travel_time*testvel.z - obspos.z + TTDELTAT*(testvel.z - obsvel.z);
	// Calculate unit vector
	vecnorm3LD(lttpos);
	// Project onto the celestial sphere.
	stateunitLD_to_celestial(lttpos, RA2, Dec2);
	timediff = TTDELTAT/SOLARDAY; // TTDELTAT is in seconds, timediff is in days
	distradec02(RA1, Dec1, RA2, Dec2, &dist, &pa);
	RAvel = dist*sin(pa/DEGPRAD)/timediff; // Degrees per day
	Decvel = dist*cos(pa/DEGPRAD)/timediff; // Degrees per day
	// Write the point to the output file
	outstream1 << fixed << setprecision(6) << mjdend << " " << RA1 << " " << Dec1 << " " << RAvel << " "  << Decvel << " " << newchi << "\n";
	if(itct%PRINTNUM==0) cout << "ACCEPTED" << " itct " << itct << "\n";
      } else {
	// Reject the point
	normvec[changepar]+=1.0;
	if(itct%PRINTNUM==0) cout << "REJECTED" << " itct " << itct << "\n";
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
    // Cycle changepar	  
    if(changepar<5 && changepar>=0) changepar++;
    else changepar=0;
  }
  
  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromrms << "\n";
  cout << "Best ephemeris prection at MJD " << fixed << setprecision(6) << mjdend << " was position " << bestephRA << " " << bestephDec << " and angular velocity " << bestephRAvel << " "  << bestephDecvel << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  cout << "State vectors correspond to reference point " << planetfile_refpoint << " in the input planet files\n";

  outstream1.close();
  
  return(0);
}

