// February 14, 2025: arctrace01g.cpp
// Like arctrace01f.cpp, but takes the input observation file in
// heliolinc hldet format, and requires an observer code file.
//
// Description of ancestor program arctrace01f.cpp:
// Given an input file giving observer positions and measured RA and Dec
// for an object at various times, given as MJD, fit the object's orbit.
// Start by identifying a short sequence of observations that can be
// reasonably fit with a Keplerian orbit, and performing this fit using
// the Method of Herget. Then, using the approximate state vectors output
// from the Keplerian fit to perform a full, 6-dimensional orbit fit
// including planetary perturbations. This version tests the refactored
// arctrace02 function.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define MINCOLS 9
#define DEBUG 0

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain observations in heliolinc's hldet format.
static void show_usage()
{
  cerr << "Usage: arctrace01g -cfg configfile -observations obsfile -kepspan time_span_for_Keplerian_fit(day) -minchi min_chi_change -obscode obscodefile -outfile outfile -statefile statefile \n";
}

int main(int argc, char *argv[])
{
  ifstream instream1;
  int i, j, status, planetnum, planetct, pctSun, pctEarth, configread;
  int badread,reachedeof,obsnum,obsct;
  vector <hldet> detvec = {};
  vector <observatory> observatory_list = {};
  string obscodefile;
  i=j=status=configread=badread=reachedeof=obsnum=obsct=0;
  int polyorder = 3;
  string stest;
  string configfile;
  string obsfile;
  string planetfile;
  string outfile;
  string statefile;
  double kepspan=15.0;
  long double ldval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> planetpos_reverse;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> observer_heliopos;
  vector <long double> obsMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> sigastrom;
  string lnfromfile;
  vector <string> linestringvec;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  vector <long double> scalestate;
  make_LDvec(60,scalestate);
  vector <long double> newstate;
  make_LDvec(60,newstate);
  vector <long double> chivec;
  make_LDvec(70,chivec);
  vector <long double> beststatevec;
  long double bestchi=0L;
  vector <vector <long double>> simplex;
  make_LDmat(70,60,simplex);
  vector <double> bestRA;
  vector <double> bestDec;
  vector <double> bestresid;
  int refpoint = 0;
  long double minchichange = 0.001l;
  long double astromrms = 1.0;
  double obslon = 0.0l;
  double plxcos = 0.0l;
  double plxsin = 0.0l;
  point3LD obspos = point3LD(0,0,0);
  int verbose = 0;
  
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
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  cout.precision(17);  
  cout << "input configuration file " << configfile << "\n";
  cout << "input observation file " << obsfile << "\n";
  cout << "input timespan for Keplerian fit " << kepspan << "\n";
  cout << "general output file " << outfile << "\n";
  cout << "output state vector file " << statefile << "\n";

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
  observer_heliopos = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    status = obscode_lookup(observatory_list,detvec[obsct].obscode,obslon,plxcos,plxsin);
    if(status>0) {
      cerr << "ERROR: obscode_lookup failed for observatory code " << detvec[obsct].obscode << "\n";
      return(3);
    }
    // Calculate observer's exact barycentric position and velocity.
    observer_barycoords01LD(obsMJD[obsct], 5, (long double)obslon, (long double)plxcos, (long double)plxsin, planetmjd, Earthpos, obspos);
    // Calculate the sun's position at the same time.
    planetposvel01LD(obsMJD[obsct], polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
    // Convert obspos and obsvel from barycentric to heliocentric coords.
    outpos.x = obspos.x - outpos.x;
    outpos.y = obspos.y - outpos.y;
    outpos.z = obspos.z - outpos.z;
    observer_heliopos.push_back(outpos);
  }
  cout << obsMJD.size() << " " << obsMJD[obsMJD.size()-1] << " " << obsRA[obsMJD.size()-1] << " " << obsDec[obsMJD.size()-1] << "\n";
  if(obsnum<2 || long(obsRA.size())!=obsnum || long(obsDec.size())!=obsnum || long(observer_heliopos.size())!=obsnum) {
    cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << obsnum << " " << obsRA.size() << " " << obsDec.size() << " " << observer_heliopos.size() << "\n";
    return(1);
  }

  arctrace02(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_heliopos,obsMJD,obsRA,obsDec,sigastrom,kepspan,minchichange, bestRA, bestDec, bestresid, outpos, outvel, &bestchi, &astromrms, &refpoint);
  
  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromrms << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  cout << "State vectors correspond to reference point " << refpoint << " in the input planet files\n";
  
  ofstream outstream1 {outfile};
  for(obsct=0;obsct<obsnum;obsct++) {
    outstream1 << fixed << setprecision(8) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << bestRA[obsct] << " " << bestDec[obsct] << fixed << setprecision(8) << " " << (obsRA[obsct]-bestRA[obsct])*cos(obsDec[obsct]/DEGPRAD)*3600.0l << " " << (obsDec[obsct]-bestDec[obsct])*3600.0l << " " << bestresid[obsct] << " " << detvec[obsct].obscode << " " << detvec[obsct].mag << " " << detvec[obsct].band << "\n";
  }
  outstream1.close();

  cout << "Writing state vector filed called " << statefile << "\n";
  outstream1.open(statefile);
  outstream1 << fixed << setprecision(10) << planetmjd[refpoint] << " " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << " " << fixed << setprecision(5) << bestchi << " " << astromrms << "\n";
  outstream1.close();
  
  return(0);
}

