// March 13, 2025: tracklet_arctrace01c.cpp:
// Like tracklet_arctrace01b.cpp, but interpolates the input
// observation file in time, in order to track objects with high
// angular velocities more accurately and efficiently.

// Description of ancestor program tracklet_arctrace01b.cpp
// Like tracklet_arctrace01a.cpp, but takes the input observation file in
// heliolinc hldet format, and requires an observer code file.
// Also writes an output file in hldet format.
//
// Description of ancestor program tracklet_arctrace01a.cpp
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
#define TIMERANGEFAC 2000000000

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain
// MJD, observer X, Y, Z, VX, VY, VZ, RA, and Dec. The observer coordinates and velocities must
// be in km and km/sec, relative to the Sun. The RA and Dec must be in decimal degrees.
static void show_usage()
{
  cerr << "Usage: tracklet_arctrace01b -cfg configfile -observations obsfile  -obscode obscodefile -mjd mjd_for_statevecs -statevec x y z vx vy vz -minchi min_chi_change -mjdstart mjdstart -mjdend mjdend -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file  -timetol MJD_matching_tolerance(days) -skytol sky_matching_radius(deg) -veltol velocity_matching_radius(deg/day) -max_astrom_rms max astrometric RMS (arcsec) -outfile outfile -logfile logfile\n";
}

int main(int argc, char *argv[])
{
  ifstream instream1;
  long i, j, k;
  int status, planetnum, planetct, pctSun, pctEarth, configread;
  long badread,reachedeof,obsnum,obsct,tct;
  i=j=k=status=configread=badread=reachedeof=obsnum=obsct=tct=0;
  int polyorder = 3;
  vector <observatory> observatory_list = {};
  string obscodefile;
  vector <hldet> obsdetvec = {};
  vector <hldet> outdetvec = {};
  vector <hldet> outdetvec2 = {};
  string stest;
  string configfile;
  string obsfile;
  string planetfile;
  string logfile;
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
  vector <point3LD> observer_heliopos;
  vector <point3LD> observer_barypos;
  vector <long double> obsMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> sigastrom;
  vector <point3LD> observer_barypos2;
  vector <long double> obsMJD2;
  vector <double> obsRA2;
  vector <double> obsDec2;
  vector <double> sigastrom2;
  string lnfromfile;
  vector <string> linestringvec;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
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
  long double minchichange = 0.001l;
  long double astromRMS = 1.0;
  long double mjd_statevecs = 0.0l;
  long planetfile_refpoint=0;
  point3LD tppos = point3LD(0,0,0);
  long double mjdstart = 0;
  long double mjdend = 0;
  long planetfile_startpoint,planetfile_endpoint;
  long bestmatch=0;
  double bestRMS=LARGERR;
  vector <point3LD> orbit05pos;
  vector <point3LD> orbit05vel;
  vector <long double> orbit05MJD;
  long double light_travel_time = 0.0l;
  long double timediff = 0.0l;

  string pairdetfile="pairdetfile01.csv";
  string trackletfile="trackletfile01.csv";
  string trk2detfile="trk2detfile01.csv";
  string imfile;
  vector <hldet> detvec = {};
  vector <hlimage> image_log;
  vector <tracklet> tracklets;
  vector <longpair> trk2det;
  long imnum,pairnum,pairct,i1,i2,detnum;
  double dist,pa;

  vector <double> trkMJD;
  vector <double> trkRA;
  vector <double> trkDec;
  vector <double> trkRAvel;
  vector <double> trkDecvel;
  vector <double> ephMJD;
  vector <double> ephRA;
  vector <double> ephDec;
  vector <double> ephRAvel;
  vector <double> ephDecvel;
  double RAvel,Decvel;
  double delem = 0.0;
  long ind = 0;
  double_index one_dind = double_index(delem,ind);
  vector <double_index> trk_dindvec;

  vector <long> matching_trkind;
  vector <long> matching_trkind2;
  vector <long> trkvec;
  vector <long> trkvec_temp;
  long tracklets_added = 0;
  vector <long> detections_added;
  double timetol=1.0;
  double skytol=1.0;
  double veltol=0.1;
  double vdist=0.0;
  ofstream outstream1;
  int verbose=0;
  double obslon = 0.0l;
  double plxcos = 0.0l;
  double plxsin = 0.0l;
  point3LD obspos = point3LD(0,0,0);
  point3LD targpos = point3LD(0,0,0);
  point3LD obsvel = point3LD(0,0,0);
  point3LD targvel = point3LD(0,0,0);
  long ephct=0;
  double max_astrom_rms = 1.0;
  long queryrange;
  point6ix2 onepoint = point6ix2(0,0,0,0,0,0,0,0);
  point6ix2 querypoint = point6ix2(0,0,0,0,0,0,0,0);
  vector <point6ix2> alltrackvec;
  double earliest_mjd,latest_mjd,central_mjd;
  double timeintscale,posintscale,velintscale;
  point3d projpt = point3d(0.0,0.0,0.0);
  long unsigned int splitpoint = 0;
  long unsigned int kdroot = 0;
  KD_point6ix2 kdpoint = KD_point6ix2(onepoint,-1,-1,1,-1);
  vector <KD_point6ix2> kdtree;
  vector <long> queryout;
  
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
	mjd_statevecs=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Keplerian time-span keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-mjdstart") {
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
    } else if(string(argv[i]) == "-mjdend") {
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
    } else if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
      if(i+1 < argc) {
	//There is still something to read;
        pairdetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output paired detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	//There is still something to read;
	imfile=argv[++i];
	i++;
      }
      else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-tf" || string(argv[i]) == "-trkfile" || string(argv[i]) == "-trackletfile" || string(argv[i]) == "--trackletfile") {
      if(i+1 < argc) {
	//There is still something to read;
	trackletfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
      if(i+1 < argc) {
	//There is still something to read;
	trk2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timetol") {
      if(i+1 < argc) {
	//There is still something to read;
	timetol=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Time tolerance keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-skytol") {
      if(i+1 < argc) {
	//There is still something to read;
	skytol=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Sky tolerance keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-veltol") {
      if(i+1 < argc) {
	//There is still something to read;
	veltol=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Velocity tolerance keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-max_astrom_rms" || string(argv[i]) == "-maxastromrms" || string(argv[i]) == "-astrom_rms" || string(argv[i]) == "-arms" || string(argv[i]) == "-marms" || string(argv[i]) == "-max_astrometric_rms" || string(argv[i]) == "--max_arms" ) {
      if(i+1 < argc) {
	//There is still something to read;
	max_astrom_rms=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Max astrometric RMS keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-log" || string(argv[i]) == "-logfile" || string(argv[i]) == "--logfile") {
      if(i+1 < argc) {
	//There is still something to read;
	logfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output log file keyword supplied with no corresponding argument\n";
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
	cerr << "Output hldet file keyword supplied with no corresponding argument\n";
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
  cout << "output file " << logfile << "\n";

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
  obsdetvec={};
  status=read_hldet_file(obsfile, obsdetvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << obsfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  obsnum = obsdetvec.size();
  outdetvec = obsdetvec;
  cout << "Read " << obsnum << " data lines from observation file " << obsfile << "\n";

  // Load the MJD, RA, and Dec vectors needed by arctrace02.
  obsMJD = {};
  obsRA = obsDec = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD.push_back(obsdetvec[obsct].MJD);
    obsRA.push_back(obsdetvec[obsct].RA);
    obsDec.push_back(obsdetvec[obsct].Dec);
    sigastrom.push_back(obsdetvec[obsct].sig_across);
  }
  // Calculate the observer's heliocentric position at the time of each observation.
  observer_heliopos = observer_barypos = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    status = obscode_lookup(observatory_list,obsdetvec[obsct].obscode,obslon,plxcos,plxsin);
    if(status>0) {
      cerr << "ERROR: obscode_lookup failed for observatory code " << obsdetvec[obsct].obscode << "\n";
      return(3);
    }
    // Calculate observer's exact barycentric position and velocity.
    observer_barycoords01LD(obsMJD[obsct], 5, (long double)obslon, (long double)plxcos, (long double)plxsin, planetmjd, Earthpos, obspos);
    observer_barypos.push_back(obspos);
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

  planetfile_refpoint=-99;
  j=0;
  while(j<long(planetmjd.size()) && planetfile_refpoint==-99) {
    j++;
    if(fabs(planetmjd[j]-mjd_statevecs) <= IMAGETIMETOL/SOLARDAY) planetfile_refpoint=j;
  }
  if(planetfile_refpoint>=0) {
    cout << "Input mjd " << mjd_statevecs << " corresponds to planet file entry number " << planetfile_refpoint << "\n";
  } else {
    cerr << "ERROR: input mjd " << mjd_statevecs << " does not match any entry in the planet file\n";
    return(1);
  }
  
  planetfile_startpoint = planetfile_endpoint = -99;
  j=0;
  while(j<long(planetmjd.size()) && planetmjd[j]<=mjdstart) j++;
  planetfile_startpoint = j-1; // This is the last point in planetmjd that is before mjdstart.
  j=planetfile_startpoint;
  while(j<long(planetmjd.size()) && planetmjd[j]<=mjdend) j++;
  planetfile_endpoint = j; // This is the first point in planetmjd that is after mjdend.
  cout << "MJD range " << mjdstart << " to " << mjdend << " for ephemeris corresponds to planet file index range from " << planetfile_startpoint << " to " << planetfile_endpoint << "\n";
  
  arc6D01(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,minchichange,planetfile_refpoint,startpos,startvel,bestRA,bestDec,bestresid,outpos, outvel, &bestchi, &astromRMS);
  
  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromRMS << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  cout << "State vectors correspond to reference point " << planetfile_refpoint << " in the input planet files\n";

  // FINISHED WITH ORBIT FIT TO ORIGINAL INPUT DATA

  // CALCULATE EPHEMERIS BASED ON ORBIT FIT
  planetfile_startpoint = planetfile_endpoint = -99;
  j=0;
  while(j<long(planetmjd.size()) && planetmjd[j]<=mjdstart) j++;
  planetfile_startpoint = j-1; // This is the last point in planetmjd that is before mjdstart.
  j=planetfile_startpoint;
  while(j<long(planetmjd.size()) && planetmjd[j]<=mjdend) j++;
  planetfile_endpoint = j; // This is the first point in planetmjd that is after mjdend.
  planetfile_startpoint -= polyorder+1;
  planetfile_endpoint += polyorder+1;
  if(planetfile_startpoint<0) planetfile_startpoint=0;
  if(planetfile_endpoint>long(planetmjd.size()-1)) planetfile_endpoint=planetmjd.size()-1;
  
  startpos = outpos;
  startvel = outvel;
  integrate_orbit05LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, startpos, startvel, planetfile_startpoint, planetfile_refpoint, planetfile_endpoint, orbit05MJD, orbit05pos, orbit05vel);

  // Load ephMJD vector
  ephMJD={};
  ephMJD.push_back(planetmjd[planetfile_startpoint+polyorder+1]);
  long tmax = (planetmjd[planetfile_endpoint-polyorder-1] - planetmjd[planetfile_startpoint+polyorder+1])/timetol;
  if(tmax<0) {
    cerr << "ERROR: starting time and ending time are too close together\n";
    return(7);
  }
  tmax+=2;
  for(long tct=1;tct<=tmax;tct++) {
    ldval = planetmjd[planetfile_startpoint+polyorder+1] + (long double)tct*(planetmjd[planetfile_endpoint-polyorder-1] - planetmjd[planetfile_startpoint+polyorder+1])/((long double)tmax);
    if(ldval>planetmjd[planetfile_endpoint-polyorder-1]) ldval=planetmjd[planetfile_endpoint-polyorder-1];
    ephMJD.push_back(ldval);
  }
  
  for(ephct=0;ephct<long(ephMJD.size());ephct++) {
    // Interpolate to get coordinates of target and observer
    planetpos01LD(ephMJD[ephct], polyorder, planetmjd, Earthpos, obspos);
    planetpos01LD(ephMJD[ephct], polyorder, planetmjd, Earthvel, obsvel);
    planetpos01LD(ephMJD[ephct], polyorder, orbit05MJD, orbit05pos, targpos);
    planetpos01LD(ephMJD[ephct], polyorder, orbit05MJD, orbit05vel, targvel);

    // Initial approximation of the coordinates relative to the observer
    outpos.x = targpos.x - obspos.x;
    outpos.y = targpos.y - obspos.y;
    outpos.z = targpos.z - obspos.z;
    // Initial approximation of the observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = targpos.x - light_travel_time*targvel.x - obspos.x;
    outpos.y = targpos.y - light_travel_time*targvel.y - obspos.y;
    outpos.z = targpos.z - light_travel_time*targvel.z - obspos.z;
    // Light-travel-time corrected observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    long double RA1,RA2,Dec1,Dec2;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, RA1, Dec1);
    // Calculate the position a little bit later
    outpos.x = targpos.x - light_travel_time*targvel.x - obspos.x;
    outpos.y = targpos.y - light_travel_time*targvel.y - obspos.y;
    outpos.z = targpos.z - light_travel_time*targvel.z - obspos.z;
    outpos.x += TTDELTAT*(targvel.x - obsvel.x);
    outpos.y += TTDELTAT*(targvel.y - obsvel.y);
    outpos.z += TTDELTAT*(targvel.z - obsvel.z);
    //cout << "Earthvel: " << obsvel.x << " " << obsvel.y << " " << obsvel.z << "\n";
    //cout << "targvel: " << targvel.x << " " << targvel.y << " " << targvel.z << "\n";
    //cout << "diffs: " << (targvel.x - obsvel.x) << " " << (targvel.y - obsvel.y) << " " << (targvel.z - obsvel.z) << "\n";
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, RA2, Dec2);
    timediff = TTDELTAT/SOLARDAY; // TTDELTAT is in seconds, timediff is in days
    distradec02(RA1, Dec1, RA2, Dec2, &dist, &pa);
    RAvel = dist*sin(pa/DEGPRAD)/timediff; // Degrees per day
    Decvel = dist*cos(pa/DEGPRAD)/timediff; // Degrees per day
    ephRA.push_back(RA1);
    ephDec.push_back(Dec1);
    ephRAvel.push_back(RAvel);
    ephDecvel.push_back(Decvel);
  }
  // FINISHED WITH EPHEMERIS CALCULATION.
  outstream1.open(logfile);
  for(ephct=0;ephct<long(ephMJD.size());ephct++) {
    cout << fixed << setprecision(6) << "ephemeris point " << ephct << " " << ephMJD[ephct] << " " << ephRA[ephct] << " " << ephDec[ephct] << " " << ephRAvel[ephct] << " " << ephDecvel[ephct] << "\n";
    outstream1 << fixed << setprecision(6) << "ephemeris point " << ephct << " " << ephMJD[ephct] << " " << ephRA[ephct] << " " << ephDec[ephct] << " " << ephRAvel[ephct] << " " << ephDecvel[ephct] << "\n";
  }
  earliest_mjd = ephMJD[0];
  latest_mjd = ephMJD[ephMJD.size()-1];
  // READ TRACKLET FILES, AND LOAD TRACKLET VECTORS
  detvec={};
  status=read_pairdet_file(pairdetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  detnum = detvec.size();
  cout << "Read " << detnum << " data lines from paired detection file " << pairdetfile << "\n";
  
  image_log={};
  status=read_image_file2(imfile, image_log);
  if(status!=0) {
    cerr << "ERROR: could not successfully read image file " << imfile << "\n";
    cerr << "read_image_file2 returned status = " << status << ".\n";
   return(1);
  }
  imnum = image_log.size();
  cout << "Read " << imnum << " data lines from image file " << imfile << "\n";
  
  tracklets={};
  status=read_tracklet_file(trackletfile, tracklets, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet file " << trackletfile << "\n";
    cerr << "read_tracklet_file returned status = " << status << ".\n";
   return(1);
  }
  pairnum = tracklets.size();
  cout << "Read " << pairnum << " data lines from tracklet file " << trackletfile << "\n";
  
  trk2det={};
  status=read_longpair_file(trk2detfile, trk2det, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read trk2det file " << trk2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << trk2det.size() << " data lines from trk2det file " << trk2detfile << "\n";

  // Main loop over all tracklets, to load tracket vectors.
  trkMJD = trkRA = trkDec = trkRAvel = trkDecvel = {};
  trk_dindvec = {};
  for(pairct=0; pairct<pairnum; pairct++) {
    // Obtain indices to the image_log and heliocentric distance vectors.
    i1=tracklets[pairct].Img1;
    i2=tracklets[pairct].Img2;
    if(i1<0 || i1>=imnum || i2<0 || i2>=imnum) {
      cerr << "ERROR: image indices " << i1 << ", " << i2 << " not in range 0 to " << imnum << " spanned by image log\n";
      return(2);
    }
    // MJD, RA, and Dec for the tracklet will refer to
    // its first point.
    trkMJD.push_back(image_log[i1].MJD);
    trkRA.push_back(tracklets[pairct].RA1);
    trkDec.push_back(tracklets[pairct].Dec1);
    // Calculate time difference between the observations
    timediff = (image_log[i2].MJD - image_log[i1].MJD);
    distradec02(tracklets[pairct].RA1, tracklets[pairct].Dec1, tracklets[pairct].RA2, tracklets[pairct].Dec2, &dist, &pa);
    RAvel = dist*sin(pa/DEGPRAD)/timediff; // Degrees per day
    Decvel = dist*cos(pa/DEGPRAD)/timediff; // Degrees per day
    trkRAvel.push_back(RAvel);
    trkDecvel.push_back(Decvel);
    one_dind = double_index(image_log[i1].MJD,pairct);
    trk_dindvec.push_back(one_dind);
  }
  sort(trk_dindvec.begin(), trk_dindvec.end(), lower_double_index());
     
  // FINISHED LOADING TRACKLET VECTORS
  
  // LOAD K-D TREE FROM TRACKLET VECTORS
  
  // Establish integer scalings for k-d tree
  if(trk_dindvec[0].delem < earliest_mjd) earliest_mjd = trk_dindvec[0].delem;
  if(trk_dindvec[trk_dindvec.size()-1].delem > latest_mjd) latest_mjd = trk_dindvec[trk_dindvec.size()-1].delem;
  central_mjd = 0.5*earliest_mjd + 0.5*latest_mjd;
  timeintscale = (latest_mjd - earliest_mjd)/TIMERANGEFAC;
  posintscale = timeintscale*skytol/timetol;
  velintscale = timeintscale*veltol/timetol;
  queryrange = timetol/timeintscale+1.0;
  cout << "Integerized tolerance check: time, sky position, velocity: " << double(queryrange)*timeintscale << " " << double(queryrange)*posintscale << " " << double(queryrange)*velintscale << "\n";
  queryrange *= sqrt(3.0);
  
  // Load integerized vector alltrackvec
  alltrackvec = {};
  for(i=0;i<long(trk_dindvec.size());i++) {
    pairct = trk_dindvec[i].index;
    projpt = celeproj01(trkRA[pairct], trkDec[pairct]);
    onepoint = point6ix2(int((trkMJD[pairct] - central_mjd)/timeintscale), int(projpt.x*DEGPRAD/posintscale + 0.5), int(projpt.y*DEGPRAD/posintscale + 0.5), int(projpt.z*DEGPRAD/posintscale + 0.5), int(trkRAvel[pairct]/velintscale + 0.5), int(trkDecvel[pairct]/velintscale + 0.5),pairct,i);
    alltrackvec.push_back(onepoint);
  }
  // Form k-d tree
  kdtree={};
  kdroot = splitpoint = 0;
  splitpoint=medind_6ix2(alltrackvec,1);
  kdpoint = KD_point6ix2(alltrackvec[splitpoint],-1,-1,1,-1);
  kdtree.push_back(kdpoint);
  kdtree_6i01(alltrackvec,1,splitpoint,kdroot,kdtree);
  cout << "Created a k-d tree with " << kdtree.size() << " branches\n";

  // SEARCH FOR TRACKLETS THAT MATCH THE EPHEMERIS PREDICTIONS
  j=0;
  matching_trkind = {};
  for(i=0;i<long(ephMJD.size());i++) {
    projpt = celeproj01(ephRA[i], ephDec[i]);
    querypoint = point6ix2(int((ephMJD[i] - central_mjd)/timeintscale), int(projpt.x*DEGPRAD/posintscale + 0.5), int(projpt.y*DEGPRAD/posintscale + 0.5), int(projpt.z*DEGPRAD/posintscale + 0.5), int(ephRAvel[i]/velintscale + 0.5), int(ephDecvel[i]/velintscale + 0.5),pairct,i);
    queryout = {};
    kdrange_6i01(kdtree, querypoint, queryrange, queryout);
    // Loop over any matching tracklets
    for(k=0;k<long(queryout.size());k++) {
      pairct = kdtree[queryout[k]].point.i1; // Index to trkMJD, trkRA, trkDec, trkRAvel, and trkDecvel
      dist = distradec01(ephRA[i],ephDec[i],trkRA[pairct],trkDec[pairct]);
      vdist = sqrt(DSQUARE(ephRAvel[i]-trkRAvel[pairct])+DSQUARE(ephDecvel[i]-trkDecvel[pairct]));
      // See if points in this tracklet match any of the inputs points
      long matchinput=0;
      trkvec={};
      trkvec = tracklet_lookup(trk2det, pairct);
      if(trkvec.size()<=0) {
	cerr << "ERROR: tracklet lookup failed for  tracklet line " << pairct << ":\n";
	return(2);
      }
      for(obsct=0;obsct<long(obsMJD.size());obsct++) {
	for(tct=0;tct<long(trkvec.size());tct++) {
	  if(fabs(obsMJD[obsct]-image_log[detvec[trkvec[tct]].image].MJD)<IMAGETIMETOL) matchinput=1;
	}
      }
      if(matchinput==0) {
	// We matched a tracklet that was not included in the input linkage.
	// Report it to the user, and record it for inclusion in the orbit fit
	cout << fixed << setprecision(4) << "Match found on tracklet line " << pairct << ":" << trkRA[pairct] << " " << trkDec[pairct] << " " << trkRAvel[pairct] << " " << trkDecvel[pairct] << ":\n";
	cout << fixed << setprecision(4) << detvec[trkvec[0]].idstring << " dist = " << dist << " degrees, vdist = " << vdist << " deg/day :\n" << ephMJD[i] << " " << ephRA[i] << " " << ephDec[i] << " " << ephRAvel[i] << " " << ephDecvel[i] << "\n";
	outstream1 << fixed << setprecision(4) << "Match found on tracklet line " << pairct << ":" << trkRA[pairct] << " " << trkDec[pairct] << " " << trkRAvel[pairct] << " " << trkDecvel[pairct] << ":\n";
	outstream1 << fixed << setprecision(4) << detvec[trkvec[0]].idstring << " dist = " << dist << " degrees, vdist = " << vdist << " deg/day :\n" << ephMJD[i] << " " << ephRA[i] << " " << ephDec[i] << " " << ephRAvel[i] << " " << ephDecvel[i] << "\n";
	matching_trkind.push_back(pairct);
      }
    }
  }

  if(matching_trkind.size()<=0) {
    cout << "No matching tracklets were found\n";
    return(0);
  }
  // De-duplicate matching_trkind
  matching_trkind2 = matching_trkind;
  sort(matching_trkind2.begin(),matching_trkind2.end());
  matching_trkind = {};
  i=0;
  matching_trkind.push_back(matching_trkind2[i]);
  for(i=1;i<long(matching_trkind2.size());i++) {
    if(matching_trkind2[i]!=matching_trkind2[i-1]) {
      matching_trkind.push_back(matching_trkind2[i]);
    }
  }
  
  cout << matching_trkind.size() << " potentially matching tracklets were found\n";
  outstream1 << matching_trkind.size() << " potentially matching tracklets were found\n";
  long matchct=0;
  long matchnum=0;
  bestmatch = -1;
  bestRMS=LARGERR;
  matching_trkind2 = {};
  for(matchct=0;matchct<long(matching_trkind.size());matchct++) {
    pairct = matching_trkind[matchct];
    trkvec={};
    trkvec = tracklet_lookup(trk2det, pairct);
    cout << "Attempting orbit fit for potential match number " << matchct+1 << ", with " << trkvec.size() << " points\n";
    observer_barypos2 = {};
    obsMJD2 = {};
    obsRA2 = obsDec2 = sigastrom2 = {};
    outdetvec2 = {};
    obsct=tct=0;
      while(obsct<long(obsMJD.size()) || tct<long(trkvec.size())) {
        if(obsct<long(obsMJD.size()) && (tct>=long(trkvec.size()) || obsMJD[obsct] < image_log[detvec[trkvec[tct]].image].MJD)) {
	  // Next point in the time-ordered set is from the orginal observation vectors
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	} else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && (obsRA[obsct] != detvec[trkvec[tct]].RA || obsDec[obsct] != detvec[trkvec[tct]].Dec)) {
	  // We have two points at the same time, but they are not identical
	  // Keep the one from the original vectors
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	  // And also keep the one from the new tracklet
	  obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	  obsRA2.push_back(detvec[trkvec[tct]].RA);
	  obsDec2.push_back(detvec[trkvec[tct]].Dec);
	  sigastrom2.push_back(1.0);
	  // Convert input heliocentric observer position to barycentric.
	  planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	  tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	  observer_barypos2.push_back(tppos);
	  outdetvec2.push_back(detvec[trkvec[tct]]);
	  tct++;
	} else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && obsRA[obsct] == detvec[trkvec[tct]].RA && obsDec[obsct] == detvec[trkvec[tct]].Dec ) {
	  // We have two exactly identical points. Keep the one from the original vectors...
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	  // ...and discard the one from the new tracklet
	  tct++;
	} else if (tct<long(trkvec.size())) {
	  // Next point in the time-ordered set is from the new tracklet
	  obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	  obsRA2.push_back(detvec[trkvec[tct]].RA);
	  obsDec2.push_back(detvec[trkvec[tct]].Dec);
	  sigastrom2.push_back(1.0);
	  // Convert input heliocentric observer position to barycentric.
	  planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	  tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	  observer_barypos2.push_back(tppos);
	  outdetvec2.push_back(detvec[trkvec[tct]]);
	  tct++;
	} else {
	  cerr << "Error: logically excluded case in loading of new vectors\n";
	  return(5);
	}
      }
    cout << "Size checks: " << obsMJD.size() << " "  << obsRA.size() << " "  << obsDec.size() << " "  << sigastrom.size() << " "  << sigastrom.size() << "\n";
    // Perform orbit-fit to augmented input data
    cout << "Launching arc6D01 for match " << matchct+1  << " of " << matching_trkind.size() << " with " << obsMJD2.size() << " data points\n";
    outstream1 << "Launching arc6D01 for match " << matchct+1  << " = " << detvec[trkvec[0]].idstring << " of " << matching_trkind.size() << " at MJD " << detvec[trkvec[0]].MJD << ", with " << trkvec.size() << " tracklet points and hence " << obsMJD2.size() << " total data points\n";
    arc6D01(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_barypos2,obsMJD2,obsRA2,obsDec2,sigastrom2,minchichange,planetfile_refpoint,startpos,startvel,bestRA,bestDec,bestresid,outpos, outvel, &bestchi, &astromRMS);
    outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
    if(astromRMS<max_astrom_rms) {
      // The fit was good
      matchnum+=1;
      matching_trkind2.push_back(matching_trkind[matchct]);
      outstream1 << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
      if(astromRMS<bestRMS) {
	if(bestRMS>max_astrom_rms) {
	  outstream1 << "the first good fit to be identified\n";
	} else {
	  outstream1 << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
	}
	bestRMS = astromRMS;
	bestmatch = matchct;
      } else {
	outstream1 << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
      }
    }
  }

  cout << "A total of " << matchnum << " tracklets produced plausible astrometric fits\n";
  outstream1 << "A total of " << matchnum << " tracklets produced plausible astrometric fits\n";
  
  if(matchnum>0) {
    // Repeat the best fit
    matchct=bestmatch;
    pairct = matching_trkind[matchct];
    trkvec={};
    trkvec = tracklet_lookup(trk2det, pairct);
    trkvec_temp={};
    cout << "Attempting orbit fit for the best match: number " << matchct << ", with " << trkvec.size() << " points\n";
    observer_barypos2 = {};
    obsMJD2 = {};
    obsRA2 = obsDec2 = sigastrom2 = {};
    outdetvec2 = {};
    obsct=tct=0;   
    while(obsct<long(obsMJD.size()) || tct<long(trkvec.size())) {
      if(obsct<long(obsMJD.size()) && (tct>=long(trkvec.size()) || obsMJD[obsct] < image_log[detvec[trkvec[tct]].image].MJD)) {
	// Next point in the time-ordered set is from the orginal observation vectors
	obsMJD2.push_back(obsMJD[obsct]);
	obsDec2.push_back(obsDec[obsct]);
	obsRA2.push_back(obsRA[obsct]);
	sigastrom2.push_back(sigastrom[obsct]);
	observer_barypos2.push_back(observer_barypos[obsct]);
	outdetvec2.push_back(outdetvec[obsct]);
	obsct++;
      } else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && (obsRA[obsct] != detvec[trkvec[tct]].RA || obsDec[obsct] != detvec[trkvec[tct]].Dec)) {
	// We have two points at the same time, but they are not identical
	// Keep the one from the original vectors
	obsMJD2.push_back(obsMJD[obsct]);
	obsDec2.push_back(obsDec[obsct]);
	obsRA2.push_back(obsRA[obsct]);
	sigastrom2.push_back(sigastrom[obsct]);
	observer_barypos2.push_back(observer_barypos[obsct]);
	outdetvec2.push_back(outdetvec[obsct]);
	obsct++;
	// And also keep the one from the new tracklet
	obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	obsRA2.push_back(detvec[trkvec[tct]].RA);
	obsDec2.push_back(detvec[trkvec[tct]].Dec);
	sigastrom2.push_back(1.0);
	// Convert input heliocentric observer position to barycentric.
	planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	observer_barypos2.push_back(tppos);
	outdetvec2.push_back(detvec[trkvec[tct]]);
	trkvec_temp.push_back(trkvec[tct]);
	tct++;
      } else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && obsRA[obsct] == detvec[trkvec[tct]].RA && obsDec[obsct] == detvec[trkvec[tct]].Dec ) {
	// We have two exactly identical points. Keep the one from the original vectors...
	obsMJD2.push_back(obsMJD[obsct]);
	obsDec2.push_back(obsDec[obsct]);
	obsRA2.push_back(obsRA[obsct]);
	sigastrom2.push_back(sigastrom[obsct]);
	observer_barypos2.push_back(observer_barypos[obsct]);
	outdetvec2.push_back(outdetvec[obsct]);
	obsct++;
	// ...and discard the one from the new tracklet
	tct++;
      } else if (tct<long(trkvec.size())) {
	// Next point in the time-ordered set is from the new tracklet
	obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	obsRA2.push_back(detvec[trkvec[tct]].RA);
	obsDec2.push_back(detvec[trkvec[tct]].Dec);
	sigastrom2.push_back(1.0);
	// Convert input heliocentric observer position to barycentric.
	planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	observer_barypos2.push_back(tppos);
	outdetvec2.push_back(detvec[trkvec[tct]]);
	trkvec_temp.push_back(trkvec[tct]);
	tct++;
      } else {
	cerr << "Error: logically excluded case in loading of new vectors\n";
	return(5);
      }
    }
    if(trkvec_temp.size()>0) {
      tracklets_added++;
      for(tct=0;tct<long(trkvec_temp.size());tct++) {
	detections_added.push_back(trkvec_temp[tct]);
      }
    }
    
    // Load the augmented vectors back into the base vectors, so we can augment again (or output).
    observer_barypos = observer_barypos2;
    obsMJD = obsMJD2;
    obsRA = obsRA2;
    obsDec = obsDec2;
    sigastrom = sigastrom2;
    outdetvec = outdetvec2;
    // Perform orbit-fit to augmented input data
    cout << "Launching arc6D01 on best new tracklet " << matchct+1 << " with " << obsMJD.size() << " data points\n";
    outstream1 << "Launching arc6D01 on best new tracklet " << matchct+1  << " = " << detvec[trkvec[0]].idstring << " at MJD " << detvec[trkvec[0]].MJD << ", with " << trkvec.size() << " tracklet points and hence " << obsMJD.size() << " total data points\n";
    arc6D01(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,minchichange,planetfile_refpoint,startpos,startvel,bestRA,bestDec,bestresid,outpos, outvel, &bestchi, &astromRMS);
    outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
    // Save best-fit state vectors for new fits.
    startpos = outpos;
    startvel = outvel;
  }
  
  while(matchnum>1) {
    // Multiple possible matches existed. Try to augment the tracklet again.
    // Remove from matching_trkind the tracklet we've already used.
    pairct = matching_trkind[bestmatch];
    cout << "Best fit was to potential match " << bestmatch+1 << " = pairct " << pairct << "\n";
    outstream1 << "Best fit was to potential match " << bestmatch+1 << " = pairct " << pairct << "\n";
    matching_trkind={};
    for(matchct=0;matchct<long(matching_trkind2.size());matchct++) {
      if(matching_trkind2[matchct]!=pairct) matching_trkind.push_back(matching_trkind2[matchct]);
    }
    cout << "Eliminating this best match, " << matching_trkind.size() << " tracklets with plausible fits remain\n";
    outstream1 << "Eliminating this best match, " << matching_trkind.size() << " tracklets with plausible fits remain\n";
    // Loop over remaining potential matches.
    matchct=0;
    matchnum=0;
    bestmatch = -1;
    bestRMS=LARGERR;
    matching_trkind2 = {};
    for(matchct=0;matchct<long(matching_trkind.size());matchct++) {
      pairct = matching_trkind[matchct];
      trkvec={};
      trkvec = tracklet_lookup(trk2det, pairct);
      cout << "Attempting orbit fit for potential match number " << matchct+1 << ", with " << trkvec.size() << " points\n";
      observer_barypos2 = {};
      obsMJD2 = {};
      obsRA2 = obsDec2 = sigastrom2 = {};
      outdetvec2 = {};
      obsct=tct=0;
      while(obsct<long(obsMJD.size()) || tct<long(trkvec.size())) {
        if(obsct<long(obsMJD.size()) && (tct>=long(trkvec.size()) || obsMJD[obsct] < image_log[detvec[trkvec[tct]].image].MJD)) {
	  // Next point in the time-ordered set is from the orginal observation vectors
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	} else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && (obsRA[obsct] != detvec[trkvec[tct]].RA || obsDec[obsct] != detvec[trkvec[tct]].Dec)) {
	  // We have two points at the same time, but they are not identical
	  // Keep the one from the original vectors
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	  // And also keep the one from the new tracklet
	  obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	  obsRA2.push_back(detvec[trkvec[tct]].RA);
	  obsDec2.push_back(detvec[trkvec[tct]].Dec);
	  sigastrom2.push_back(1.0);
	  // Convert input heliocentric observer position to barycentric.
	  planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	  tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	  observer_barypos2.push_back(tppos);
	  outdetvec2.push_back(detvec[trkvec[tct]]);
	  tct++;
	} else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && obsRA[obsct] == detvec[trkvec[tct]].RA && obsDec[obsct] == detvec[trkvec[tct]].Dec ) {
	  // We have two exactly identical points. Keep the one from the original vectors...
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	  // ...and discard the one from the new tracklet
	  tct++;
	} else if (tct<long(trkvec.size())) {
	  // Next point in the time-ordered set is from the new tracklet
	  obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	  obsRA2.push_back(detvec[trkvec[tct]].RA);
	  obsDec2.push_back(detvec[trkvec[tct]].Dec);
	  sigastrom2.push_back(1.0);
	  // Convert input heliocentric observer position to barycentric.
	  planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	  tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	  observer_barypos2.push_back(tppos);
	  outdetvec2.push_back(detvec[trkvec[tct]]);
	  tct++;
	} else {
	  cerr << "Error: logically excluded case in loading of new vectors\n";
	  return(5);
	}
      }
     // Perform orbit-fit to augmented input data
      cout << "Launching arc6D01 for new-iteration match " << matchct+1 << " of " << matching_trkind.size() << ", with " << obsMJD2.size() << " data points\n";
      outstream1 << "Launching arc6D01 for new-iteration match " << matchct+1 << " = " << detvec[trkvec[0]].idstring  << " of " << matching_trkind.size() << " at MJD " << detvec[trkvec[0]].MJD << ", with " << trkvec.size() << " tracklet points and hence " << obsMJD2.size() << " total data points\n";
      arc6D01(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_barypos2,obsMJD2,obsRA2,obsDec2,sigastrom2,minchichange,planetfile_refpoint,startpos,startvel,bestRA,bestDec,bestresid,outpos, outvel, &bestchi, &astromRMS);
      outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
      if(astromRMS<max_astrom_rms) {
	// The fit was good
	matchnum+=1;
	matching_trkind2.push_back(matching_trkind[matchct]);
	outstream1 << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
	if(astromRMS<bestRMS) {
	  if(bestRMS>max_astrom_rms) {
	    outstream1 << "the first good fit to be identified in this pass\n";
	  } else {
	    outstream1 << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
	  }
	  bestRMS = astromRMS;
	  bestmatch = matchct;
	} else {
	outstream1 << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
	}
      }
    }
    if(matchnum>0) {
      // Repeat the best fit
      matchct=bestmatch;
      pairct = matching_trkind[matchct];
      trkvec={};
      trkvec = tracklet_lookup(trk2det, pairct);
      trkvec_temp={};
      cout << "Attempting orbit fit for potential match number " << matchct << ", with " << trkvec.size() << " points\n";
      observer_barypos2 = {};
      obsMJD2 = {};
      obsRA2 = obsDec2 = sigastrom2 = {};
      outdetvec2 = {};
      obsct=tct=0;
      while(obsct<long(obsMJD.size()) || tct<long(trkvec.size())) {
        if(obsct<long(obsMJD.size()) && (tct>=long(trkvec.size()) || obsMJD[obsct] < image_log[detvec[trkvec[tct]].image].MJD)) {
	  // Next point in the time-ordered set is from the orginal observation vectors
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	} else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && (obsRA[obsct] != detvec[trkvec[tct]].RA || obsDec[obsct] != detvec[trkvec[tct]].Dec)) {
	  // We have two points at the same time, but they are not identical
	  // Keep the one from the original vectors
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	  // And also keep the one from the new tracklet
	  obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	  obsRA2.push_back(detvec[trkvec[tct]].RA);
	  obsDec2.push_back(detvec[trkvec[tct]].Dec);
	  sigastrom2.push_back(1.0);
	  // Convert input heliocentric observer position to barycentric.
	  planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	  tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	  observer_barypos2.push_back(tppos);
	  outdetvec2.push_back(detvec[trkvec[tct]]);
	  trkvec_temp.push_back(trkvec[tct]);
	  tct++;
	} else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && obsRA[obsct] == detvec[trkvec[tct]].RA && obsDec[obsct] == detvec[trkvec[tct]].Dec ) {
	  // We have two exactly identical points. Keep the one from the original vectors...
	  obsMJD2.push_back(obsMJD[obsct]);
	  obsDec2.push_back(obsDec[obsct]);
	  obsRA2.push_back(obsRA[obsct]);
	  sigastrom2.push_back(sigastrom[obsct]);
	  observer_barypos2.push_back(observer_barypos[obsct]);
	  outdetvec2.push_back(outdetvec[obsct]);
	  obsct++;
	  // ...and discard the one from the new tracklet
	  tct++;
	} else if (tct<long(trkvec.size())) {
	  // Next point in the time-ordered set is from the new tracklet
	  obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
	  obsRA2.push_back(detvec[trkvec[tct]].RA);
	  obsDec2.push_back(detvec[trkvec[tct]].Dec);
	  sigastrom2.push_back(1.0);
	  // Convert input heliocentric observer position to barycentric.
	  planetposvel01LD(image_log[detvec[trkvec[tct]].image].MJD, polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
	  tppos = point3LD(image_log[detvec[trkvec[tct]].image].X + outpos.x, image_log[detvec[trkvec[tct]].image].Y + outpos.y, image_log[detvec[trkvec[tct]].image].Z + outpos.z);
	  observer_barypos2.push_back(tppos);
	  outdetvec2.push_back(detvec[trkvec[tct]]);
	  trkvec_temp.push_back(trkvec[tct]);
	  tct++;
	} else {
	  cerr << "Error: logically excluded case in loading of new vectors\n";
	  return(5);
	}
      }
      if(trkvec_temp.size()>0) {
	tracklets_added++;
	for(tct=0;tct<long(trkvec_temp.size());tct++) {
	  detections_added.push_back(trkvec_temp[tct]);
	}
      }
      // Load the augmented vectors back into the base vectors, so we can augment again (or output).
      observer_barypos = observer_barypos2;
      obsMJD = obsMJD2;
      obsRA = obsRA2;
      obsDec = obsDec2;
      sigastrom = sigastrom2;
      outdetvec = outdetvec2;
      // Perform orbit-fit to augmented input data
      cout << "Launching arc6D01 on additional good tracklet " << matchct +1 << " with " << obsMJD.size() << " data points\n";
      outstream1 << "Launching arc6D01 on additional good tracklet " << matchct+1  << " = " << detvec[trkvec[0]].idstring << " at MJD " << detvec[trkvec[0]].MJD << ", with " << trkvec.size() << " tracklet points and hence " << obsMJD2.size() << " total data points\n";
      arc6D01(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,minchichange,planetfile_refpoint,startpos,startvel,bestRA,bestDec,bestresid,outpos, outvel, &bestchi, &astromRMS);
      outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
      startpos = outpos;
      startvel = outvel;
    } else {
      // No additional good match was found. Re-do last good fit.
      cout << "Launching arc6D01, reverting to old data set with " << obsMJD.size() << " data points\n";
      outstream1 << "Launching arc6D01, reverting to old data set with " << obsMJD.size() << " data points\n";
      arc6D01(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_barypos,obsMJD,obsRA,obsDec,sigastrom,minchichange,planetfile_refpoint,startpos,startvel,bestRA,bestDec,bestresid,outpos, outvel, &bestchi, &astromRMS);
      outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
      cout << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
    }
  }

  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromRMS << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  cout << "State vectors correspond to reference point " << planetfile_refpoint << " in the input planet files\n";
  
  outstream1 << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromRMS << "\n";
  outstream1 << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  outstream1 << "State vectors correspond to reference point " << planetfile_refpoint << " in the input planet files\n";

  for(obsct=0;obsct<long(obsMJD.size());obsct++) {
    outstream1 << fixed << setprecision(8) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << bestRA[obsct] << " " << bestDec[obsct] << fixed << setprecision(8) << " " << (obsRA[obsct]-bestRA[obsct])*cos(obsDec[obsct]/DEGPRAD)*3600.0l << " " << (obsDec[obsct]-bestDec[obsct])*3600.0l << " " << bestresid[obsct] << "\n";
    }
  outstream1 << "\n\n\n";
  for(i=0; i<long(detections_added.size()); i++) {
    outstream1 << fixed << setprecision(7) << detvec[detections_added[i]].MJD << "," << detvec[detections_added[i]].RA << "," << detvec[detections_added[i]].Dec << ",";
    outstream1 << fixed << setprecision(4) << detvec[detections_added[i]].mag << ",";
    outstream1 << fixed << setprecision(2) << detvec[detections_added[i]].trail_len << "," << detvec[detections_added[i]].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << detvec[detections_added[i]].sigmag << ",";
    outstream1 << fixed << setprecision(3) << detvec[detections_added[i]].sig_across << "," << detvec[detections_added[i]].sig_along << ",";
    outstream1 << detvec[detections_added[i]].image << "," << detvec[detections_added[i]].idstring << "," << detvec[detections_added[i]].band << ",";
    outstream1 << detvec[detections_added[i]].obscode << "," << detvec[detections_added[i]].known_obj << ",";
    outstream1 << detvec[detections_added[i]].det_qual << "," << detvec[detections_added[i]].index << "\n";
  }
  
  outstream1.close();
  // Write output file in hldet format
  sort(outdetvec.begin(), outdetvec.end(), early_hldet());
  outstream1.open(outfile);
  for(i=0;i<long(outdetvec.size());i++) {
    if(i==0 || outdetvec[i].MJD!=outdetvec[i-1].MJD || outdetvec[i].RA!=outdetvec[i-1].RA || outdetvec[i].Dec!=outdetvec[i-1].Dec) {
      outstream1 << fixed << setprecision(7) << outdetvec[i].MJD << "," << outdetvec[i].RA << "," << outdetvec[i].Dec << ",";
      outstream1 << fixed << setprecision(4) << outdetvec[i].mag << ",";
      outstream1 << fixed << setprecision(2) << outdetvec[i].trail_len << "," << outdetvec[i].trail_PA << ",";
      outstream1 << fixed << setprecision(4) << outdetvec[i].sigmag << ",";
      outstream1 << fixed << setprecision(3) << outdetvec[i].sig_across << "," << outdetvec[i].sig_along << ",";
      outstream1 << outdetvec[i].image << "," << outdetvec[i].idstring << "," << outdetvec[i].band << ",";
      outstream1 << outdetvec[i].obscode << "," << outdetvec[i].known_obj << ",";
      outstream1 << outdetvec[i].det_qual << "," << outdetvec[i].index << "\n";
    }
  }
  outstream1.close();

  cout << tracklets_added << " tracklets, totalling " << detections_added.size() << " individual points,  were added to the orbit\n";

  return(0);
}
