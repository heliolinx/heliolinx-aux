// tracklet_arctrace02g.cpp: April 16, 2025:
// Like tracklet_arctrace02f.cpp, but has some innovations to
// enable it to run faster.
//
// Description of ancestor program tracklet_arctrace02f.cpp:
// Significant advance on tracklet_arctrace02e.cpp, this program
// will fit potentially matching tracklets two at a time.
// Outline of how this works: The first part of the program, in
// which all potentially matching input tracklets are orbit-fit
// one a time, proceeds exactly as before. A list of tracklets
// yielding acceptable orbit fits (if any) is constructed.
// At this point the processing diverges from previous programs
// of this type. Given a user-selectable integer twinfitnum,
// which defaults to 10, the twinfitnum best-fitting tracklets
// will be paired successively with each of the remaining tracklets.
// The previous program paired only the very best-fitting tracklet
// with additional tracklets (equivalent to twinfitnum=1 in the
// current program). This made the implicit assumption that if
// any of the input tracklets were true matches, one of these
// true matches must yield the lowest-RMS orbit fit. This assumption
// was often correct, but has been found to fail -- meaning that,
// although a good fit existed, it would not be found. The current
// program is intended to avoid this problem and ensure that any
// viable, multi-tracklet fit that exists will be found. 
//
// Description of ancestor program tracklet_arctrace02e.cpp:
// Refactored to use new add_tracklet01() function. Not a big
// change, just making the code simpler and shorter.
//
// tracklet_arctrace02d.cpp: March 31, 2025:
// testing new fitting code, e.g. arctrace03().
//
// March 26, 2025: tracklet_arctrace02c.cpp:
// Like tracklet_arctrace02b.cpp, but replaces the calls to arc6D01()
// with a calls to arctrace02(), which does a Keplerian pre-fit, and
// appears to be more accurate and equally fast.
//
// Description of ancestor program tracklet_arctrace02b.cpp;
// Like tracklet_arctrace02a.cpp, but implements three improvements
// originally developed in the tracklet_arctrace01 series b, c, and d
// programs. These improvements are, first, interpolation of the
// ephemeris vectors (tracklet_arctrace01b.cpp). Second, incremental
// de-duplication of the augmented arcs (tracklet_arctrace01c.cpp),
// and finally, the use of a k-d tree range query for the tracklet
// matching (tracklet_arctrace01d.cpp).
//
// Description of ancestor program tracklet_arctrace02a.cpp
// Probes multiple MCMC realizations, supplied via a representative
// state vector file produced by, e.g., MCMC_arctrace01g.
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
#define DEBUG 1
#define TIMERANGEFAC 2000000000
#define TWINFIT_DEFAULT 10

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain
// MJD, observer X, Y, Z, VX, VY, VZ, RA, and Dec. The observer coordinates and velocities must
// be in km and km/sec, relative to the Sun. The RA and Dec must be in decimal degrees.
static void show_usage()
{
  cerr << "Usage: tracklet_arctrace02a -cfg configfile -observations obsfile -obscode obscodefile -repvec representative_statevector_file -minchi min_chi_change -kepspan time_span_for_Keplerian_fit(day) -twinfitnum twinfitnum -allfitnum allfitnum -mjdstart mjdstart -mjdend mjdend -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file  -timetol MJD_matching_tolerance(days) -skytol sky_matching_radius(deg) -veltol velocity_matching_radius(deg/day) -max_astrom_rms max astrometric RMS (arcsec) -outfile outfile -logfile logfile -verbose verbosity\n";
}

int add_tracklet01(const vector <point3LD> &observer_heliopos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <hldet> &outdetvec, const vector <long> &trkvec, const vector <hldet> &detvec, const vector <hlimage> &image_log, vector <point3LD> &observer_heliopos2, vector <long double> &obsMJD2, vector <double> &obsRA2, vector <double> &obsDec2, vector <double> &sigastrom2, vector <hldet> &outdetvec2, vector <long> &trkvec_temp);

int main(int argc, char *argv[])
{
  ifstream instream1;
  long i, j, k;
  int status, planetnum, planetct, pctSun, pctEarth, configread;
  long badread,reachedeof,obsnum,obsct,tct;
  i=j=status=configread=badread=reachedeof=obsnum=obsct=tct=0;
  int polyorder = 3;
  vector <observatory> observatory_list = {};
  string obscodefile;
  vector <hldet> obsdetvec = {};
  vector <hldet> outdetvec = {};
  vector <hldet> outdetvec2 = {};
  vector <hldet> outdetvec3 = {};
  string stest;
  string configfile;
  string obsfile;
  string planetfile;
  string logfile;
  string outfile;
  string repfile;
  double kepspan=50.0;
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
  vector <long double> obsMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> sigastrom;
  vector <point3LD> observer_heliopos2;
  vector <long double> obsMJD2;
  vector <double> obsRA2;
  vector <double> obsDec2;
  vector <double> sigastrom2;
  vector <point3LD> observer_heliopos3;
  vector <long double> obsMJD3;
  vector <double> obsRA3;
  vector <double> obsDec3;
  vector <double> sigastrom3;
  string lnfromfile;
  vector <string> linestringvec;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  vector <point3LD> startposvec;
  vector <point3LD> startvelvec;
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
  long double max_astrom_rms = 1.0l;
  long double mjd_statevecs = 0.0l;
  vector <long double> mjd_sv_vec;
  long planetfile_refpoint=0;
  long double mjdstart = 0;
  long double mjdend = 0;
  long planetfile_startpoint,planetfile_endpoint;
  vector <point3LD> orbit05pos;
  vector <point3LD> orbit05vel;
  vector <long double> orbit05MJD;
  long double light_travel_time = 0.0l;
  long double timediff = 0.0l;
  long ephct=0;

  string pairdetfile="pairdetfile01.csv";
  string trackletfile="trackletfile01.csv";
  string trk2detfile="trk2detfile01.csv";
  string imfile;
  vector <hldet> detvec = {};
  vector <hlimage> image_log;
  vector <tracklet> tracklets;
  vector <longpair> trk2det;
  long imnum,pairnum,pairct,pairct2,i1,i2,detnum;
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
  vector <vector <double>> ephRAmat;
  vector <vector <double>> ephDecmat;
  vector <vector <double>> ephRAvelmat;
  vector <vector <double>> ephDecvelmat;
  double RAvel,Decvel;
  double delem = 0.0;
  long ind = 0;
  double_index one_dind = double_index(delem,ind);
  vector <double_index> trk_dindvec;

  vector <long> matching_trkind;
  vector <long> matching_trkind2;
  vector <long> trkvec;
  vector <long> trkvec2;
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
  vector <point3LD> obsposvec;
  vector <point3LD> obsvelvec;
  long startpoint,endpoint;
  long repnum,repct;
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
  int refpoint=0;
  long twinfitnum = TWINFIT_DEFAULT;
  long twinfitct = 0;
  long allfitnum = 0;
  int allfitnum_default = 1;
  double_index one_dindex = double_index(0.0, 0);
  vector <double_index> track_rmsvec;
  vector <double_index> track_rmsvec2;
  vector <long double> statevec;
  vector <vector <long double>> statevecmat;
  long double bestRMS=LARGERR;
  long bp1,bp2,bp3,bestprimary;
  bp1=bp2=bp3=bestprimary=-1;
  vector <vector <int>> twinfit_indices;
  vector <int> onefit_indices;
  
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
    } else if(string(argv[i]) == "-repstate" || string(argv[i]) == "-rep"  || string(argv[i]) == "-repvec" || string(argv[i]) == "-repfile") {
      if(i+1 < argc) {
	//There is still something to read;
	repfile=argv[++i];
	i++;
      }
      else {
	cerr << "Representative state vector file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-twinfitnum") {
      if(i+1 < argc) {
	//There is still something to read;
	twinfitnum=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Twin-tracklet fitting number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-allfitnum") {
      if(i+1 < argc) {
	//There is still something to read;
	allfitnum=stol(argv[++i]);
	allfitnum_default=0;
	i++;
      }
      else {
	cerr << "Maximum re-fitting number keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
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

  // Read the representative state vector file
  instream1.open(repfile);
  if(!instream1) {
    cerr << "can't open input file " << repfile << "\n";
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
	try { mjd_statevecs = stold(stest); }
	catch(...) { cerr << "ERROR: cannot read MJD " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the state vector x component
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { startpos.x = stold(stest); }
	catch(...) { cerr << "ERROR: cannot read state vector x component " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the state vector y component
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { startpos.y = stold(stest); }
	catch(...) { cerr << "ERROR: cannot read state vector y component " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the state vector z component
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { startpos.z = stold(stest); }
	catch(...) { cerr << "ERROR: cannot read state vector z component " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the state vector vx component
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { startvel.x = stold(stest); }
	catch(...) { cerr << "ERROR: cannot read state vector vx component " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the state vector vy component
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { startvel.y = stold(stest); }
	catch(...) { cerr << "ERROR: cannot read state vector vy component " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the state vector vz component
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { startvel.z = stold(stest); }
	catch(...) { cerr << "ERROR: cannot read state vector vz component " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      if(badread==0) {
	mjd_sv_vec.push_back(mjd_statevecs);
	startposvec.push_back(startpos);
	startvelvec.push_back(startvel);
      }
    } else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(badread!=0) {
      cerr << "ERROR reading representative state vector file " << repfile << "\n";
      cerr << "Last point was " << startposvec.size() << "; last file line was " << lnfromfile << "\n";
      return(badread);
    }
  }
  instream1.close();

  if(badread!=0) {
    cerr << "ERROR reading representative state vector file " << repfile << "\n";
    return(badread);
  } 
  if(reachedeof==1) { 
    if(verbose>=1) cout << "Input file " << repfile << " read successfully to the end.\n";
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << repfile << " before the end\n";
  } else if(reachedeof==-1) {
    cerr << "ERROR: file read failed\n";
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
  }
  repnum = startposvec.size();
  if(repnum!=long(startvelvec.size())) {
    cerr << "ERROR: vector mismatch 1 from representative state vector file\n";
    cerr << repnum << " != " << startvelvec.size() << "\n";
    return(2);
  }
  if(repnum!=long(mjd_sv_vec.size())) {
    cerr << "ERROR: vector mismatch 2 from representative state vector file\n";
    cerr << repnum << " != " << mjd_sv_vec.size() << "\n";
    return(2);
  }
  cout << "Read " << repnum << " lines from representative state vector file " << repfile << "\n";
 
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
  observer_heliopos = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    status = obscode_lookup(observatory_list,obsdetvec[obsct].obscode,obslon,plxcos,plxsin);
    if(status>0) {
      cerr << "ERROR: obscode_lookup failed for observatory code " << obsdetvec[obsct].obscode << "\n";
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
  cout << "MJD range " << mjdstart << " to " << mjdend << " for ephemeris corresponds to planet file index range from " << planetfile_startpoint << " to " << planetfile_endpoint << "\n";
  
  // Attempt orbit-fit using arctrace02()
  arctrace02(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_heliopos,obsMJD,obsRA,obsDec,sigastrom,kepspan,minchichange, bestRA, bestDec, bestresid, outpos, outvel, &bestchi, &astromRMS, &refpoint, verbose);
  
  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromRMS << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  cout << "State vectors correspond to reference point " << refpoint << " in the input planet files\n";

  // FINISHED WITH ORBIT FIT TO ORIGINAL INPUT DATA

  // CALCULATE EPHEMERIDES BASED ON INPUT STATE VECTORS
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
  earliest_mjd = ephMJD[0];
  latest_mjd = ephMJD[ephMJD.size()-1];

  // Load observer (Earth) positions corresponding to ephMJD
  obsposvec = obsvelvec = {};
  for(ephct=0;ephct<long(ephMJD.size());ephct++) {
    // Interpolate to get coordinates of observer (Earth).
    planetposvel01LD(ephMJD[ephct], polyorder, planetmjd, Earthpos, Earthvel, obspos, obsvel);
    obsposvec.push_back(obspos);
    obsvelvec.push_back(obsvel);
  }

  for(repct=0;repct<repnum;repct++) {
    cout << "Calculating representative ephemeris " << repct << "\n";
    ephRA = {};
    ephDec = {};
    ephRAvel = {};
    ephDecvel = {};

    mjd_statevecs = mjd_sv_vec[repct];
    startpos = startposvec[repct];
    startvel = startvelvec[repct];
    planetfile_refpoint=-99;
    j=0;
    while(j<long(planetmjd.size()) && planetfile_refpoint==-99) {
      j++;
      if(fabs(planetmjd[j]-mjd_statevecs) <= IMAGETIMETOL/SOLARDAY) planetfile_refpoint=j;
    }
    if(planetfile_refpoint>=0) {
      //cout << "Input mjd " << mjd_statevecs << " corresponds to planet file entry number " << planetfile_refpoint << "\n";
    } else {
      cerr << "ERROR: input mjd " << mjd_statevecs << " does not match any entry in the planet file\n";
      return(1);
    }

    integrate_orbit05LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, startpos, startvel, planetfile_startpoint, planetfile_refpoint, planetfile_endpoint, orbit05MJD, orbit05pos, orbit05vel);
  
    for(ephct=0;ephct<long(ephMJD.size());ephct++) {
      // Load observer (Earth) coords from previously calculated vectors
      obspos = obsposvec[ephct];
      obsvel = obsvelvec[ephct];
      // Interpolate to get coordinates of target object
      planetposvel01LD(ephMJD[ephct], polyorder, orbit05MJD, orbit05pos, orbit05vel, targpos, targvel);
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
    ephRAmat.push_back(ephRA);
    ephDecmat.push_back(ephDec);
    ephRAvelmat.push_back(ephRAvel);
    ephDecvelmat.push_back(ephDecvel);
  }
  // FINISHED WITH EPHEMERIS CALCULATION.
  cout << "Full matrix sizes: " << ephRAmat.size() << " " << ephDecmat.size() << " " << ephRAvelmat.size() << " " << ephDecvelmat.size() << "\n";
  cout << "Final row sizes: " << ephRAmat[repnum-1].size() << " " << ephDecmat[repnum-1].size() << " " << ephRAvelmat[repnum-1].size() << " " << ephDecvelmat[repnum-1].size() << " " << ephMJD.size() << "\n";

  for(i=0;i<long(ephMJD.size());i++) {
    double maxdist=0;
    double maxvdist=0;
    double maxra,minra,maxravel,minravel;
    long maxrapt,minrapt,maxvelpt,minvelpt;
    maxrapt = minrapt = maxvelpt = minvelpt = 0;
    maxra = minra = ephRAmat[0][i];
    maxravel = minravel = ephRAvelmat[0][i];
    for(repct=1;repct<repnum;repct++) {
      if(ephRAmat[repct][i]>maxra) {
	maxra = ephRAmat[repct][i];
	maxrapt = repct;
      }
      if(ephRAmat[repct][i]<minra) {
	minra = ephRAmat[repct][i];
	minrapt = repct;
      }
       if(ephRAvelmat[repct][i]>maxravel) {
	maxravel = ephRAvelmat[repct][i];
	maxvelpt = repct;
      }
      if(ephRAvelmat[repct][i]<minravel) {
	minravel = ephRAvelmat[repct][i];
	minvelpt = repct;
      }
    }
    for(repct=0;repct<repnum;repct++) {
      if(repct!=maxrapt && repct!=minrapt) {
	dist = distradec01(ephRAmat[maxrapt][i],ephDecmat[maxrapt][i],ephRAmat[repct][i],ephDecmat[repct][i]);
	if(dist>maxdist) maxdist=dist;
	dist = distradec01(ephRAmat[minrapt][i],ephDecmat[minrapt][i],ephRAmat[repct][i],ephDecmat[repct][i]);
	if(dist>maxdist) maxdist=dist;
      }
      if(repct!=maxvelpt && repct!=minvelpt) {
	dist = sqrt(DSQUARE(ephRAvelmat[maxvelpt][i]-ephRAvelmat[repct][i]) + DSQUARE(ephDecvelmat[maxvelpt][i]-ephDecvelmat[repct][i]));
	if(dist>maxvdist) maxvdist=dist;
	dist = sqrt(DSQUARE(ephRAvelmat[minvelpt][i]-ephRAvelmat[repct][i]) + DSQUARE(ephDecvelmat[minvelpt][i]-ephDecvelmat[repct][i]));
	if(dist>maxvdist) maxvdist=dist;
      }
    }
    cout << "Ephemeris point " << i << ", MJD " << ephMJD[i] << ", sky spread is " << maxdist << " deg, and velocity spread is " << maxvdist << " deg/day\n";
  }

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
  outstream1.open(logfile);
  cout << "Preparing to search for matches among " << trk_dindvec.size() << " tracklets at " << ephMJD.size() << " different times, for " << repnum << " orbit clones\n";
  for(i=0;i<long(ephMJD.size());i++) {
    cout << "i,j,trk_dindvec[j].delem,ephMJD[i],timetol " << i << " " << j << " " << trk_dindvec[j].delem << " " << ephMJD[i] << " " << timetol << "\n";
    for(repct=0;repct<repnum;repct++) {
    projpt = celeproj01(ephRAmat[repct][i], ephDecmat[repct][i]);
    querypoint = point6ix2(int((ephMJD[i] - central_mjd)/timeintscale), int(projpt.x*DEGPRAD/posintscale + 0.5), int(projpt.y*DEGPRAD/posintscale + 0.5), int(projpt.z*DEGPRAD/posintscale + 0.5), int(ephRAvelmat[repct][i]/velintscale + 0.5), int(ephDecvelmat[repct][i]/velintscale + 0.5),pairct,i);
    queryout = {};
    kdrange_6i01(kdtree, querypoint, queryrange, queryout);
    // Loop over any matching tracklets
    for(k=0;k<long(queryout.size());k++) {
      pairct = kdtree[queryout[k]].point.i1; // Index to trkMJD, trkRA, trkDec, trkRAvel, and trkDecvel
      dist = distradec01(ephRAmat[repct][i],ephDecmat[repct][i],trkRA[pairct],trkDec[pairct]);
      vdist = sqrt(DSQUARE(ephRAvelmat[repct][i]-trkRAvel[pairct])+DSQUARE(ephDecvelmat[repct][i]-trkDecvel[pairct]));
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
	cout << fixed << setprecision(4) << "Match found for on tracklet line " << pairct << " for representative ephemeris " << repct << " : " << trkRA[pairct] << " " << trkDec[pairct] << " " << trkRAvel[pairct] << " " << trkDecvel[pairct] << ":\n";
	cout << fixed << setprecision(4) << detvec[trkvec[0]].idstring << " dist = " << dist << " degrees, vdist = " << vdist << " deg/day :\n" << ephMJD[i] << " " << ephRAmat[repct][i] << " " << ephDecmat[repct][i] << " " << ephRAvelmat[repct][i] << " " << ephDecvelmat[repct][i] << "\n";
	outstream1 << fixed << setprecision(4) << "Match found on tracklet line " << pairct << " for representative ephemeris " << repct << " : " << trkRA[pairct] << " " << trkDec[pairct] << " " << trkRAvel[pairct] << " " << trkDecvel[pairct] << ":\n";
	outstream1 << fixed << setprecision(4) << detvec[trkvec[0]].idstring << " dist = " << dist << " degrees, vdist = " << vdist << " deg/day :\n" << ephMJD[i] << " " << ephRAmat[repct][i] << " " << ephDecmat[repct][i] << " " << ephRAvelmat[repct][i] << " " << ephDecvelmat[repct][i] << "\n";
	matching_trkind.push_back(pairct);
      }
    }
    }
  }

  cout << "Prior to de-duplication, we have " << matching_trkind.size() << " total matching tracklets\n";
  
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
  track_rmsvec = {};
  for(matchct=0;matchct<long(matching_trkind.size());matchct++) {
    pairct = matching_trkind[matchct];
    trkvec={};
    trkvec = tracklet_lookup(trk2det, pairct);
    status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
    if(status!=0) {
      cerr << "ERROR: add_tracklet01 returned status " << status << "\n";
      return(status);
    }
    if(trkvec_temp.size()>0) {
      cout << "Attempting orbit fit for potential match number " << matchct+1 << ", with " << trkvec.size() << " points\n";  
      // Perform orbit-fit to augmented input data
      cout << "Launching arctrace03() for match " << matchct+1  << " of " << matching_trkind.size() << " with " << obsMJD2.size() << " data points\n";
      outstream1 << "Launching arctrace03() for match " << matchct+1  << " = " << detvec[trkvec[0]].idstring << " of " << matching_trkind.size() << " at MJD " << detvec[trkvec[0]].MJD << ", with " << trkvec.size() << " tracklet points and hence " << obsMJD2.size() << " total data points\n";
      arctrace03(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_heliopos2,obsMJD2,obsRA2,obsDec2,sigastrom2,kepspan,minchichange, bestRA, bestDec, bestresid, outpos, outvel, &bestchi, &astromRMS, &refpoint, verbose);
    
      outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
      if(astromRMS<bestRMS) bestRMS=astromRMS;
      if(astromRMS<max_astrom_rms) {
	// The fit was good
	matchnum+=1;
	one_dindex=double_index(astromRMS,matching_trkind[matchct]);
	track_rmsvec.push_back(one_dindex);
	outstream1 << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
      }
    } else {
      cout << "Potential match number " << matchct+1 << " had " << trkvec.size() << " points,\nbut all were redundant with observations already in the input linkage\n";  
    }
  }
  cout << "A total of " << matchnum << " individual tracklets produced plausible astrometric fits\n";
  outstream1 << "A total of " << matchnum << " individual tracklets produced plausible astrometric fits\n";
  cout << "Lowest astrometric RMS was " << bestRMS << "\n";
  outstream1 << "Lowest astrometric RMS was " << bestRMS << "\n";

  if(matchnum>1) {
    bestRMS=LARGERR;
    // Sort the tracklet RMS vector
    sort(track_rmsvec.begin(),track_rmsvec.end(),lower_double_index());
    for(matchct=0;matchct<matchnum;matchct++) {
      cout << "Match of rank " << matchct << ", index = " << track_rmsvec[matchct].index << ", RMS = " << track_rmsvec[matchct].delem << "\n";
    }
    bp1=bp2=bp3=bestprimary=-1;
    if(twinfitnum>matchnum) twinfitnum = matchnum;
    // Load empty index vector for each primary tracklet
    onefit_indices = {};
    twinfit_indices = {};
    for(twinfitct=0;twinfitct<twinfitnum;twinfitct++) twinfit_indices.push_back(onefit_indices);
    // Loop over the twinfitnum best matches.    
    for(twinfitct=0;twinfitct<twinfitnum;twinfitct++) {
      // Load the tracklet of rank twinfitct in the RMS-sorted list
      pairct = track_rmsvec[twinfitct].index;
      // Add this tracklet to the provisional fitting vectors.
      trkvec={};
      trkvec = tracklet_lookup(trk2det, pairct);
      status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
      if(status!=0) {
	cerr << "ERROR: add_tracklet01 returned status " << status << "\n";
	return(status);
      }
      if(allfitnum>matchnum || allfitnum_default==1) allfitnum = matchnum;
      // Loop over the allfitnum best tracklets after twinfit in the ordered list.
      for(matchct=twinfitct+1;matchct<allfitnum;matchct++) {
	pairct2 = track_rmsvec[matchct].index;
	trkvec2={};
	trkvec2 = tracklet_lookup(trk2det, pairct2);
	// Add this additional tracklet to the provisional fitting vectors
	status = add_tracklet01(observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec2, detvec, image_log, observer_heliopos3, obsMJD3, obsRA3, obsDec3, sigastrom3, outdetvec3, trkvec_temp);
	if(trkvec_temp.size()<=0) {
	  cout << "Tracklets " << twinfitct << " and " << matchct << " are apparently redundant\n";
	  outstream1 << "Tracklets " << twinfitct << " and " << matchct << " are apparently redundant\n";
	} else {
	  cout << "Attempting orbit fit for tracklet pair " << twinfitct << " (" << pairct << "), with " << trkvec.size() << " points and " << matchct << " (" << pairct2 << "), with " << trkvec2.size() << " points\n";
	  outstream1 << "Attempting orbit fit for tracklet pair " << twinfitct << " (" << pairct << "), with " << trkvec.size() << " points and " << matchct << " (" << pairct2 << "), with " << trkvec2.size() << " points\n";
	  arctrace03(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_heliopos3,obsMJD3,obsRA3,obsDec3,sigastrom3,kepspan,minchichange, bestRA, bestDec, bestresid, outpos, outvel, &bestchi, &astromRMS, &refpoint, verbose);
	  cout << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
	  outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
	  if(astromRMS<max_astrom_rms) {
	    // The fit is good.
	    twinfit_indices[twinfitct].push_back(pairct2);
	    if(matchct<twinfitnum) twinfit_indices[matchct].push_back(pairct);
	    cout << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
	    outstream1 << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
	    if(DEBUG>0) {
	      cout << "Loaded index " << pairct2 << " into twinfit_indices[" << twinfitct << "], which now has a length of " << twinfit_indices[twinfitct].size() << "\n";
	      outstream1 << "Loaded index " << pairct2 << " into twinfit_indices[" << twinfitct << "], which now has a length of " << twinfit_indices[twinfitct].size() << "\n";
	      if(matchct<twinfitnum) {
		cout << "Also loaded index " << pairct << " into twinfit_indices[" << matchct << "], which now has a length of " << twinfit_indices[matchct].size() << "\n";
		outstream1 << "Also loaded index " << pairct << " into twinfit_indices[" << matchct << "], which now has a length of " << twinfit_indices[matchct].size() << "\n";
	      }
	    }
	    if(astromRMS<bestRMS) {
	      if(bestRMS>max_astrom_rms) {
		cout << "the first good fit to be identified in this pass\n";
		outstream1 << "the first good fit to be identified in this pass\n";
	      } else {
		cout << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
		outstream1 << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
	      }
	      bestRMS = astromRMS;
	      bp1 = pairct;
	      bp2 = pairct2;
	      bestprimary=twinfitct;
	    } else {
	      cout << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
	      outstream1 << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
	    }
	    // Close loop on if-statement confirming we have a good fit.
	  } else if(astromRMS<bestRMS) {
	    // Record this as the best RMS so far, even though it
	    // does not qualify as a good fit.
	    bestRMS = astromRMS;
	  }
	  // Close loop on case checking that tracklets twinfitct and matchct were non-redundant.
	}
	// Close loop over secondary tracklets.
      }
      // Close loop over primary tracklets.
    }
    if(bp1>=0 && bp2>=0) {
      if(DEBUG>0) {
	cout << "bp1, bp2, bestprimary: " << bp1 << " " << bp2 << " " << bestprimary << "\n";
	outstream1 << "bp1, bp2, bestprimary: " << bp1 << " " << bp2 << " " << bestprimary << "\n";
	for(i=0;i<long(twinfit_indices.size());i++) {
	  cout << "Fitting indices for primary tracklet " << i << ": ";
	  for(j=0;j<long(twinfit_indices[i].size());j++) cout << twinfit_indices[i][j] << " ";
	  cout << "\n";
	  outstream1 << "Fitting indices for primary tracklet " << i << ": ";
	  for(j=0;j<long(twinfit_indices[i].size());j++) outstream1 << twinfit_indices[i][j] << " ";
	  outstream1 << "\n";
	}
      }
      // At least one good fit was found in the pair-search.
      // Recover the best pair of tracklets
      trkvec={};
      trkvec = tracklet_lookup(trk2det, bp1);
      trkvec2={};
      trkvec2 = tracklet_lookup(trk2det, bp2);
      // Augment the original observation vectors with the first tracklet, as before
      status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
      // Since we now know this tracklet is good, also add it to the
      // detections_added vector that will be output
      tracklets_added++;
      for(tct=0;tct<long(trkvec_temp.size());tct++) {
	detections_added.push_back(trkvec_temp[tct]);
      }
      // Augment the original observation vectors with the second tracklet, as before
      status = add_tracklet01(observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec2, detvec, image_log, observer_heliopos3, obsMJD3, obsRA3, obsDec3, sigastrom3, outdetvec3, trkvec_temp);
      // Since we now know this tracklet is good, also add it to the
      // detections_added vector that will be output
      tracklets_added++;
      for(tct=0;tct<long(trkvec_temp.size());tct++) {
	detections_added.push_back(trkvec_temp[tct]);
      }
      // Load the augmented vectors back into the base vectors, so we can augment again (or output).
      observer_heliopos = observer_heliopos3;
      obsMJD = obsMJD3;
      obsRA = obsRA3;
      obsDec = obsDec3;
      sigastrom = sigastrom3;
      outdetvec = outdetvec3;
      // Construct a new version of track_rmsvec without the
      // already-added tracklets, and including only tracklets
      // that were successfully paired with bp1=bestprimary
      onefit_indices = twinfit_indices[bestprimary];
      track_rmsvec2 = track_rmsvec;
      track_rmsvec = {};
      for(i=0;i<long(track_rmsvec2.size());i++) {
	if(track_rmsvec2[i].index!=bp1 && track_rmsvec2[i].index!=bp2) {
	  // This tracklet has NOT already been added as bp1 or bp2.
	  // See if it's in onefit_indices, which here is a list of all
	  // the tracklets that produced a viable orbit fit when paired
	  // with primary tracklet bp1=bestprimary.
	  for(j=0;j<long(onefit_indices.size());j++) {
	    if(track_rmsvec2[i].index==onefit_indices[j]) {
	      // Tracklet is on the list: add it
	      track_rmsvec.push_back(track_rmsvec2[i]);
	    }
	  }
	}
      }
      // Attempt to add a single tracklet, which would be the third  
      bestRMS=LARGERR;
      matchnum=0;
      bp3=-1;
      for(matchct=0;matchct<long(track_rmsvec.size());matchct++) {
	pairct = track_rmsvec[matchct].index;
	// Add this tracklet to the current vectors
	trkvec={};
	trkvec = tracklet_lookup(trk2det, pairct);
	status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
	if(status!=0) {
	  cerr << "ERROR: add_tracklet01 returned status " << status << "\n";
	  return(status);
	}
	if(trkvec.size()>0) {
	  // The tracklet was new and non-redundant
	  cout << "Attempting orbit fit for a potential 3rd tracklet " << matchct << " (" << pairct << "), of " << track_rmsvec.size() << " with " << trkvec.size() << " points.\n";
	  outstream1 << "Attempting orbit fit for a potential 3rd tracklet " << matchct << " (" << pairct << "), of " << track_rmsvec.size() << " with " << trkvec.size() << " points.\n";
	  arctrace03(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_heliopos2,obsMJD2,obsRA2,obsDec2,sigastrom2,kepspan,minchichange, bestRA, bestDec, bestresid, outpos, outvel, &bestchi, &astromRMS, &refpoint, verbose);
	  cout << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
	  outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
	  if(astromRMS<max_astrom_rms) {
	    matchnum++;
	    // The fit is good.
	    cout << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
	    outstream1 << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
	    if(astromRMS<bestRMS) {
	      if(bestRMS>max_astrom_rms) {
		cout << "the first good fit to be identified in this pass\n";
		outstream1 << "the first good fit to be identified in this pass\n";
	      } else {
		cout << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
		outstream1 << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
	      }
	      bestRMS = astromRMS;
	      bp3 = pairct;
	    } else {
	      cout << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
	      outstream1 << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
	    }
	    // Close loop on if-statement confirming we have a good fit.
	  }
	  // Close loop on case checking that tracklet matchct was non-redundant.
	}
	// Close loop over tracklets.
      }
      if(bp3>=0) {
	// At least one good fit was found. Load the best one available.
	pairct = bp3;
	// Add this tracklet to the current vectors
	trkvec={};
	trkvec = tracklet_lookup(trk2det, pairct);
	status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
	// Since we now know this tracklet is good, also add it to the
	// detections_added vector that will be output
	tracklets_added++;
	for(tct=0;tct<long(trkvec_temp.size());tct++) {
	  detections_added.push_back(trkvec_temp[tct]);
	}
	// Load the augmented vectors back into the base vectors, so we can augment again (or output).
	observer_heliopos = observer_heliopos2;
	obsMJD = obsMJD2;
	obsRA = obsRA2;
	obsDec = obsDec2;
	sigastrom = sigastrom2;
	outdetvec = outdetvec2;
	while(matchnum>1) {
	  // Multiple possible matches existed. Try to augment the linkage again
	  // Remove from matching_trkind the tracklet we've already used.
	  // Construct a new version of track_rmsvec without the
	  // already-added tracklets.
	  track_rmsvec2 = track_rmsvec;
	  track_rmsvec = {};
	  for(i=0;i<long(track_rmsvec2.size());i++) {
	    if(track_rmsvec2[i].index!=bp3) {
	      track_rmsvec.push_back(track_rmsvec2[i]);
	    }
	  }
	  // Attempt to add a single additional tracklet      
	  bestRMS=LARGERR;
	  matchnum=0;
	  bp3=-1;
	  for(matchct=0;matchct<long(track_rmsvec.size());matchct++) {
	    pairct = track_rmsvec[matchct].index;
	    // Add this tracklet to the current vectors
	    trkvec={};
	    trkvec = tracklet_lookup(trk2det, pairct);
	    status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
	    if(status!=0) {
	      cerr << "ERROR: add_tracklet01 returned status " << status << "\n";
	      return(status);
	    }
	    if(trkvec.size()>0) {
	      // The tracklet was new and non-redundant
	      cout << "Attempting orbit fit for beyond-3rd-tracklet candidate " << matchct << " (" << pairct << "), of " << track_rmsvec.size() << " with " << trkvec.size() << " points.\n";
	      outstream1 << "Attempting orbit fit for beyond-3rd-tracklet candidate " << matchct << " (" << pairct << "), of " << track_rmsvec.size() << " with " << trkvec.size() << " points.\n";
	      arctrace03(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_heliopos2,obsMJD2,obsRA2,obsDec2,sigastrom2,kepspan,minchichange, bestRA, bestDec, bestresid, outpos, outvel, &bestchi, &astromRMS, &refpoint, verbose);
	      cout << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
	      outstream1 << "Fit complete, chisq = " << bestchi << ", astromRMS = " << astromRMS << "\n";
	      if(astromRMS<max_astrom_rms) {
		matchnum++;
		// The fit is good.
		cout << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
		outstream1 << "With astromRMS = " << astromRMS << ", this is a good fit: \n";
		if(astromRMS<bestRMS) {
		  if(bestRMS>max_astrom_rms) {
		    cout << "the first good fit to be identified in this pass\n";
		    outstream1 << "the first good fit to be identified in this pass\n";
		  } else {
		    cout << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
		    outstream1 << "this fit supplants the previous best, which had astroRMS = " << bestRMS << "\n";
		  }
		  bestRMS = astromRMS;
		  bp3 = pairct;
		} else {
		  cout << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
		  outstream1 << "however, it does not supplant the previous best, which had astroRMS = " << bestRMS << "\n";
		}
		// Close loop on if-statement confirming we have a good fit.
	      }
	      // Close loop on case checking that tracklet matchct was non-redundant.
	    }
	    // Close loop over tracklets.
	  }
	  if(bp3>=0) {
	    // At least one good fit was found. Load the best one available.
	    pairct = bp3;
	    // Add this tracklet to the current vectors
	    trkvec={};
	    trkvec = tracklet_lookup(trk2det, pairct);
	    status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
	    // Since we now know this tracklet is good, also add it to the
	    // detections_added vector that will be output
	    tracklets_added++;
	    for(tct=0;tct<long(trkvec_temp.size());tct++) {
	      detections_added.push_back(trkvec_temp[tct]);
	    }
	    // Load the augmented vectors back into the base vectors, so we can augment again (or output).
	    observer_heliopos = observer_heliopos2;
	    obsMJD = obsMJD2;
	    obsRA = obsRA2;
	    obsDec = obsDec2;
	    sigastrom = sigastrom2;
	    outdetvec = outdetvec2;
	    // Close if-statement checking we had at least one good fit in the inner loop.
	  } else {
	    // There were no good fits in the inner loop. We have added at least
	    // three tracklets, and maybe several more, but no more can be added now.
	    cout << "Exiting inner loop after successfully adding " << tracklets_added << " tracklets totalling " <<  detections_added.size() << " points\n";
	    outstream1 << "Exiting inner loop after successfully adding " << tracklets_added << " tracklets totalling " <<  detections_added.size() << " points\n";
	  }
	  // Close while loop on possibility of multiple fits
	}
	// Close if-statement checking we had at least one good fit in the outer loop.
      } else {
	// We added two tracklets in the pair fit, but could not add any more.
	cout << "Exiting outer loop: could only add 2 = " << tracklets_added << " tracklets totalling " <<  detections_added.size() << " points\n";
	outstream1 << "Exiting outer loop: could only add 2 = " << tracklets_added << " tracklets totalling " <<  detections_added.size() << " points\n";	
      }
      // Close if-statement checking that the tracklet pair fit succeeded
    } else {
      // We could not find any valid pair of tracklets to be added
      cout << "Exiting pair-fitting loop without success: best fit identified was still unacceptable, with RMS = " << bestRMS << " arcsec\n";
      outstream1 << "Exiting pair-fitting loop without success: best fit identified was still unacceptable, with RMS = " << bestRMS << " arcsec\n";
      // There were, however, valid single-tracklet extensions.
      // Add the best of these to the output vectors:
      pairct = track_rmsvec[0].index;
      trkvec={};
      trkvec = tracklet_lookup(trk2det, pairct);
      status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
      // Load the augmented vectors back into the base vectors for output
      observer_heliopos = observer_heliopos2;
      obsMJD = obsMJD2;
      obsRA = obsRA2;
      obsDec = obsDec2;
      sigastrom = sigastrom2;
      outdetvec = outdetvec2;
      // Also add this tracklet to the detections_added vector that will be output
      tracklets_added++;
      for(tct=0;tct<long(trkvec_temp.size());tct++) {
	detections_added.push_back(trkvec_temp[tct]);
      }
    }
  } else if(matchnum==1) {
    // Only one plausibly-fitting tracklet was found, so there is no
    // point in performing the pair-search or anything else.
    pairct = track_rmsvec[0].index;
    // Add this tracklet to the provisional fitting vectors.
    trkvec={};
    trkvec = tracklet_lookup(trk2det, pairct);
    status = add_tracklet01(observer_heliopos, obsMJD, obsRA, obsDec, sigastrom, outdetvec, trkvec, detvec, image_log, observer_heliopos2, obsMJD2, obsRA2, obsDec2, sigastrom2, outdetvec2, trkvec_temp);
    // Load the augmented vectors back into the base vectors for output
    observer_heliopos = observer_heliopos2;
    obsMJD = obsMJD2;
    obsRA = obsRA2;
    obsDec = obsDec2;
    sigastrom = sigastrom2;
    outdetvec = outdetvec2;
    // Also add this tracklet to the detections_added vector that will be output
    tracklets_added++;
    for(tct=0;tct<long(trkvec_temp.size());tct++) {
      detections_added.push_back(trkvec_temp[tct]);
    }
  } else if(matchnum<=0) {
    // No plausibly-fitting tracklets were found: exit without producing
    // an augmented file.
    cout << "No plausibly-fitting tracklets were found (best RMS was " << bestRMS << ")\n";
    cout << "Hence, no augmented file will be written\n";
    outstream1 << "No plausibly-fitting tracklets were found (best RMS was " << bestRMS << ")\n";
    outstream1 << "Hence, no augmented file will be written\n";
    outstream1.close();
    return(0);
  }
  // If we did not return in the matchnum<=0 case, the linkage
  // was successfully augmented to some extent. Perform a fit
  // to the final, augmented linkage.

  cout << "Attempting orbit fit to final linkage, augmented with " << tracklets_added << " tracklets totalling " << detections_added.size() << " additional points\n";
  outstream1 << "Attempting orbit fit to final linkage, augmented with " << tracklets_added << " tracklets totalling " << detections_added.size() << " additional points\n";
  arctrace03(polyorder,planetnum,planetmjd,planetmasses,planetpos,Sunpos,Sunvel,observer_heliopos,obsMJD,obsRA,obsDec,sigastrom,kepspan,minchichange, bestRA, bestDec, bestresid, outpos, outvel, &bestchi, &astromRMS, &refpoint, verbose);
  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromRMS << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  cout << "State vectors correspond to reference point " << refpoint << " in the input planet files\n";
  
  outstream1 << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromRMS << "\n";
  outstream1 << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  outstream1 << "State vectors correspond to reference point " << refpoint << " in the input planet files\n";

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
      outstream1 <<  outdetvec[i].image << "," << outdetvec[i].idstring << "," << outdetvec[i].band << ",";
      outstream1 << outdetvec[i].obscode << "," << outdetvec[i].known_obj << ",";
      outstream1 << outdetvec[i].det_qual << "," << outdetvec[i].index << "\n";
    }
  }
  outstream1.close();

  cout << tracklets_added << " tracklets, totalling " << detections_added.size() << " individual points,  were added to the orbit\n";

  return(0);
}

int add_tracklet01(const vector <point3LD> &observer_heliopos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <hldet> &outdetvec, const vector <long> &trkvec, const vector <hldet> &detvec, const vector <hlimage> &image_log, vector <point3LD> &observer_heliopos2, vector <long double> &obsMJD2, vector <double> &obsRA2, vector <double> &obsDec2, vector <double> &sigastrom2, vector <hldet> &outdetvec2, vector <long> &trkvec_temp)
{
  long obsnum,obsct,tct;
  point3LD tppos = point3LD(0,0,0);
  
  obsnum = obsMJD.size();
  if(obsnum != long(observer_heliopos.size()) || obsnum != long(obsRA.size()) || obsnum != long(obsDec.size()) || obsnum != long(sigastrom.size()) || obsnum != long(outdetvec.size())) {
    cerr << "ERROR: add_tracklet01 finds length mismatch in input observation vectors: " << obsnum << " " << observer_heliopos.size() << " " << obsRA.size() << " " << obsDec.size() << " " << sigastrom.size() << " " << outdetvec.size() << "\n";
    return(1);
  }
  obsct=tct=0;
  trkvec_temp={};
  observer_heliopos2 = {};
  obsMJD2 = {};
  obsRA2 = obsDec2 = sigastrom2 = {};
  outdetvec2 = {};
  while(obsct<long(obsMJD.size()) || tct<long(trkvec.size())) {
    if(obsct<long(obsMJD.size()) && (tct>=long(trkvec.size()) || obsMJD[obsct] < image_log[detvec[trkvec[tct]].image].MJD)) {
      // Next point in the time-ordered set is from the orginal observation vectors
      obsMJD2.push_back(obsMJD[obsct]);
      obsDec2.push_back(obsDec[obsct]);
      obsRA2.push_back(obsRA[obsct]);
      sigastrom2.push_back(sigastrom[obsct]);
      observer_heliopos2.push_back(observer_heliopos[obsct]);
      outdetvec2.push_back(outdetvec[obsct]);
      obsct++;
    } else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && (obsRA[obsct] != detvec[trkvec[tct]].RA || obsDec[obsct] != detvec[trkvec[tct]].Dec)) {
      // We have two points at the same time, but they are not identical
      // Keep the one from the original vectors
      obsMJD2.push_back(obsMJD[obsct]);
      obsDec2.push_back(obsDec[obsct]);
      obsRA2.push_back(obsRA[obsct]);
      sigastrom2.push_back(sigastrom[obsct]);
      observer_heliopos2.push_back(observer_heliopos[obsct]);
      outdetvec2.push_back(outdetvec[obsct]);
      obsct++;
      // And also keep the one from the new tracklet
      obsMJD2.push_back(image_log[detvec[trkvec[tct]].image].MJD);
      obsRA2.push_back(detvec[trkvec[tct]].RA);
      obsDec2.push_back(detvec[trkvec[tct]].Dec);
      sigastrom2.push_back(1.0);
      // Load input heliocentric observer position.
      tppos = point3LD(image_log[detvec[trkvec[tct]].image].X, image_log[detvec[trkvec[tct]].image].Y, image_log[detvec[trkvec[tct]].image].Z);
      observer_heliopos2.push_back(tppos);
      outdetvec2.push_back(detvec[trkvec[tct]]);
      trkvec_temp.push_back(trkvec[tct]);
      tct++;
    } else if(obsct<long(obsMJD.size()) && tct<long(trkvec.size()) && obsMJD[obsct] == image_log[detvec[trkvec[tct]].image].MJD && obsRA[obsct] == detvec[trkvec[tct]].RA && obsDec[obsct] == detvec[trkvec[tct]].Dec ) {
      // We have two exactly identical points. Keep the one from the original vectors...
      obsMJD2.push_back(obsMJD[obsct]);
      obsDec2.push_back(obsDec[obsct]);
      obsRA2.push_back(obsRA[obsct]);
      sigastrom2.push_back(sigastrom[obsct]);
      observer_heliopos2.push_back(observer_heliopos[obsct]);
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
      // Load input heliocentric observer position.
      tppos = point3LD(image_log[detvec[trkvec[tct]].image].X, image_log[detvec[trkvec[tct]].image].Y, image_log[detvec[trkvec[tct]].image].Z);
      observer_heliopos2.push_back(tppos);
      outdetvec2.push_back(detvec[trkvec[tct]]);
      trkvec_temp.push_back(trkvec[tct]);
      tct++;
    } else {
      cerr << "Error: logically excluded case in loading of new vectors\n";
      return(5);
    }
  }
  return(0);
}
    
