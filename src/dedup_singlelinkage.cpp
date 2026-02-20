// February 03, 2026: dedup_singlelinkage.cpp
// Attempt to remove duplicates from a single linkage,
// performing an orbit-fit if necessary.
// The input data must already be time-sorted, at least
// to matchtol accuracy.

#define DAY_TO_DEG_CONV 24.0 // Treats one day as equivalent to 24 degrees:
                             // in other words, one second of time is
                             // equivalent to one arcsecond on the sky.
#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define OBSCODECOL 5
#define COLS_TO_READ 5
#define HNUM 4

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: dedup_singlelinkage.cpp -cfg configfile -infile infile -colformat column_format_file -matchtol matchtol_in_seconds -kepspan time_span_for_Keplerian_fit(day) -minchi min_chi_change -rmsthresh astrometric_rms_threshold -obscode obscodefile -maxiter maxiter -outfile outfile -isTAI isTAI -verbose verbosity\n";
}

int main(int argc, char *argv[])
{
  double kepspan=15.0;
  long double minchichange = 0.001l;
  double timestep = 5.0;
  long maxiter = 10;
  ofstream outstream1;
  ifstream instream1;
  string configfile,colformatfile,infile,stest,idstring;
  vector <string> ObsCodevec;
  string obscodefile;
  vector <observatory> observatory_list = {};
  string obscodestring,Sunfile,Earthfile,planetfile;
  char obscodechar[MINSTRINGLEN];
  string detfile, uniquefile,dupfile,lnfromfile,outfile;
  double ldval;
  vector <double> MJDvec;
  vector <double> TDBvec;
  vector <double> RAvec;
  vector <double> Decvec;
  vector <double> mjdtest;
  vector <vector <double>> statevecs;
  vector <double> starting_statevec;
  vector <double> out_statevec;
  vector <double> Sun_mjd;
  vector <vector <double>> Sun_statevecs;
  vector <double> Earth_mjd;
  vector <vector <double>> Earth_statevecs;
  vector <double> planetmasses;
  vector <double> planet_backward_mjd;
  vector <double> planet_forward_mjd;
  vector <vector <double>> planet_backward_statevecs;
  vector <vector <double>> planet_forward_statevecs;
  int configread=0;
  int planetnum,planetct;
  long detnum;
  vector <string> linevec;
  vector <string> IDvec;
  double MJD,RA,Dec;
  double timescale = DAY_TO_DEG_CONV;
  double matchtol = 1.0;
  long i,j,k,bandlen,status;
  i = j = k = bandlen = status = 0;
  long colreadct,lct;
  int reachedeof = 0;
  int startpoint,endpoint;
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int obscodecol = OBSCODECOL;
  int colformatfile_set=0;
  vector <string> linklines;
  vector <long> linkdets;
  long isdup=0;
  vector <long> dupivec;
  vector <vector <long>> dup_indexvec;
  long detct=0;
  long dupnum,goodnum;
  
  int hnum = HNUM;
  vector <double> hspace;
  int verbose = 0;
  vector <double> obsstate;
  vector <double> sunstate;
  vector <vector <double>> observer_statevecs;
  vector <vector <double>> observer_heliostate;
  long refpoint;
  double obslon,plxcos,plxsin,dist;
  double astromrms=0.0;
  double astromrmsthresh = 0.1;
  double mjdref=0.0;

  vector <vector <double>> Kepobserverpos;
  vector <double> KepMJD;
  vector <double> KepRA;
  vector <double> KepDec;
  vector <double> Kepsig;
  double geodist1 = 1.0;
  double geodist2 = 1.1;
  double ftol = FTOL_HERGET_SIMPLEX;
  double simplex_scale = SIMPLEX_SCALEFAC;
  int bestpoint,point1;
  int simptype,point2,kepnum,kepmax;
  bestpoint=point1=0;
  simptype=point2=kepnum=kepmax=1;
  double kepmetric,kepmetbest;
  kepmetric=kepmetbest=1.0;
  double stateMJD = 0.0;
  long itnum = 0;
  double chisq=0.0l;
  vector <long> keepvec;
  int isTAI=0;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> fitresid;
  vector <double> sigastrom;
  long bp;
  double bestresid=LARGERR2;

  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    //cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
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
	} //else cout << "HNUM (number of substeps per big step) for Everhart integrator read as " << hnum << "\n";
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
	} //else cout << "Timestep for Everhart integrator (big step in days) read as " << timestep << "\n";
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
	} //else cout << "Evenly-sampled ephemeris file for the sun is named " << Sunfile << "\n";
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
	} //else cout << "Evenly-sampled ephemeris file for Earth is named " << Earthfile << "\n";
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
	} else if(verbose>0) cout << "Number of planets read as " << planetnum << "\n";
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
	    if(verbose>0) cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
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
	  } else if(verbose>0) cout << "Backward ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
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
		cerr << "ERROR: backward time vectors do not match for input planet files\n";
		cerr << planetct+1 << " and 1!\n";
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
	  } else if(verbose>0) cout << "Forward ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
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
		cerr << "ERROR: forward time vectors do not match for input planet files\n";
		cerr << planetct+1 << " and 1!\n";
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
    } else if(verbose>0) cout << "Configuration file read successfully\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-in" || string(argv[i]) == "-filelist" || string(argv[i]) == "-file_list" || string(argv[i]) == "--infile" || string(argv[i]) == "--file_list") {
      if(i+1 < argc) {
	//There is still something to read;
	infile=argv[++i];
	i++;
      }
      else {
	cerr << "Input file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format"  || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile=argv[++i];
	colformatfile_set = 1;
	i++;
      }
      else {
	cerr << "Column format file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-matchtol" || string(argv[i]) == "-match_tol" || string(argv[i]) == "-mtol" || string(argv[i]) == "-matching_tolerance" || string(argv[i]) == "--matching_tolerance" || string(argv[i]) == "-tolerance") {
      if(i+1 < argc) {
	//There is still something to read;
	matchtol=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Matching tolerance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timescale" || string(argv[i]) == "-time_scale" || string(argv[i]) == "-tscale" || string(argv[i]) == "-day2deg" || string(argv[i]) == "-day2degscale" || string(argv[i]) == "-day2deg_scale"|| string(argv[i]) == "--timescale" || string(argv[i]) == "--time_scale" || string(argv[i]) == "--day2deg" || string(argv[i]) == "--day2degscale" || string(argv[i]) == "--day2deg_scale") {
      if(i+1 < argc) {
	//There is still something to read;
	timescale=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Timescale keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-rmsthresh") {
      if(i+1 < argc) {
	//There is still something to read;
	astromrmsthresh=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Astrometric RMS threshold keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-maxiter") {
      if(i+1 < argc) {
	//There is still something to read;
	maxiter=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "max iterations keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--output_file") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output unique file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-isTAI" || string(argv[i]) == "-TAI" || string(argv[i]) == "-istai" || string(argv[i]) == "-TAItime" || string(argv[i]) == "-taitime" || string(argv[i]) == "--isTAI" || string(argv[i]) == "--input_TAI") {
      if(i+1 < argc) {
	//There is still something to read;
	isTAI=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "keyword for whether input tims is TAI supplied with no corresponding argument";
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
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }
  cout.precision(17);  

  if(verbose>0) cout << "Input file is called " << infile << "\n";
  if(colformatfile_set==0) {
    cout << "WARNING: No column formatting file supplied\n";
  } else if(verbose>0) cout << "Column formatting file for is " << colformatfile << "\n";
  if(verbose>0) cout << "Matching tolerance will be " << matchtol << " seconds\n";
  if(verbose>0) cout << "Timescale for converting 1 day of time to equivalent degrees will be " << timescale << "\n";
  if(verbose>0) cout << "Output unique file will be called " << uniquefile << "\n";
  if(verbose>0) cout << "Output duplicate file will be called " << dupfile << "\n";
  
  // Read the column formatting file
  if(colformatfile.size()>0)
    {
      instream1.open(colformatfile);
      if(!instream1)  {
	cerr << "ERROR: unable to open input file " << colformatfile << "\n";
	return(1);
      }
      colreadct=0;
      while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct<COLS_TO_READ) {
	instream1 >> stest;
	if(stest == "MJDCOL") {
	  instream1 >> mjdcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "RACOL") {
	  instream1 >> racol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DECCOL") {
	  instream1 >> deccol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if (stest == "IDCOL") {
	  instream1 >> idcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if (stest == "OBSCODECOL") {
	  instream1 >> obscodecol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} 
      }
      instream1.close();
      if(colreadct<COLS_TO_READ) {
	cout << "WARNING: only " << colreadct << " column specifications, of " << COLS_TO_READ << " expected, were read from column format file " << colformatfile << ".\n";
      }
    }
  
  // Read input file
  instream1.open(infile);
  if(!instream1) {
    cerr << "can't open input detection file " << infile << "\n";
    return(1);
  }
  if(verbose>0) cout << "Successfully opened file " << infile << "\n";
  reachedeof = lct = 0;
  // We assume there is no header
  MJDvec = Decvec = RAvec = sigastrom = {};
  linevec = IDvec = {};
  ObsCodevec = {};
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(verbose>1) cout << "Just read line " << lnfromfile << "\n";
    if(verbose>1) cout << "reachedeof = " << reachedeof << " " << instream1.eof() << " " << instream1.fail() << " " << instream1.bad() << " " << lct << "\n";
	
    MJD = 0.0l;
    RA = Dec = -999.9;
    idstring = "";
    startpoint = endpoint = j = 0;
    while(endpoint<long(lnfromfile.size()) && startpoint<long(lnfromfile.size()) && lnfromfile.size()>=15 && reachedeof == 0) {
      endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      j++;
      if(j==idcol) {
	if(endpoint>0) {
	  idstring = stest;
	}
      } else if(j==mjdcol) {
	if(endpoint>0) {
	  try { MJD = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read MJD string " << stest << " from line " << lnfromfile << "\n"; }
	}
      } else if(j==racol) {
	if(endpoint>0) {
	  try { RA = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read RA string " << stest << " from line " << lnfromfile << "\n"; }
	}
      } else if(j==deccol) {
	if(endpoint>0) {
	  try { Dec = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read DEC string " << stest << " from line " << lnfromfile << "\n"; }
	}
      } else if(j==obscodecol) {
	if(endpoint>0) {
	  obscodestring = stest;
	}
      }
      startpoint = endpoint+1;
    }
    if(reachedeof == 0 && lnfromfile.size()>=15) {
      if(MJD==0.0) {
	cerr << "ERROR: MJD not read from line " << lct << " of input file " << detfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(RA==-999.9) {
	cerr << "ERROR: RA not read from line " << lct << " of input file " << detfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(Dec==-999.9) {
	cerr << "ERROR: Dec not read from line " << lct << " of input file " << detfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(idstring.size()<=0) {
	cerr << "ERROR: Object ID not read from line " << lct << " of input file " << detfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(obscodestring.size()<=0 || obscodestring.size()>=MINSTRINGLEN) {
	cerr << "ERROR: ObsCode string " << obscodestring << " of length " << obscodestring.size() << " is invalid\n";
	return(2);
      }
      MJDvec.push_back(MJD);
      RAvec.push_back(RA);
      Decvec.push_back(Dec);
      IDvec.push_back(idstring);
      ObsCodevec.push_back(obscodestring);
      linevec.push_back(lnfromfile);
      sigastrom.push_back(astromrmsthresh);
    }
  }
  instream1.close();
  detnum = MJDvec.size();
  if(detnum!=long(RAvec.size()) || detnum!=long(Decvec.size()) || detnum!=long(IDvec.size()) || detnum!=long(linevec.size())) {
    cerr << "ERROR: vector mismatch after reading file " << detfile << " : " << MJDvec.size() << " " << RAvec.size() << " " << Decvec.size() << " " << IDvec.size() << " " << linevec.size() << "\n";
    return(2);
  }
  if(verbose>0) cout << "Successfully read " << detnum << " lines from input file " << infile << "\n";

  // Load keepvec, and the duplicate index vectors.
  dup_indexvec={};
  keepvec={};
  for(detct=0;detct<detnum;detct++) {
    dupivec={};
    dup_indexvec.push_back(dupivec);
    keepvec.push_back(1);
  }
  // See if there are any time-duplicates
  isdup = 0;
  for(detct=0;detct<detnum;detct++) {
    if(keepvec[detct]>=1) {
      for(j=detct+1;j<detnum;j++) {
	if(fabs(MJDvec[j] - MJDvec[detct])<matchtol/SOLARDAY) {
	  isdup=1;
	  // We have a duplicate. Add it to dup_indexvec
	  dup_indexvec[detct].push_back(j);
	  // Mark it as a duplicate in keepvec
	  keepvec[j]=0;
	}
      }
    }
  }
  if(isdup==0) {
    // There were no duplicates. Simply copy all the input lines to output
    cout << "No duplicates found, writing " << detnum << " lines to output file\n";
    outstream1.open(outfile);
    for(detct=0;detct<detnum;detct++) {
      outstream1 << linevec[detct] << "\n";
    }
    outstream1.close();
  } else {
    // PERFORM ASTROMETRIC FIT
    cout.precision(17);
    if(verbose>0) cout << "input configuration file " << configfile << "\n";
    if(verbose>0) cout << "input observatory code file " << obscodefile << "\n";
    if(verbose>0) cout << "output file " << outfile << "\n";
    // Load a bunch of preliminary stuff

    // Load hspace vector based on selected value of hnum
    // Note that all of these numerical values are from solarsyst_dyn_geo01.h

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
    // Read observatory code file
    status = read_obscode_file2(obscodefile, observatory_list, verbose);
    if(status!=0) {
      cerr << "ERROR reading observatory code file " << obscodefile << "\n";
      return(1);
    }
    if(verbose>0) cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";
    if(verbose>=2) {
      for(i=0;i<long(observatory_list.size());i++) {
	cout << observatory_list[i].obscode << " " << observatory_list[i].obslon << " " << observatory_list[i].plxcos << " " << observatory_list[i].plxsin << "\n";
      }
    }
    // Load TDB time
    TDBvec={};
    if(isTAI==1) {
      // Assume input MJD values are TAI
      for(detct=0;detct<detnum;detct++) TDBvec.push_back(MJDvec[detct] -TAI_MINUS_UTC/SOLARDAY + TTDELTAT/SOLARDAY);
    } else {
      // Assume input MJD values are UTC
      for(detct=0;detct<detnum;detct++) TDBvec.push_back(MJDvec[detct] + TTDELTAT/SOLARDAY);
    }

    // Calculate the observer's heliocentric position at the time of each observation.
    observer_statevecs = observer_heliostate = {};
    for(detct=0;detct<detnum;detct++) {
      strcpy(obscodechar,ObsCodevec[detct].c_str());
      status = obscode_lookup(observatory_list,obscodechar,obslon,plxcos,plxsin);
      if(verbose>0) cout << "Detection " << detct << " is from obscode " << obscodechar << ", coords " << obslon << " " << plxcos << " " << plxsin << "\n";
      if(status>0) {
	cerr << "ERROR: obscode_lookup failed for observatory code " << obscodechar << "\n";
	return(3);
      }
      // Calculate observer's exact barycentric position and velocity.
      // Note that observer_barystate01 assumes input is UTC, so we apply a correction to TDBvec
      if(verbose>0) cout << "About to call observer_barystate01 " << Earth_mjd.size() << " " << Earth_statevecs.size() << "\n";
      observer_barystate01(TDBvec[detct]-TTDELTAT/SOLARDAY, 5, obslon, plxcos, plxsin, Earth_mjd, Earth_statevecs, obsstate, verbose);
      observer_statevecs.push_back(obsstate);
      // Calculate the sun's position at the same time.
      // Note that planetpos02 assumes input is TDB
      planetpos02(TDBvec[detct],5,Sun_mjd,Sun_statevecs,sunstate);
      // Convert obsstate from barycentric to heliocentric coords.
      for(k=0;k<6;k++) obsstate[k] -= sunstate[k];
      observer_heliostate.push_back(obsstate);
      if(verbose>0) {
	cout << "Observation " << detct << ", barystate : " << observer_statevecs[detct][0] << " " << observer_statevecs[detct][1] << " " << observer_statevecs[detct][2] << " " << observer_statevecs[detct][3] << " " << observer_statevecs[detct][4] << " " << observer_statevecs[detct][5] << "\n";
	cout << "heliostate : " << observer_heliostate[detct][0] << " " << observer_heliostate[detct][1] << " " << observer_heliostate[detct][2] << " " << observer_heliostate[detct][3] << " " << observer_heliostate[detct][4] << " " << observer_heliostate[detct][5] << "\n";
      }
    }
    if(verbose>0) cout << MJDvec.size() << " " << MJDvec[MJDvec.size()-1] << " " << RAvec[MJDvec.size()-1] << " " << Decvec[MJDvec.size()-1] << "\n";
    if(detnum<2 || long(RAvec.size())!=detnum || long(Decvec.size())!=detnum || long(observer_heliostate.size())!=detnum) {
      cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
      cerr  << fixed << setprecision(6) << detnum << " " << RAvec.size() << " " << Decvec.size() << " " << observer_heliostate.size() << "\n";
      return(1);
    }
  
    // SET UP FOR PRELIMINARY KEPLERIAN FIT
    // Find the time period of length kepspan
    // with the largest number of observations
    bestpoint=kepnum=kepmax=0;
    kepmetric=kepmetbest=0.0;
    for(i=0;i<detnum;i++) {
      j=i;
      kepnum=0;
      kepmetric=0.0;
      while(j<detnum && MJDvec[j]-MJDvec[i] < kepspan) {
	kepnum++;
	j++;
      }
      kepmetric = double(kepnum)*(MJDvec[j-1]-MJDvec[i]);
      if(kepmetric>kepmetbest) {
	kepmetbest=kepmetric;
	kepmax=kepnum;
	bestpoint=i;
      }
    }
    kepnum=kepmax;
    if(verbose>0) cout << "Best sequence contains " << kepnum << " points, and runs from MJDvec[" << bestpoint << "]=" << MJDvec[bestpoint] << " to MJDvec[" << bestpoint+kepnum-1 << "]=" << MJDvec[bestpoint+kepnum-1] << "\n";
    KepMJD = KepRA = KepDec = Kepsig = {};
    Kepobserverpos = {};
    for(i=bestpoint;i<bestpoint+kepnum;i++) {
      KepMJD.push_back(MJDvec[i]);
      Kepobserverpos.push_back(observer_statevecs[i]);
      KepRA.push_back(RAvec[i]);
      KepDec.push_back(Decvec[i]);
      Kepsig.push_back(sigastrom[i]);
    }
    if(verbose>0) cout << "kepnum = " << kepnum << " " << KepMJD.size() << "\n";
    for(i=0;i<kepnum;i++) {
      if(verbose>0) cout << Kepobserverpos[i][1] << " " << Kepobserverpos[i][2] << " " << Kepobserverpos[i][3] << " " << KepMJD[i] << " " << KepRA[i] << " " << KepDec[i] << "\n";
    }
    fitDec = fitRA = fitresid = out_statevec = {};
    if(verbose>0) cout << "Launching Hergetfit_vstarSV\n";
    chisq = Hergetfit_vstarSV(geodist1, geodist2, simplex_scale, simptype, ftol, 1, kepnum, Kepobserverpos, KepMJD, KepRA, KepDec, Kepsig, fitRA, fitDec, fitresid, out_statevec, stateMJD, itnum, verbose);
    if(verbose>0) cout << "Keplerian fit produced chisq = " << chisq << "\n";

    // Choose a suitable reference MJD near the middle of the Keplerian fit.
    // Require that it also fall exactly on an integer time step for the Everhart
    // sampling read from the configuration file.
    mjdref = 0.5*(KepMJD[0]+KepMJD[kepnum-1]);
    mjdref = round(mjdref/timestep)*timestep;
    // Make sure that mjdref matches a point in the Everhart-sampled planet files
    refpoint = -99;
    for(j=0;j<long(planet_forward_mjd.size());j++) {
      if(fabs(planet_forward_mjd[j]-mjdref) < STATEMJD_TIMETOL) refpoint = j;
    }
    if(refpoint<0) {
      cerr << "ERROR: reference mjd " << mjdref << " from the Keplerian fit did not match any point in the input Everhart-sampled ephemeris vectors\n";
    }
    if(verbose>0) cout << "Reference MJD " << mjdref << " corresponds to point " << refpoint << " in the Everhart-sampled ephemeris vectors\n";

    // Integrate the Keplerian state vectors from their current reference MJD
    // to the MJD corresponding to planetfile_refpoint
    Kepler_univ_int_SV(GMSUN_KM3_SEC2, stateMJD, out_statevec, mjdref, starting_statevec, verbose);
    if(verbose>0) cout << "Kepler_univ_int_SV output statevec for MJD " << mjdref << ":\n";
    if(verbose>0) cout << starting_statevec[0] << " " << starting_statevec[1] << " " << starting_statevec[2] << " " << starting_statevec[3] << " " << starting_statevec[4] << " " << starting_statevec[5] << "\n";

    out_statevec={};
    fitDec = fitRA = fitresid = out_statevec = {};
    status = evertrace01(planetnum, planetmasses, planet_backward_mjd, planet_backward_statevecs, planet_forward_mjd, planet_forward_statevecs, starting_statevec, mjdref, MJDvec, observer_statevecs, RAvec, Decvec, sigastrom, fitRA, fitDec, out_statevec, timestep, hnum, hspace, minchichange, astromrmsthresh, maxiter, itnum, astromrms, chisq, verbose);
    if(status!=0) {
      cerr << "ERROR: evertrace01 returned error status " << status << "\n";
      return(status);
    }
    for(detct=0;detct<detnum;detct++) {
      dist = 3600.0*distradec01(RAvec[detct],Decvec[detct],fitRA[detct],fitDec[detct]);
      fitresid.push_back(dist);
    }
    // Now go through the detections and for every case of overlap,
    // ouput the one with the lowest residual.
    outstream1.open(outfile);
    dupnum=goodnum=0;
    for(detct=0;detct<detnum;detct++) {
      if(keepvec[detct]>=1) {
	// We need to output some line corresponding to this detection.
	if(dup_indexvec[detct].size()<=0) {
	  // There were no duplicates: just write the line as-is.
	  outstream1 << linevec[detct] << "\n";
	  goodnum++;
	} else {
	  // There are some duplicates. Find the one with the lowest astromeric residual.
	  bp = detct;
	  bestresid = fitresid[detct];
	  dupnum += dup_indexvec[detct].size();
	  for(j=0;j<long(dup_indexvec[detct].size());j++) {
	    if(fitresid[dup_indexvec[detct][j]] < bestresid) {
	      bp = dup_indexvec[detct][j];
	      bestresid = fitresid[dup_indexvec[detct][j]];
	    }
	  }
	  // Output the line corresponding to the best residual.
	  outstream1 << linevec[bp] << "\n";
	  goodnum++;
	}
      }
    }
    outstream1.close();
    cout << "Found " << dupnum << " duplicates out of " << detnum << " total detections; writing " << goodnum << " lines to output file " << outfile << "\n";

  }
  return(0);
}
