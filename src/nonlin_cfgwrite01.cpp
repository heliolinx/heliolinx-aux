/// April 24, 2024: nonlin_cfgwrite01.cpp:
// Combines the functionality of thorhyp01.cpp and wisetrack01.cpp,
// to write full parameter files for nonlinear digital tracking with
// Ari Heinze's hacked version of Joe Trollo's search.x code.

// Given observer coordinates in the barycentric ecliptic
// system, the corresponding MJD, a celestial position in RA, Dec,
// a distance from the observer in AU, and velocity sampling intervals
// in km/sec in counterclockwise circular speed relative to the Sun,
// the +z dimension of the tangential velocity relative to the Sun,
// and the radial velocity relative to the Sun -- plus the number of
// sample steps to probe in each of these three dimensions, calculate
// a state vector at a set earlier time that will pass
// through the exact specified RA, Dec position, as viewed from the
// observer's location, at the required target-observer distance and
// with the specified heliocentric velocity.

// Use these state vectors to predict the pixel positions corresponding
// to the hypothetical orbit on a sequence of images, and write the
// search.x cfg file based on these pixel positions.
// (spacecraft or otherwise), 
// 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3
#define NAMECOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MINOBSINTERVAL 1.0 // Minimum time-between-images in seconds
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME 1.5 // Default max inter-image time interval
                    // for tracklets, in hours (will be converted
                    // to days before use).
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define DEBUG 0
#define PHASE_G 0.15
#define IDNUMLEN 7

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
  cerr << "Usage: nonlin_cfgwrite01 -cfg configfile -images imfile -mjdref mjdref -mjdstart mjdstart -mjdend mjdend -projcen cenRA cenDec -platescale platescale -obsdist dist_from_observer(AU) -circvelsamp circvelsamp -zvelsamp zvelsamp -radvelsamp radvelsamp -circvelnum circvelnum -zvelnum zvelnum -radvelnum radvelnum -chunk1 output_file_chunk_1 -chunk2 output_file_chunk_2 -outroot outroot -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string imfile;
  string planetfile;
  string spaceobsfile;
  string outroot,outfile,chunk1,chunk2,cfgout;
  string cfgprefix = "nltrack_";
  string detlogprefix = "detlog_";
  string orbid,outlog,runid,detlog;
  long double mjdref = 0.0L;
  long double mjdstart = 0.0L;
  long double mjdend = 0.0L;
  long double circvel,circvelsamp,zvelsamp,radvelsamp;
  circvel = circvelsamp = zvelsamp = radvelsamp = 0L;
  long circvelnum,zvelnum,radvelnum;
  circvelnum = zvelnum = radvelnum = 0;
  int planetfile_startpoint=0;
  int planetfile_endpoint=0;
  int configread=0;
  int polyorder=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
  long double ldval=0.0L;
  long double xynorm,a,e,incl;
  xynorm=a=e=incl=0.0L;
  double dval=0.0L;
  vector <string> imagenames;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> obsmjd;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> targpos;
  vector <point3LD> targvel;
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  long double startX,startY,startZ,startVX,startVY,startVZ;
  startX = startY = startZ = startVX = startVY = startVZ = 0L;
  point3LD targ_to_obs = point3LD(0,0,0);
  long double obsdist=1.0L;
  long double sundist=1.0L;
  long double light_travel_time;
  vector <point3LD> observerpos;
  vector <long double> targMJD;
  vector <double> imageRA;
  vector <double> imageDec;
  vector <long double> imageMJD;
  point3LD imagepos = point3LD(0,0,0);
  point3LD imagevel = point3LD(0,0,0);
  int imnum=0;
  int imct=0;
  string lnfromfile;
  vector <string> linestringvec;
  int badread=0;
  int reachedeof=0;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  long double outRA,outDec;
  double racenter,deccenter,dist,pa,xpix,ypix,platescale;
  racenter = deccenter = dist = pa = xpix = ypix = 0.0;
  platescale = 1.0;
  int fieldnum=0;
  int NAMEcol = NAMECOL;
  int MJDcol = MJDCOL;
  int RAcol = RACOL;
  int Deccol = DECCOL;
  ofstream outstream1;
  ofstream outstream2;
  point3LD obsrefpos = point3LD(0,0,0);
  point3LD refpos = point3LD(0,0,0);
  point3LD refvel = point3LD(0,0,0);
  point3LD unitcirc = point3LD(0,0,0);
  point3LD unitrad = point3LD(0,0,0);
  point3LD unitz = point3LD(0,0,0);


  if(argc<19) {
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
	}
	// Read the ephemeris file for the space observatory (e.g., WISE)
	status=readconfigstring(instream1,spaceobsfile);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigstring(instream1,spaceobsfile);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Ephemeris file for space observatory is named " << spaceobsfile << "\n";
	obsmjd={};
	obspos={};
	obsvel={};
	read_horizons_fileLD(spaceobsfile,obsmjd,obspos,obsvel);
	cout << "Finished reading ephemeris file " << spaceobsfile << "\n";
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
    if(string(argv[i]) == "-i" || string(argv[i]) == "-im" || string(argv[i]) == "-imf" || string(argv[i]) == "-images" || string(argv[i]) == "-imfile" || string(argv[i]) == "--imagefile" || string(argv[i]) == "--imfile" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	// There is still something to read;
	imfile=argv[++i];
	i++;
      } else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdref" || string(argv[i]) == "-me" || string(argv[i]) == "-Mjdref" || string(argv[i]) == "--mjdref" || string(argv[i]) == "--Mjdref" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-projcen" || string(argv[i]) == "-pcen" || string(argv[i]) == "-projcenter" || string(argv[i]) == "--projcen" || string(argv[i]) == "--projection_center" || string(argv[i]) == "-projcenter" || string(argv[i]) == "-pc") {
      if(i+1 < argc) {
	//There is still something to read;
	racenter=stod(argv[++i]);
	if(i+1 >= argc) {
	  cerr << "projection center keyword supplied with fewer than the required 2 arguments\n";
	  show_usage();
	  return(1);
	}
	deccenter=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "projection center keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-platescale" || string(argv[i]) == "-pixscale" || string(argv[i]) == "--platescale" || string(argv[i]) == "--pixscale") {
      if(i+1 < argc) {
	//There is still something to read;
	platescale=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Plate scale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obsdist" || string(argv[i]) == "-observerdist" || string(argv[i]) == "--obsdist" || string(argv[i]) == "-geodist") {
      if(i+1 < argc) {
	//There is still something to read;
	obsdist=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Observer distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
     
    } else if(string(argv[i]) == "-circvelsamp") {
      if(i+1 < argc) {
	//There is still something to read;
	circvelsamp = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Circular velocity sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-zvelsamp") {
      if(i+1 < argc) {
	//There is still something to read;
	zvelsamp = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Z-direction velocity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-radvelsamp") {
      if(i+1 < argc) {
	//There is still something to read;
	radvelsamp = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Radial velocity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-circvelnum") {
      if(i+1 < argc) {
	//There is still something to read;
	circvelnum = stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Circular velocity step number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-zvelnum") {
      if(i+1 < argc) {
	//There is still something to read;
	zvelnum = stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Z-direction velocity step number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-radvelnum") {
      if(i+1 < argc) {
	//There is still something to read;
	radvelnum = stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Radial velocity step number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-chunk1" || string(argv[i]) == "--chunk1") {
      if(i+1 < argc) {
	//There is still something to read;
	chunk1=argv[++i];
	i++;
      }
      else {
	cerr << "File chunk1 keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-chunk2" || string(argv[i]) == "--chunk2") {
      if(i+1 < argc) {
	//There is still something to read;
	chunk2=argv[++i];
	i++;
      }
      else {
	cerr << "File chunk2 keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-out" || string(argv[i]) == "-outroot" || string(argv[i]) == "-o" || string(argv[i]) == "--outroot" || string(argv[i]) == "-root" || string(argv[i]) == "--root") {
      if(i+1 < argc) {
	//There is still something to read;
	outroot=argv[++i];
	i++;
      }
      else {
	cerr << "Output root keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--file" || string(argv[i]) == "-file" || string(argv[i]) == "-outputfile" || string(argv[i]) == "--outputfile") {
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
  cout << "input image file " << imfile << "\n";
  cout << "input starting MJD " << mjdstart << "\n";
  cout << "input ending MJD " << mjdend << "\n";
  cout << "output file " << outroot << "\n";

  // SET UP FOR LATER WISETRRACK-STYLE PORTION: Given an input
  // position and velocity at mjdstart, calculate the corresponding
  // pixel offset for every image in the image file.
  // Match mjdstart and mjdend to planet file.
  planetfile_startpoint = planetfile_endpoint = -99;
  for(j=0;j<long(planetmjd.size());j++) {
    if(fabs(planetmjd[j]-mjdstart) < IMAGETIMETOL/SOLARDAY) planetfile_startpoint = j;
    if(fabs(planetmjd[j]-mjdend) < IMAGETIMETOL/SOLARDAY) planetfile_endpoint = j;
  }
  if(planetfile_startpoint<0 || planetfile_endpoint<0) {
    cerr << "ERROR: mjdstart " << mjdstart << " and/or mjdend " << mjdend << " could not be matched\nto any timestep in the planet files\n";
    return(1);
  } else {
    cout << "mjdstart " << mjdstart << " corresponds to timestep " << planetfile_startpoint << " in the planet files\n";
    cout << "mjdend " << mjdend << " corresponds to timestep " << planetfile_endpoint << " in the planet files\n";
  }
  
  // Read input image file.
  instream1.open(imfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << imfile << "\n";
    return(1);
  }
  badread=0;
  while(reachedeof==0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(lnfromfile.size()>10 && (reachedeof==0 || reachedeof==1)) {
      fieldnum = get_col_vector01(lnfromfile, linestringvec);
      if(fieldnum>=NAMEcol) {
	imagenames.push_back(linestringvec[NAMEcol-1]);
      }
      if(fieldnum>=MJDcol) {
	ldval = stold(linestringvec[MJDcol-1]);
	imageMJD.push_back(ldval);
      } else badread=1;
      if(fieldnum>=RAcol) {
	dval = stod(linestringvec[RAcol-1]);
	imageRA.push_back(dval);
      } else badread=1;
      if(fieldnum>=Deccol) {
	dval = stod(linestringvec[Deccol-1]);
	imageDec.push_back(dval);
      } else badread=1;
    }
    cout << imageMJD.size() << " " << imageMJD[imageMJD.size()-1] << " " << imageRA[imageMJD.size()-1] << " " << imageDec[imageMJD.size()-1] << "\n";
  }
  imnum=imageMJD.size();
  if(imnum<2 || long(imagenames.size())!=imnum || long(imageRA.size())!=imnum || long(imageDec.size())!=imnum) {
    cerr  << fixed << setprecision(6) << "Error: image time vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << imnum << " " << imagenames.size() << " " << imageRA.size() << " " << imageDec.size() << "\n";
    return(1);
  }
  if(badread>=1){
    cerr << "Error reading " << imfile << ": apparent short line\n";
    return(1);
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();
  
  // Calculate the exact position of the observer at the time of each image.
  observerpos={};
  for(imct=0;imct<imnum;imct++) {
    //    observer_barycoords01LD(imageMJD[imct], 5, 0L, 0L, 0L, obsmjd, obspos, outpos);
    planetpos01LD(imageMJD[imct], polyorder, obsmjd, obspos, outpos);
    observerpos.push_back(outpos);
  }

  // Set up for thorhyp-style portion: a grid of hypothetical
  // orbits following specific sampling rules, all of which
  // pass through a specified observer-centric RA, Dec
  // position at mjdref.
  // Calculate observer's position at the reference time
  planetpos01LD(mjdref, polyorder, obsmjd, obspos, obsrefpos);
  
  cout << "mjdref = " << mjdref << ",Observer position: " << obsrefpos.x << " " << obsrefpos.y << " " << obsrefpos.z << "\n";
  // Calculate reference position for hypothetical asteroids
  celestial_to_statevecLD(racenter,deccenter,obsdist*AU_KM,refpos);
  // Add in observer's position to get absolute coordinates.
  refpos.x += obsrefpos.x;
  refpos.y += obsrefpos.y;
  refpos.z += obsrefpos.z;
  // Calculate reference velocity: a circular counterclockwise orbit with no Z-component.
  sundist = vecabs3LD(refpos);
  circvel = sqrt(GMSUN_KM3_SEC2/sundist);
  refvel.z = 0L;
  xynorm = sqrt(refpos.x*refpos.x + refpos.y*refpos.y);
  if(xynorm==0L) {
    // We are on the polar axis of the solar system: not responsible
    // for weird behavior. Just make it go full speed in the +X direction.
    refvel = point3LD(circvel,0L,0L);
  } else {
    // Circular counterclockwise orbit with no Z-component to the velocity at this instant.
    refvel = point3LD(-circvel*refpos.y/xynorm,circvel*refpos.x/xynorm,0L);
  }
  // Calculate orthogonal unit vectors for velocity system.
  if(xynorm==0L) {
    unitcirc = point3LD(1.0L,0L,0L);
    unitz = point3LD(0L,1.0L,0L);
    unitrad = point3LD(0L,0L,1.0L);
  } else {
    unitcirc = point3LD(-refpos.y/xynorm,refpos.x/xynorm,0L);
    unitz = point3LD(-refpos.x/xynorm*refpos.z/sundist,-refpos.y/xynorm*refpos.z/sundist,xynorm/sundist);
    unitrad = point3LD(refpos.x/sundist,refpos.y/sundist,refpos.z/sundist);
  }
  // Check orthogonality and normalization:
  cout << "Check 1 1 1 0 0 0: " << vecabs3LD(unitcirc) << " " << vecabs3LD(unitz) << " " << vecabs3LD(unitrad) << " " << dotprod3LD(unitcirc,unitz) << " " << dotprod3LD(unitcirc,unitrad) << " " << dotprod3LD(unitz,unitrad) << "\n";
  
  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "circvelnum zvelnum radvelnum refvelX refvelY refvelZ X Y Z VX VY VZ a e incl\n";
  // Loop on velocities
  for(long cvi = -circvelnum ; cvi <= circvelnum ; cvi++) {
    for(long zvi = -zvelnum ; zvi <= zvelnum ; zvi++) {
      for(long rvi = -radvelnum ; rvi <= radvelnum ; rvi++) {
	vector <long double> targMJDH = {};
	vector <point3LD> outposH = {};
	vector <point3LD> outvelH = {};
	point3LD targposH = point3LD(0,0,0);
	point3LD targvelH = point3LD(0,0,0);
	targvelH.x = refvel.x + double(cvi)*circvelsamp*unitcirc.x + double(zvi)*zvelsamp*unitz.x + double(rvi)*radvelsamp*unitrad.x;
	targvelH.y = refvel.y + double(cvi)*circvelsamp*unitcirc.y + double(zvi)*zvelsamp*unitz.y + double(rvi)*radvelsamp*unitrad.y;
	targvelH.z = refvel.z + double(cvi)*circvelsamp*unitcirc.z + double(zvi)*zvelsamp*unitz.z + double(rvi)*radvelsamp*unitrad.z;
	targposH = refpos;
	// Apply a light travel time correction to the reference position
	cout << "Time offset : " << (obsdist*AU/CLIGHT) << "\n";
	cout << "Targpos: " << targposH.x << " "  << targposH.y << " "  << targposH.z << "\n"; 
	targposH.x += targvelH.x*(obsdist*AU/CLIGHT);
	targposH.y += targvelH.y*(obsdist*AU/CLIGHT);
	targposH.z += targvelH.z*(obsdist*AU/CLIGHT); 
	cout << "Targpos: " << targposH.x << " "  << targposH.y << " "  << targposH.z << "\n";
	targMJDH.push_back(mjdstart);
	// Perform n-body integration beginning with the input state vector.
	// Note that here, the integration starts at mjdref, targposH, targvelH,
	// and proceeds (backwards in time) to targMJDH. The function integrate_orbit03LD()
	// is designed to calculate a vector of positions and velocities corresponding
	// to a vector of input MJD values -- hence, targMJDH, outposH, and outvelH are
	// vectors even though in this case each contains only one element.
	integrate_orbit03LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, targMJDH, targposH, targvelH, mjdref, outposH, outvelH);
	outstream1.precision(12);  
	outstream1 << cvi << " " << zvi << " " << rvi << " " << targvelH.x << " " << targvelH.y << " " << targvelH.z << " ";
	outstream1 << outposH[0].x << " " << outposH[0].y << " " << outposH[0].z << " ";
	outstream1 << outvelH[0].x << " " << outvelH[0].y << " " << outvelH[0].z << " ";

	long double v0 = vecabs3LD(targvelH);
	long double u = dotprod3LD(targvelH,targposH);
 	a = sundist*GMSUN_KM3_SEC2/(2.0l*GMSUN_KM3_SEC2-v0*v0*sundist);
	long double n = sqrt(GMSUN_KM3_SEC2/a/a/a);
	long double alpha = GMSUN_KM3_SEC2/a;
	if(alpha>0.0l) {
	  // Bound, elliptical orbit
	  long double EC = 1.0L - sundist/a;
	  long double ES = u/n/a/a;
	  e = sqrt(EC*EC + ES*ES);
	} else if (alpha<0.0l) {
	  // Unbound, hyperbolic orbit
	  long double CH = 1.0L - sundist/a;
	  long double SH = u/sqrt(-GMSUN_KM3_SEC2*a);
	  e = sqrt(CH*CH - SH*SH);
	}
	point3LD lvec = crossprod3LD(targposH,targvelH);
	vecnorm3LD(lvec);
	incl = DEGPRAD*acos(lvec.z);
	if(lvec.z > 1.0L) incl = 90.0L;
	else if(lvec.z < -1.0L) incl = -90.0L;
	outstream1.precision(6);  
	outstream1 << a/AU_KM << " " << e << " " << incl << "\n";
	// This concludes the part taken from thorhyp01.cpp. Now we
	// move to a part taken from wisetrack01.cpp. For this portion,
	// integration starts at mjdstart,outpos,outvel, and the position
	// is calculated at the time of each image.

	outpos = outposH[0];
	outvel = outvelH[0];
	// Perform n-body integration beginning with the input state vectors.
	integrate_orbit04LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, outpos, outvel, planetfile_startpoint, planetfile_endpoint, targMJD, targpos, targvel);
	// Construct the name of the output file
	orbid="";
	if(cvi>=0) {
	  orbid.push_back('p');
	  string s = to_string(cvi);
	  orbid += s;
	} else {
	  orbid.push_back('m');
	  string s = to_string(-cvi);
	  orbid += s;
	}
	if(zvi>=0) {
	  orbid.push_back('p');
	  string s = to_string(zvi);
	  orbid += s;
	} else {
	  orbid.push_back('m');
	  string s = to_string(-zvi);
	  orbid += s;
	}
	if(rvi>=0) {
	  orbid.push_back('p');
	  string s = to_string(rvi);
	  orbid += s;
	} else {
	  orbid.push_back('m');
	  string s = to_string(-rvi);
	  orbid += s;
	}
	cfgout = cfgprefix + outroot + orbid + ".cfg";
	cout << "Writing output cfg file called " << cfgout << "\n";
	// Construct the names for the output log, runID, and detectionlog
	outlog = outroot + orbid + ".log";
	runid = outroot + orbid;
	detlog = detlogprefix + outroot + orbid + ".txt";
	outstream2.open(cfgout);
	// Write first lines involving opLog and runID
	outstream2 << "opLog = \"" << outlog << "\";\n";
	outstream2 << "runID = \"" << runid << "\";\n";
	// concatenate on file chunk1
	instream1.open(chunk1);
	while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
	  lnfromfile = "";
	  getline(instream1,lnfromfile);
	  outstream2 << lnfromfile << "\n";
	}
	instream1.close();
	// Write out the part with detlog and source size.
	outstream2 << "detLog = \"" << detlog << "\";\n";
	outstream2 << "outputEvery = 0;\n";
	outstream2 << "stopAfter = 0;\n";
	outstream2 << "sourceSizeLog = \"sourceSize" << runid << ".txt\";\n";
	outstream2 << "stackSizeLog = \"stackSize" << runid << ".txt\";\n";
	// concatenate on file chunk2
	instream1.open(chunk2);
	while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
	  lnfromfile = "";
	  getline(instream1,lnfromfile);
	  outstream2 << lnfromfile << "\n";
	}
	instream1.close();
	// Write start of image section
	outstream2 << "images = (\n";
	for(imct=0;imct<imnum;imct++)
	  {
	    planetposvel01LD(imageMJD[imct],polyorder,targMJD,targpos,targvel,imagepos,imagevel);
	    startX = imagepos.x;
	    startY = imagepos.y;
	    startZ = imagepos.z;
	    targ_to_obs.x = observerpos[imct].x - imagepos.x;
	    targ_to_obs.y = observerpos[imct].y - imagepos.y;
	    targ_to_obs.z = observerpos[imct].z - imagepos.z;
	    long double imobsdist = sqrt(targ_to_obs.x*targ_to_obs.x + targ_to_obs.y*targ_to_obs.y + targ_to_obs.z*targ_to_obs.z);
	    // Instantaneous velocity
	    planetpos01LD(imageMJD[imct],polyorder,targMJD,targvel,imagevel);
	    // Initial approximation of the coordinates relative to the observer
	    imagepos.x -= observerpos[imct].x;
	    imagepos.y -= observerpos[imct].y;
	    imagepos.z -= observerpos[imct].z;
	    light_travel_time = imobsdist*1000.0/CLIGHT; // Factor of 1000 converts imobsdist to meters
	    // Light-travel-time corrected version of coordinates relative to the observer
	    imagepos.x -= light_travel_time*imagevel.x;
	    imagepos.y -= light_travel_time*imagevel.y;
	    imagepos.z -= light_travel_time*imagevel.z;
	    // Light-travel-time corrected observer-target distance
	    imobsdist = sqrt(imagepos.x*imagepos.x + imagepos.y*imagepos.y + imagepos.z*imagepos.z);
	    // Calculate unit vector
	    imagepos.x /= imobsdist;
	    imagepos.y /= imobsdist;
	    imagepos.z /= imobsdist;
	    // Project onto the celestial sphere.
	    stateunitLD_to_celestial(imagepos, outRA, outDec);
	    // Carry out stereographic projection
	    distradec02(racenter,deccenter,outRA,outDec,&dist,&pa);
	    dist = 2.0l*tan(dist/DEGPRAD/2.0l);
	    xpix = -dist*sin(pa/DEGPRAD)*ASECPRAD/platescale; // East is left
	    ypix = dist*cos(pa/DEGPRAD)*ASECPRAD/platescale;
	    
	    cout << fixed << setprecision(6) << imageMJD[imct] << " " << outRA << " " << outDec << " ";
	    cout << fixed << setprecision(3) << xpix << " " << ypix << " ";
	    cout << fixed << setprecision(6) << observerpos[imct].x << " " << observerpos[imct].y << " " << observerpos[imct].z << " ";
	    cout << fixed << setprecision(6) << startX << " " << startY << " " << startZ << " ";
	    cout << fixed << setprecision(6) << imagevel.x << " " << imagevel.y << " " << imagevel.z << "\n";
	    // Write image data to output cfgfile
	    outstream2 << fixed << setprecision(6) << "{ name = \"" << imagenames[imct] << "\" time = " << 24.0L*(imageMJD[imct]-mjdref);
	    outstream2 << fixed << setprecision(3) << " xcen = " << xpix << " ycen = " << ypix;
	    if(imct<imnum-1) outstream2 << "},\n";
	    else outstream2 << "}\n";
	  }
	outstream2 << ");\n";
	outstream2.close();
      }
    }
  }
  outstream1.close();
  
  cout << "mjdref = " << mjdref << ",Observer position: " << obsrefpos.x << " " << obsrefpos.y << " " << obsrefpos.z << "\n";
  for(imct=0;imct<imnum;imct++) {
    cout << imagenames[imct] << " " << imageMJD[imct] << " " << observerpos[imct].x << " " << observerpos[imct].y << " " << observerpos[imct].z << "\n";
  }
  return(0);
}
