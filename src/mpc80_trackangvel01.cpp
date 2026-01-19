// January 13, 2026: mpc80_trackangvel01.cpp:
// Given an input file in mpc80 format, identify distinct tracklets
// and calculate their angular velocities in ecliptic coordinates
// (ecliptic longitude and latitude), and also calculate their
// sun-relative ecliptic coordinates

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: mpc80_trackangvel -mpcfile mpcfile -minpts minpts -minarc minarc -mintime mintime -maxGCR maxGCR -outfile outfile -verbose verbosity\n";
}
    
int main(int argc, char *argv[])
{
  string mpcfile, outfile, lnfromfile;
  int minpts = 3;
  double minarc = 5.0; // Min tracklet arc in arcseconds
  double mintime = 0.02; // Min tracklet timespan in days
  double maxGCR = 1.0;
  string objname,oldname,obscode,oldcode,band;
  double MJD,RA,Dec,mag,oldmjd;
  vector <double> MJDvec;
  vector <double> RAvec;
  vector <double> Decvec;
  double poleRA,poleDec,angvel,pa,crosstrack,alongtrack,arc,timespan,GCR;
  double ecliplon,ecliplat,lonmed,lonmean,lonrms,latmean,MJDmean;
  int status=0;
  double sunra,sundec,sunlon,sunlat,ecliplon_relative_to_sun,ecliplon_relative_to_opposition;
  long i=0;
  int verbose=0;
  
  
  if(argc<5) {
    show_usage();
    return(1);
  }

  // Parse input values
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-mpcfile" || string(argv[i]) == "-mpc" || string(argv[i]) == "-MPC" || string(argv[i]) == "-mpc80" || string(argv[i]) == "--mpcfile" || string(argv[i]) == "--mpc80" || string(argv[i]) == "--MPC80") {
      if(i+1 < argc) {
	//There is still something to read;
	mpcfile=argv[++i];
	i++;
      } else {
	cerr << "mpcfile keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minpts" || string(argv[i]) == "-minpoints" || string(argv[i]) == "-min_points") {
      if(i+1 < argc) {
	//There is still something to read;
	minpts=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum point number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minarc" || string(argv[i]) == "-min_arc" || string(argv[i]) == "-min_arclength") {
      if(i+1 < argc) {
	//There is still something to read;
	minarc=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum arc length keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mintime" || string(argv[i]) == "-min_time" || string(argv[i]) == "-min_timespan") {
      if(i+1 < argc) {
	//There is still something to read;
	mintime=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum time span keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxGCR" || string(argv[i]) == "-maxgcr" || string(argv[i]) == "-max_GCR") {
      if(i+1 < argc) {
	//There is still something to read;
	maxGCR=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Maximum GCR keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-of" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "-out" || string(argv[i]) == "-ofile" || string(argv[i]) == "-outf" || string(argv[i]) == "--OutputFile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile = argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "verbose keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Unrecognized input keyword " << argv[i] << "\n";
      show_usage();
      return(1);
    }
  }

  // Echo input
  cout << "Input mpc file: " << mpcfile << "\n";
  cout << "Minimum points per tracklet: " << minpts << "\n";
  cout << "Minimum tracklet arc length: " << minarc << " arcsec\n";
  cout << "Minimum tracklet time span: " << mintime << " days = " << mintime*1440.0 << " minutes\n";
  cout << "Maximum GCR: " << maxGCR << " arcsec\n";
  cout << "Output file name: " << outfile << "\n";
  
  // Read input detection file
  ifstream instream1 {mpcfile};
  ofstream outstream1 {outfile};
  outstream1 << "ecliplon_relative_to_sun ecliplon_relative_to_opposition ecliplat angvel_ecliplon angvel_ecliplat MJD object obscode\n";
  objname=oldname=obscode=oldcode=band="";
  oldmjd = arc = timespan = GCR = 0.0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the mpc file
    getline(instream1,lnfromfile);
    objname=obscode=band="";
    MJD = 0.0;
    RA = Dec = -LARGERR2;
    if(verbose>0) cout << "Read line " << lnfromfile << "\n";
    status = mpc80_parseline(lnfromfile, objname, &MJD, &RA, &Dec, &mag, band, obscode);
    if(verbose>0) cout << "status, objname, MJD, RA, Dec, mag, band, obscode: " << status << " " << objname << " " << MJD << " " << RA << " " << Dec << " " << mag << " " << band << " " << obscode << "\n";
    if(status==0 && objname.size()>0 && obscode.size()>0 && MJD>0.0 && RA>=0.0 && RA<360.0 && Dec>=-90.0 && Dec<90.0) {
      // This line seems to have been read successfully
      // Convert RA and Dec to ecliptic coordinates.
      poleswitch02(RA, Dec, NEPRA, NEPDEC, 90.0l, ecliplon, ecliplat);
      if(verbose>0) cout << "RA, Dec, ecliplon,ecliplat = " << RA << " " << Dec << " " << ecliplon << " " << ecliplat << "\n";
      if(verbose>0) cout << "objname, oldname, obscode, oldcode: " << objname << " " << oldname << " " << obscode << " " << oldcode << " " << MJD << " " << oldmjd << "\n";
      if(objname==oldname && obscode==oldcode && fabs(MJD-oldmjd)<NIGHTSTEP) {
	// This line appears to continue a previously started tracklet.
	// Add the new point to the tracklet vectors.
	MJDvec.push_back(MJD);
	RAvec.push_back(ecliplon);
	Decvec.push_back(ecliplat);
      } else {
	// We appear to have transitioned to a new tracklet.
	if(verbose>0) cout << "Finishing out a tracklet with " << MJDvec.size() << " points\n";
	if(long(MJDvec.size())>=minpts) {
	// Analyze the newly completed tracklet.
	  arc = timespan = GCR = 0.0;
	  status = greatcircfit(MJDvec, RAvec, Decvec, poleRA, poleDec, angvel, pa, crosstrack, alongtrack);
	  if(status==0) {
	    // Great Circle fit succeeded. Calculate angular arc, timespan, and GCR.
	    timespan = MJDvec[MJDvec.size()-1] - MJDvec[0];
	    arc = angvel*timespan*3600.0;
	    GCR = sqrt(crosstrack*crosstrack + alongtrack*alongtrack);
	    if(verbose>0) cout << "angvel, pa, timespan, arc, GCR = " << angvel << " " << pa << " " << timespan << " " << arc << " " << GCR << "\n";
	    if(timespan>=mintime && arc>=minarc && GCR<=maxGCR) {
	      // This tracklet is good. Calculate the mean ecliptic coordinates
	      lonmed = lonmean = lonrms = latmean = MJDmean = 0.0;
	      anglevec_meanrms(RAvec, 360.0, &lonmed, &lonmean, &lonrms);
	      latmean = dmean01(Decvec);
	      // Calculate the mean MJD and the solar ecliptic coordinates at that time.
	      MJDmean = dmean01(MJDvec);
	      sunradec_approx01(MJDmean, sunra, sundec);
	      poleswitch02(sunra, sundec, NEPRA, NEPDEC, 90.0l, sunlon, sunlat);
	      // Calculate the object's distance from the sun in ecliptic longitude
	      ecliplon_relative_to_sun = lonmean - sunlon;
	      while(ecliplon_relative_to_sun <= -180.0) ecliplon_relative_to_sun += 360.0;
	      while(ecliplon_relative_to_sun > 180.0) ecliplon_relative_to_sun -= 360.0;
	      // Calculate the object's distance from opposition in ecliptic longitude
	      ecliplon_relative_to_opposition = lonmean - sunlon-180.0;
	      while(ecliplon_relative_to_opposition <= -180.0) ecliplon_relative_to_opposition += 360.0;
	      while(ecliplon_relative_to_opposition > 180.0) ecliplon_relative_to_opposition -= 360.0;

	      if(verbose>0) cout << "MJDmean, sunra, sundec, sunlon, sunlat, sunelong: " << MJDmean << " " << sunra << " " << sundec << " " << sunlon << " " << sunlat << " " << ecliplon_relative_to_sun << "\n";
	      // Write sun-relative coordinates and ecliptic angular velocities to output file
	      outstream1 << fixed << setprecision(6) << ecliplon_relative_to_sun << " " << ecliplon_relative_to_opposition << " " << latmean << " " << angvel*sin(pa/DEGPRAD) << " " << angvel*cos(pa/DEGPRAD) << " " << MJDmean << " " << oldname << " " << oldcode << "\n";
	    }
	  }
	}
	// Set up to load the new tracklet.
	MJDvec={};
	RAvec={};
	Decvec={};
	oldname = objname;
	oldcode = obscode;
	oldmjd = MJD;
	// Load the current point as the first point of the new tracklet.
	MJDvec.push_back(MJD);
	RAvec.push_back(ecliplon);
	Decvec.push_back(ecliplat);
      }
    }
  }
  instream1.close();
  outstream1.close();
  return(0);
}
