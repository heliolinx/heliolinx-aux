// April 16, 2024: match_itf01.cpp.
// Read a tracklet catalog with columns
// idstring MJD RA Dec angvel PA RAvel Decvel avgmag pnum GCR alongtrack crosstrack
// and an emphemeris file with columns
// MJD RA Dec RAvel Decvel mag,
// and find potential matches in the tracklet catalog.

#include "solarsyst_dyn_geo01.h"
#include "cmath"


static void show_usage()
{
  cerr << "Usage: match_itf01 -trk tracklet_file -eph ephemeris_file -outfile outfile -timetol MJD_matching_tolerance(days) -skytol sky_matching_radius(deg) -veltol velocity_matching_radius(deg/day)\n";
}
int main(int argc, char *argv[])
{
  string incat,ineph,outfile,lnfromfile,sdat;
  vector <string> trkline;
  
  int reachedeof=0;
  ifstream instream1;
  ofstream outstream1;
  double MJD,RA,Dec,RAvel,Decvel,mag;
  MJD = RA = Dec = 0.0l;
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
  vector <double> ephmag;
  char c;
  double timetol=1.0;
  double skytol=1.0;
  double veltol=0.1;
  double dist=0.0;
  double vdist=0.0;
  int goodline=1;
  long i=0;
  long j=0;

  if(argc<7) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-trk" || string(argv[i]) == "-incat" ) {
      if(i+1 < argc) {
	//There is still something to read;
        incat=argv[++i];
	i++;
      }
      else {
	cerr << "Input tracklet catalog keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-eph" || string(argv[i]) == "-ineph" ) {
      if(i+1 < argc) {
	//There is still something to read;
	ineph=argv[++i];
	i++;
      }
      else {
	cerr << "Input ephmeris file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" ) {
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
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  instream1.open(incat);
  getline(instream1,lnfromfile); // Skip header line
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    goodline=1;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    // idstring
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      i++;
    }
    // MJD
    sdat = {};
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      if(c!=' ' && c!=',') sdat.push_back(c);
      i++;
    }
    try {
      MJD = stod(sdat);
    } catch(std::invalid_argument const& ex) {
      cerr << "WARNING: unable to read MJD from line " << lnfromfile << "\n";
      goodline=0;
    }
    // RA
    sdat = {};
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      if(c!=' ' && c!=',') sdat.push_back(c);
      i++;
    }
    try {
      RA = stod(sdat);
    } catch(std::invalid_argument const& ex) {
      cerr << "WARNING: unable to read RA from line " << lnfromfile << "\n";
      goodline=0;
    }
    // Dec
    sdat = {};
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      if(c!=' ' && c!=',') sdat.push_back(c);
      i++;
    }
    try {
      Dec = stod(sdat);
    } catch(std::invalid_argument const& ex) {
      cerr << "WARNING: unable to read Dec from line " << lnfromfile << "\n";
      goodline=0;
    }
    // angvel
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      i++;
    }
    // PA
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      i++;
    }
    // RAvel
    sdat = {};
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      if(c!=' ' && c!=',') sdat.push_back(c);
      i++;
    }
    try {
      RAvel = stod(sdat);
    } catch(std::invalid_argument const& ex) {
      cerr << "WARNING: unable to read RAvel from line " << lnfromfile << "\n";
      goodline=0;
    }
    // Decvel
    sdat = {};
    c='a';
    while(i<long(lnfromfile.size()) && c!=' ' && c!=',') {
      c = lnfromfile[i];
      if(c!=' ' && c!=',') sdat.push_back(c);
      i++;
    }
    try {
      Decvel = stod(sdat);
    } catch(std::invalid_argument const& ex) {
      cerr << "WARNING: unable to read Decvel from line " << lnfromfile << "\n";
      goodline=0;
    }
    if(goodline==1) {
      trkline.push_back(lnfromfile);
      trkMJD.push_back(MJD);
      trkRA.push_back(RA);
      trkDec.push_back(Dec);
      trkRAvel.push_back(RAvel);
      trkDecvel.push_back(Decvel);
    }
  }
  cout << "Vector sizes: " << trkline.size() << " "  << trkMJD.size() << " "   << trkRA.size() << " "   << trkDec.size() << " "   << trkRAvel.size() << " "   << trkDecvel.size() << "\n";

  instream1.close();
  getline(instream1,lnfromfile); // Skip header line
  instream1.open(ineph);
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    MJD = 0.0l;
    mag = 0.0l;
    instream1 >> MJD >> RA >> Dec >> RAvel >> Decvel >> mag;
    if(MJD>0.0l) {
      ephMJD.push_back(MJD);
      ephRA.push_back(RA);
      ephDec.push_back(Dec);
      ephRAvel.push_back(RAvel);
      ephDecvel.push_back(Decvel);
      ephmag.push_back(mag);
    }
  }
  instream1.close();
  cout << "Vector sizes: " << ephMJD.size() << " "   << ephRA.size() << " "   << ephDec.size() << " "   << ephRAvel.size() << " " << ephDecvel.size() << " " << ephmag.size() << "\n";

  outstream1.open(outfile);
  // Loop on the ephemeris vectors
  j=0;
  for(i=0;i<long(ephMJD.size());i++) {
    while(j<long(trkMJD.size()) && trkMJD[j] < ephMJD[i]-timetol) j++;
    if(j<long(trkMJD.size())) {
      // We are set to explore a region of the tracklet vectors
      // that is appropriate for finding matches to ephemeris point i.
      long k=j;
      while(k<long(trkMJD.size()) && trkMJD[k] < ephMJD[i]+timetol) {
	// Tracklet catalog entry k is in the time range for matching
	// to ephemeris point i. Check its match in other respects.
	dist = distradec01(ephRA[i],ephDec[i],trkRA[k],trkDec[k]);
	if(dist<=skytol) {
	  // Tracklet is in the right part of the sky. Test its velocity.
	  vdist = sqrt(DSQUARE(ephRAvel[i]-trkRAvel[k])+DSQUARE(ephDecvel[i]-trkDecvel[k]));
	  if(vdist<=veltol) {
	    // Tracklet is a perfect match. Write it to the output file
	    cout << "Match found on tracklet line " << k << ":\n";
	    cout << trkline[k] << " :\n" << ephMJD[i] << " " << ephRA[i] << " " << ephDec[i] << " " << ephRAvel[i] << " " << ephDecvel[i] << " " << ephmag[i] << "\n";
	    outstream1 << trkline[k] << " : " << ephMJD[i] << " " << ephRA[i] << " " << ephDec[i] << " " << ephRAvel[i] << " " << ephDecvel[i] << " " << ephmag[i] << "\n";
	  }
	}
	k++;
      }
    }
  }
  outstream1.close();
  return(0);
}
