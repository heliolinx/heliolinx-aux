// February 20, 2025: find_detections01a.cpp
// Given an input detection file in hldet format, and an input
// ephemeris file in MJD RA Dec format, find detections corresponding
// to ephemeris points to within input tolerances for time and on-sky distance.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: find_detections01a -obs file of observations to be matched -eph ephemeris file -timetol MJD_matching_tolerance(days) -skytol sky_matching_radius(deg) -outfile output file\n";
}
  
int main(int argc, char *argv[])
{
  long i,detnum,detct,ephnum,ephct;
  string obsfile,ephfile,outfile,stest,lnfromfile;
  double timetol=0.0001;
  double skytol=2.0;
  double MJD,RA,Dec,dist;
  MJD = RA = Dec = dist = 0.0;
  vector <hldet> detvec = {};
  vector <hldet> outdetvec = {};
  vector <hldet> outdettemp = {};
  vector <det_bsc> ephvec = {};
  det_bsc onedet = det_bsc(0,0,0);
  ifstream instream1;
  ofstream outstream1;
  int badread=0;
  int reachedeof=0;
  int startpoint=0;
  int endpoint=0;
  int status=0;
  int verbose=0;

  if(argc<9) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-det" || string(argv[i]) == "-d" || string(argv[i]) == "-o" || string(argv[i]) == "-obs" || string(argv[i]) == "-obsfile" || string(argv[i]) == "-observations" || string(argv[i]) == "--obsfile" || string(argv[i]) == "--observation" || string(argv[i]) == "--observationfile" || string(argv[i]) == "--observations") {
      if(i+1 < argc) {
	// There is still something to read;
	obsfile=argv[++i];
	i++;
      } else {
	cerr << "Observation file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-eph" || string(argv[i]) == "-e" || string(argv[i]) == "-eff" || string(argv[i]) == "-ephem" || string(argv[i]) == "-ephemeris" || string(argv[i]) == "-ephemerisfile" || string(argv[i]) == "--eph" || string(argv[i]) == "--ephemeris_file" || string(argv[i]) == "--ephemeris_file" || string(argv[i]) == "--ephemeris") {
      if(i+1 < argc) {
	// There is still something to read;
	ephfile=argv[++i];
	i++;
      } else {
	cerr << "Observation file keyword supplied with no corresponding argument\n";
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
    }
  }
  
  // Read input observation file.
  detvec={};
  status=read_hldet_file(obsfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << obsfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  detnum = detvec.size();
  outdetvec = detvec;
  cout << "Read " << detnum << " data lines from observation file " << obsfile << "\n";

  instream1.open(ephfile);
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Read on.
      // Read the MJD
      startpoint=0;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { MJD = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read MJD string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      else badread=1;
      // Read the RA
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { RA = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read RA string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the Dec
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { Dec = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read Dec string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      onedet = det_bsc(MJD,RA,Dec);
      ephvec.push_back(onedet);
    } else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(badread!=0) {
      cerr << "ERROR reading observation file " << obsfile << "\n";
      cerr << "Last point was " << detvec.size() << "; last file line was " << lnfromfile << "\n";
    }
  }
  instream1.close();
  ephnum = ephvec.size();
  cout << "Read " << ephnum << " lines from ephemeris file " << ephfile << "\n";
  
  // Sort both input vectors
  sort(detvec.begin(),detvec.end(),early_hldet());
  sort(ephvec.begin(),ephvec.end(),early_det());

  // Perform search, sorting matches in outdetvec
  outdetvec={};
  ephct=detct=0;
  while(ephct<ephnum && detct<detnum) {
    while(detct<detnum && detvec[detct].MJD < ephvec[ephct].MJD-timetol) {
      // Skip forward to find detections that are in the right time range
      detct++;
    }
    // Search through the detection vector for possible matches
    i=detct;
    while(i<detnum && detvec[i].MJD <= ephvec[ephct].MJD+timetol) {
      dist = 3600.0*distradec01(ephvec[ephct].RA, ephvec[ephct].Dec, detvec[i].RA, detvec[i].Dec);
      //cout << fixed << setprecision(7) << "skytol = " << skytol << ", dist = " << dist << " for comp " << ephvec[ephct].MJD << " " << ephvec[ephct].RA << " " << ephvec[ephct].Dec << " " << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << "\n";
      if(dist<=skytol) {
	// Load this detection to the output list of matches.
	outdetvec.push_back(detvec[i]);
      }
      i++;
    }
    ephct++;
  }
  outstream1.close();
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
      cout << fixed << setprecision(7) << outdetvec[i].MJD << " " << outdetvec[i].RA << " " << outdetvec[i].Dec << " ";
      cout << fixed << setprecision(4) << outdetvec[i].mag << " ";
      cout << outdetvec[i].image << " " << outdetvec[i].idstring << " " << outdetvec[i].index << "\n";
    }
  }
  outstream1.close();
  cout << outdetvec.size() << " matching detections were found\n";
  return(0);
}
