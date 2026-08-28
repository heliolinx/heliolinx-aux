// July 16, 2026: Convert an ADES psv file into hldet format

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define MAX_READ_ATTEMPTS 100

static void show_usage()
{
  cerr << "Usage: ADESpsv2hldet -infile ADESfile -timeshift timeshift -out outfile\n";
}

int main(int argc, char *argv[])
{
  hldet onedet;
  vector <hldet> detvec = {};
  string ADESfile,outfile;
  string lnfromfile,stest;
  long i=0;
  long read_attempts=0;
  int verbose=0;
  int badread=0;
  int reachedeof=0;
  int startpoint=0;
  int endpoint=0;
  int colct=0;
  int year, month, day, hour, minute;
  year = month = day = hour = minute = 0;
  double second, MJD;
  second = MJD = 0.0;
  double RA, Dec, mag;
  RA = Dec = mag = 0.0;
  string signstring;
  long double timeshift = 0.0l;
  ofstream outstream1;
  ifstream instream1;
  string permID,provID,artSat;
  string trkSub="a000000";
  int trkSub_col,obsSubID_col,mode_col,stn_col,obsTime_col,ra_col,dec_col,rmsRA_col,rmsDec_col,astCat_col,mag_col,rmsMag_col,band_col;
  trkSub_col = obsSubID_col = mode_col = stn_col = obsTime_col = ra_col = dec_col = rmsRA_col = rmsDec_col = astCat_col = mag_col = rmsMag_col = band_col = 0;
  string trkSub_string, time_string, ra_string, dec_string, mag_string, band_string, ObsCode_string, rmsRA_string, rmsDec_string, rmsMag_string;
  double trail_len, trail_PA, sigmag, sig_across, sig_along;
  trail_len = trail_PA = sigmag = sig_across = sig_along = 0.0;
  int image = 0; 
  long known_obj,det_qual,index;
  known_obj = det_qual = index = 0;
  int trkSub_found,ra_found,dec_found,obsTime_found;
  trkSub_found = ra_found = dec_found = obsTime_found = 0;
  
  if(argc<7) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-hlfile" || string(argv[i]) == "-hldet" || string(argv[i]) == "-ADESfile") {
      if(i+1 < argc) {
	//There is still something to read;
	ADESfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input hldet file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timeshift" || string(argv[i]) == "-time" || string(argv[i]) == "--timeshift" || string(argv[i]) == "-timediff") {
      if(i+1 < argc) {
	//There is still something to read;
	timeshift = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--output_ADES_file") {
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

  // Catch required parameters if missing
  if(ADESfile.size()<=0) {
    cerr << "\nERROR: input ADES psv file is required\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cerr << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  cout << "input ADES psv file " << ADESfile << "\n";
  cout << "time shift " << timeshift << " seconds\n";
  cout << "output file " << outfile << "\n";

  // Read input ADS psv file.
  instream1.open(ADESfile);
  if(!instream1) {
    cerr << "can't open input file " << ADESfile << "\n";
    return(1);
  }
  // Find the data header line that includes the strings trkSub, ra, dec, and obsTime;
  reachedeof = badread = 0;
  trkSub_found = ra_found = dec_found = obsTime_found = 0;
  stest="";
  while(reachedeof==0 && (trkSub_found==0 || ra_found == 0 || dec_found == 0 || obsTime_found==0)) {
    getline(instream1,lnfromfile);  
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Look for the string "trkSub" anywhere in this line
      stest="";
      i=0;
      while(stest!="trkSub" && i<long(lnfromfile.size()-6)) {
	stest = lnfromfile.substr(i,6);
	if(stest=="trkSub") trkSub_found=1;
	i++;
      }
      // Look for the string "ra" anywhere in this line
      stest="";
      i=0;
      while(stest!="ra" && i<long(lnfromfile.size()-2)) {
	stest = lnfromfile.substr(i,2);
	if(stest=="ra") ra_found=1;
	i++;
      }
      // Look for the string "dec" anywhere in this line
      stest="";
      i=0;
      while(stest!="dec" && i<long(lnfromfile.size()-3)) {
	stest = lnfromfile.substr(i,3);
	if(stest=="dec") dec_found=1;
	i++;
      }
      // Look for the string "obsTime" anywhere in this line
      stest="";
      i=0;
      while(stest!="obsTime" && i<long(lnfromfile.size()-7)) {
	stest = lnfromfile.substr(i,7);
	if(stest=="obsTime") obsTime_found=1;
	i++;
      }
    } else badread=1;
  }
  if(badread==1 || trkSub_found==0 || ra_found == 0 || dec_found == 0 || obsTime_found==0) {
    cerr << "ERROR: could not find valid column-header line\n";
    return(1);
  } else if(trkSub_found==1 && ra_found == 1 && dec_found == 1 && obsTime_found==1) {
    cout << "Found a valid header line:\n";
    cout << lnfromfile << "\n";
  }
  // We assume we have found a line that includes the trkSub field.
  colct = endpoint = startpoint = 0;
  while(startpoint<int(lnfromfile.size()) && endpoint<=int(lnfromfile.size())) {
    endpoint = get_psv_string01(lnfromfile,stest,startpoint);
    startpoint = endpoint+1;
    colct++;
    if(stest=="trkSub") trkSub_col = colct;
    else if(stest=="obsSubID") obsSubID_col = colct;
    else if(stest=="mode") mode_col = colct;
    else if(stest=="stn") stn_col = colct;
    else if(stest=="obsTime") obsTime_col = colct;
    else if(stest=="ra") ra_col = colct;
    else if(stest=="dec") dec_col = colct;
    else if(stest=="rmsRA") rmsRA_col = colct;
    else if(stest=="rmsDec") rmsDec_col = colct;
    else if(stest=="astCat") astCat_col = colct;
    else if(stest=="mag") mag_col = colct;
    else if(stest=="rmsMag") rmsMag_col = colct;
    else if(stest=="band") band_col = colct;
  }
  // Check that required columns were found
  if(stn_col<=0 || obsTime_col<=0 || ra_col<=0 || dec_col<=0 || stn_col<=0) {
    cerr << "One of the required columns is missing from the input file column header\n";
    cerr << "stn_col = " << stn_col << " obsTime_col = " << obsTime_col << " ra_col = " << ra_col << " dec_col = " << dec_col  << " stn_col = " << stn_col << "\n";
    return(2);
  }
  if(verbose>0) cout << "stn_col = " << stn_col << " obsTime_col = " << obsTime_col << " ra_col = " << ra_col << " dec_col = " << dec_col  << " stn_col = " << stn_col << "\n";
  // If we got here, we assume we have successfully read the column header line.
  // Read and convert all of the other lines.
  reachedeof = badread = 0;
  while(reachedeof==0) {
    badread=0;
    getline(instream1,lnfromfile);
    year = month = day = hour = minute = 0;
    second = MJD = mag = 0.0;
    RA = -999.999;
    Dec = -99.999;
    sigmag = -9.99;
    sig_across = -9.99;
    sig_along = -9.99;
    if(verbose>0) cout << "Line read: " << lnfromfile << "\n";
    colct=startpoint=endpoint=0;
    read_attempts=0;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      while(startpoint<int(lnfromfile.size()) && endpoint<=int(lnfromfile.size()) && read_attempts<MAX_READ_ATTEMPTS) {
	endpoint = get_psv_string01(lnfromfile,stest,startpoint);
	startpoint = endpoint+1;
	colct++;
	read_attempts++;
	if(colct==trkSub_col) trkSub_string = stest;
	else if(colct==stn_col) ObsCode_string = stest;
	else if(colct==obsTime_col) time_string = stest;
	else if(colct==ra_col) ra_string = stest;
	else if(colct==dec_col) dec_string = stest;
	else if(colct==mag_col) mag_string = stest;
	else if(colct==band_col) band_string = stest;
	else if(colct==rmsRA_col) rmsRA_string = stest;
	else if(colct==rmsDec_col) rmsDec_string = stest;
	else if(colct==rmsMag_col) rmsMag_string = stest;
      }
      if(verbose>0) cout << "String read: " << trkSub_string << " " << ObsCode_string  << " " << time_string  << " " << ra_string  << " " <<  dec_string  << " " <<  mag_string << " " <<  band_string << "\n";
      if(time_string.size()<24 || ra_string.size()<5 || dec_string.size()<4 || ObsCode_string.size()<=0) {
	cerr << "Warning: read an invalid line\n";
	continue;
      }
      if(verbose>0) cout << "Timestrings: " << time_string.substr(0,4) << " " << time_string.substr(5,2)  << " " << time_string.substr(8,2)  << " " << time_string.substr(11,2) << " " << time_string.substr(14,2) << " " << time_string.substr(17,6) << "\n";
      // 2025-06-05T07:10:40.802Z
      try { year = stoi(time_string.substr(0,4)); }
	catch(...) { cerr << "ERROR: cannot read obsTime string " << time_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      try { month = stoi(time_string.substr(5,2)); }
	catch(...) { cerr << "ERROR: cannot read obsTime string " << time_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      try { day = stoi(time_string.substr(8,2)); }
	catch(...) { cerr << "ERROR: cannot read obsTime string " << time_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      try { hour = stoi(time_string.substr(11,2)); }
	catch(...) { cerr << "ERROR: cannot read obsTime string " << time_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      try { minute = stoi(time_string.substr(14,2)); }
	catch(...) { cerr << "ERROR: cannot read obsTime string " << time_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      try { second = stod(time_string.substr(17,6)); }
	catch(...) { cerr << "ERROR: cannot read obsTime string " << time_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      // Read the RA
      try { RA = stod(ra_string); }
	catch(...) { cerr << "ERROR: cannot read RA string " << ra_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      // Read the Dec
      try { Dec = stod(dec_string); }
	catch(...) { cerr << "ERROR: cannot read Dec string " << dec_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      if(mag_col>0) {
	// Read the mag
	try { mag = stod(mag_string); }
	catch(...) { cerr << "ERROR: cannot read mag string " << mag_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      if(rmsRA_col>0 && rmsRA_string.size()>0) {
	// Read the RA rms, load into hldet sig_along
	try { sig_along = stod(rmsRA_string); }
	catch(...) { cerr << "ERROR: cannot read rmsRA string " << rmsRA_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      if(rmsDec_col>0 && rmsDec_string.size()>0) {
	// Read the Dec RMS, load into hldet sig_across
	try { sig_across = stod(rmsDec_string); }
	catch(...) { cerr << "ERROR: cannot read rmsDec string " << rmsDec_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      if(rmsMag_col>0 && rmsMag_string.size()>0) {
	// Read the mag RMS, load into hldet sigmag
	try { sigmag = stod(rmsMag_string); }
	catch(...) { cerr << "ERROR: cannot read rmsMag string " << rmsMag_string << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      if(badread==0) {
	// Convert calendar date to MJD
	MJD = cal2MJD(year, month, day, hour, minute, second);
	MJD += timeshift/SOLARDAY;
	// Change two-character LSST bands to one-character, if needed.
	if(band_string=="Lu") band_string="u";
	else if(band_string=="Lg") band_string="g";
	else if(band_string=="Lr") band_string="r";
	else if(band_string=="Li") band_string="i";
	else if(band_string=="Lz") band_string="z";
	else if(band_string=="Ly") band_string="y";
	// Convert available quantities to hldet format
	onedet = hldet(MJD, RA, Dec, mag, trail_len, trail_PA, sigmag, sig_across, sig_along, image, trkSub_string, band_string, ObsCode_string, known_obj, det_qual, index);
	detvec.push_back(onedet);
      } else {
	cerr << "ERROR: bad read on line " << lnfromfile << "\n";
	return(3);
      }
    } else if(instream1.eof() || instream1.fail() || instream1.bad()) reachedeof=1;
  }
  instream1.close();
  cout << "Read " << detvec.size() << " lines from input file " << ADESfile << "\n";
  
  outstream1.open(outfile);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(i=0;i<long(detvec.size());i++) {
    outstream1 << fixed << setprecision(9) << detvec[i].MJD << "," << detvec[i].RA << "," << detvec[i].Dec << ",";
    outstream1 << fixed << setprecision(6) << detvec[i].mag << ",";
    outstream1 << fixed << setprecision(2) << detvec[i].trail_len << "," << detvec[i].trail_PA << ",";
    outstream1 << fixed << setprecision(6) << detvec[i].sigmag << ",";
    outstream1 << fixed << setprecision(6) << detvec[i].sig_across << "," << detvec[i].sig_along << ",";
    outstream1 << detvec[i].image << "," << detvec[i].idstring << "," << detvec[i].band << ",";
    outstream1 << detvec[i].obscode << "," << detvec[i].known_obj << ","; 
    outstream1 << detvec[i].det_qual << "," << detvec[i].index << "\n"; 
  }
  outstream1.close();
  return(0);
}
