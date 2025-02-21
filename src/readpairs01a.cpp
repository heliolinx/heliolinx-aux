// February 09, 2023: readpairs01a.cpp:
// Re-written from original, November 09, 2021, version.
// Read pair files output by make_tracklets, and calculate
// the angular velocities, etc.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define TYPECODE_PAIR "P"
#define TYPECODE_TRACKLET "T"

static void show_usage()
{
  cerr << "Usage: readpairs01a -dets detfile -pairs pairfile -out outfile\n";
}
    
int main(int argc, char *argv[])
{
  string indetfile;
  string inpairfile;
  string outfile;
  string lnfromfile;
  
  double MJD,RA, Dec;
  int reachedeof=0;
  det_obsmag_indvec o1 = det_obsmag_indvec(0L,0l,0l,0L,0L,0L,"null",0l,"V","I11",0,{});
  vector <det_obsmag_indvec> detvec = {};
  vector <long> pindexvec={};
  vector <vector <long>> pairvec={};
  double mag=0l;
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  string typecode;
  long double X=0.0L;
  long double Y=0.0L;
  long double Z=0.0L;
  char detid[SHORTSTRINGLEN];
  long origind=0;
  long i1=0;
  long i2=0;
  ifstream instream1;
  long detfilelinect;
  int badread,startpoint,endpoint;
  int ipt,trackpointct,trackpointnum;
  string stest;
  double dist,pa,timediff;
  
  if(argc<3)
    {
      show_usage();
      return(1);
    }
  
  int i=1;
  int j=0;
  int c='0';
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" || string(argv[i]) == "--detections") {
      if(i+1 < argc) {
	//There is still something to read;
	indetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-p" || string(argv[i]) == "-pair" || string(argv[i]) == "-pairs" || string(argv[i]) == "--pairs" || string(argv[i]) == "--pair" || string(argv[i]) == "--pairfile" || string(argv[i]) == "--pairsfile") {
      if(i+1 < argc) {
	//There is still something to read;
	inpairfile=argv[++i];
	i++;
      }
      else {
	cerr << "Pair file keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outpair" || string(argv[i]) == "--outpairs") {
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
    }
  }
    
  cout << "input detection file " << indetfile << "\n";
  cout << "input pair file " << inpairfile << "\n";
  cout << "output file " << outfile << "\n";
  
  // Read input detection file.
  instream1.open(indetfile,ios_base::in);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << indetfile << "\n";
    return(1);
  }
  detvec={};
  // Skip header line
  getline(instream1,lnfromfile);
  cout << "Header line from input paired detection file " << indetfile << ":\n";
  cout << lnfromfile << "\n";
  // Read body of the file
  detfilelinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the paired detections file, and load an object of class det_obsmag_indvec
    getline(instream1,lnfromfile);
    detfilelinect++;
    badread=0;
    if(lnfromfile.size()>60) {
      // Read MJD, RA, Dec, observer x, y, z
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) MJD = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) RA = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Dec = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) X = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Y = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Z = stold(stest);
      else badread=1;
      // Read the IDstring
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(detid,stest,SHORTSTRINGLEN);
      else badread=1;
      // Read the magnitude
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) mag = stod(stest);
      else badread=1;
      // Read the band and observatory code
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(band,stest,MINSTRINGLEN);
      else badread=1;
       startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(obscode,stest,MINSTRINGLEN);
      else badread=1;
      // Read the original detection index
       startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) origind = stol(stest);
      else badread=1;

      // If there was a file read error, abort.
      if(badread==1) {
	cerr << "ERROR reading line " << detfilelinect << " of paired detection file" << indetfile << "\n";
	return(1);
      }
      // If we reach this point, the line was read OK. Write it to detvec.
      o1=det_obsmag_indvec(MJD,RA,Dec,X,Y,Z,detid,mag,band,obscode,origind,{});
      detvec.push_back(o1);
    } else if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      cerr << "WARNING: line " << detfilelinect << " of paired detection file " << indetfile << " was too short\n";
    }
  }
  instream1.close();
  cout << detvec.size() << " detection records read from " << indetfile << ".\n";
  
  // Read input image pair file
  pairvec={};
  instream1.open(inpairfile);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << inpairfile << "\n";
    return(1);
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    typecode = ""; // Wipe previously read typecode
    // Read the current type code
    instream1 >> typecode;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad() && typecode == TYPECODE_PAIR) {
      // Read a single, isolated pair
      pindexvec={};
      instream1 >> i1 >> i2;
      pindexvec.push_back(i1);
      pindexvec.push_back(i2);
      pairvec.push_back(pindexvec);
    } else if(!instream1.eof() && !instream1.fail() && !instream1.bad() && typecode == TYPECODE_TRACKLET) {
      // Read the representative pair
      pindexvec={};
      instream1 >> i1 >> i2;
      pindexvec.push_back(i1);
      pindexvec.push_back(i2);
      // Read revised RA, Dec for representative detection 1.
      instream1 >> RA >> Dec;
      // Re-assign paired detection RA, Dec
      detvec[i1].RA = RA;
      detvec[i1].Dec = Dec;
      // Read revised RA, Dec for representative detection 2.
      instream1 >> RA >> Dec;
      // Re-assign paired detection RA, Dec
      detvec[i2].RA = RA;
      detvec[i2].Dec = Dec;
      // Read number of points in the tracklet
      instream1 >> trackpointnum;
      for(trackpointct=0; trackpointct<trackpointnum; trackpointct++) {
	instream1 >> ipt;
	if(ipt!=i1 && ipt!=i2) pindexvec.push_back(ipt);
      }
      pairvec.push_back(pindexvec);
    } else if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      cerr << "ERROR: unrecognized pair type code " << typecode << "\n";
      return(1);
    }
  }
  instream1.close();
  cout << "Read " << detvec.size() << " detections and " << pairvec.size() << "pairs.\n";

  ofstream outstream1 {outfile};
  outstream1 << "MJD,RA,Dec,timediff,arc,angvelE,angvelN\n";
  for(i=0;i<pairvec.size();i++) {
    i1=pairvec[i][0];
    i2=pairvec[i][1];
    distradec02(detvec[i1].RA,detvec[i1].Dec,detvec[i2].RA,detvec[i2].Dec,&dist,&pa);
    timediff = detvec[i2].MJD - detvec[i1].MJD;
    outstream1 << fixed << setprecision(6) << detvec[i1].MJD << "," << detvec[i1].RA << "," << detvec[i1].Dec << "," << timediff << "," << dist << "," << dist/timediff*sin(pa/DEGPRAD) << "," << dist/timediff*cos(pa/DEGPRAD) << "\n";
    // cout << fixed << setprecision(6) << detvec[i1].MJD << "," << detvec[i1].RA << "," << detvec[i1].Dec << "," << timediff << "," << dist << "," << dist/timediff*sin(pa/DEGPRAD) << "," << dist/timediff*cos(pa/DEGPRAD) << "\n";
  }
  outstream1.close();
  
  return(0);
}
