// January 25, 2025: MJDcode2statevec.cpp:
// 
// 
// 
// 

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: MJDcode2statevec -infile infile -earth earthfile -obscode obscodefile -outfile outfile\n";
}

int main(int argc, char *argv[])
{
  long i=0;
  int verbose=0;
  int reachedeof=0;
  char c = '0';
  string infile,outfile;
  string stest;
  ifstream instream1;
  ofstream outstream1;
  string lnfromfile;
  double MJD, RA, Dec;
  int status=0;
  string teststring;
  string earthfile;
  string obscodefile;
  vector <point3d> Earthpos;
  vector <point3d> Earthvel;
  vector <double> EarthMJD;
  point3d obspos = point3d(0.0,0.0,0.0);
  point3d obsvel = point3d(0.0,0.0,0.0);
  vector <observatory> observatory_list = {};
  char obscode[MINSTRINGLEN];
  double obslon,plxcos,plxsin;

  if(argc<0 && argc<1) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-in" || string(argv[i]) == "-filein" || string(argv[i]) == "-input" || string(argv[i]) == "--infile" || string(argv[i]) == "-MJDinfile" || string(argv[i]) == "--mjdfile" || string(argv[i]) == "-mjdin" || string(argv[i]) == "--MJD") {
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
    } else if(string(argv[i]) == "-earth" || string(argv[i]) == "-e" || string(argv[i]) == "-Earth" || string(argv[i]) == "--earthfile" || string(argv[i]) == "--Earthfile" || string(argv[i]) == "--earth" || string(argv[i]) == "--Earth") {
      if(i+1 < argc) {
	//There is still something to read;
	earthfile=argv[++i];
	i++;
      }
      else {
	cerr << "Earth file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--out" || string(argv[i]) == "-o") {
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

  // Check for required arguments
  if(infile.size()<=0) {
    cerr << "ERROR: you have not specified an input file\n";
    show_usage();
    return(1);
  }
  if(outfile.size()<=0) {
    cerr << "ERROR: you have not specified an output file\n";
    show_usage();
    return(1);
  }

  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, verbose);
  if(status!=0) {
    cerr << "ERROR reading observatory code file " << obscodefile << "\n";
    return(1);
  }
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";
 
  EarthMJD={};
  Earthpos={};
  Earthvel={};
  read_horizons_csv(earthfile,EarthMJD,Earthpos,Earthvel);
  cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";
  
  // Read input file of observations
  instream1.open(infile,ios_base::in);
  cout << "Reading file " << infile << "\n";
  outstream1.open(outfile);
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the input file
    getline(instream1,lnfromfile);
    if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      teststring = "";
      c='0';
      i=0;
      while(i<long(lnfromfile.size()) && reachedeof == 0 && c!='=' && c!=' ' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=' ' && c!='=' && c!='\n' && c!=EOF) teststring.push_back(c);
	if(c==EOF) reachedeof=1;
	i++;
      }
      MJD=stod(teststring);
      teststring = "";
      c='0';
      while(i<long(lnfromfile.size()) && reachedeof == 0 && c!='=' && c!=' ' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=' ' && c!='=' && c!='\n' && c!=EOF) teststring.push_back(c);
	if(c==EOF) reachedeof=1;
	i++;
      }
      RA=stod(teststring);
      teststring = "";
      c='0';
      while(i<long(lnfromfile.size()) && reachedeof == 0 && c!='=' && c!=' ' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=' ' && c!='=' && c!='\n' && c!=EOF) teststring.push_back(c);
	if(c==EOF) reachedeof=1;
	i++;
      }
      Dec=stod(teststring);
      teststring = "";
      c='0';
      while(i<long(lnfromfile.size()) && reachedeof == 0 && c!='=' && c!=' ' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=' ' && c!='=' && c!='\n' && c!=EOF) teststring.push_back(c);
	if(c==EOF) reachedeof=1;
	i++;
      }
      stringncopy01(obscode,teststring,MINSTRINGLEN);
      status = obscode_lookup(observatory_list,obscode,obslon,plxcos,plxsin);
      if(status>0) {
	cerr << "WARNING: obscode_lookup failed for observatory code " << obscode << "\n";
      } else {
	// Calculate observer's exact heliocentric position and velocity.
	observer_baryvel01(MJD, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, Earthvel, obspos, obsvel);
	outstream1 << fixed << setprecision(8) << MJD << " ";
	outstream1 << fixed << setprecision(2) << obspos.x << " " << obspos.y << " " << obspos.z << " ";
	outstream1 << fixed << setprecision(8) << obsvel.x << " " << obsvel.y << " " << obsvel.z << " " << RA << " " << Dec << "\n";
      }	
    }
  }
  instream1.close();
  outstream1.close();
  return(0);
}

