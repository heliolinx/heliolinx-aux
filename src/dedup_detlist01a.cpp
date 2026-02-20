// January 30, 2026: dedup_detlist01a.cpp
// Given an input file and column formatting enabling the reading
// of MJD, RA, and Dec, produce a new output file with all duplicates
// eliminated. Warning: the code uses an n^2 search, which means the
// input files don't have to be sorted, on the one hand, but on the
// other hand, it will be very slow on large files.

#define DAY_TO_DEG_CONV 24.0 // Treats one day as equivalent to 24 degrees:
                             // in other words, one second of time is
                             // equivalent to one arcsecond on the sky.
#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define COLS_TO_READ 4

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: dedup_detlist01a -infile infile -colformat column_format_file -matchrad match_in_arcsec -timescale day_to_deg_conversion -outfile outfile\n";
}

int main(int argc, char *argv[])
{
  ofstream outstream1;
  ifstream instream1;
  string colformatfile,idstring,stest;
  string infile,outfile,lnfromfile;
  vector <double> MJDvec;
  vector <double> RAvec;
  vector <double> Decvec;
  vector <long> linesizevec;
  long detnum;
  vector <string> linevec;
  vector <string> IDvec;
  double MJD,RA,Dec;
  double timescale = DAY_TO_DEG_CONV;
  double matchrad = 1.0;
  double dist=0.0;
  double timediff=0.0;
  long i,j,bandlen,status,startpoint,endpoint,reachedeof,lct;
  i = j = bandlen = status = startpoint = endpoint = reachedeof = lct = 0;
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int colformatfile_set=0;
  long detct=0;
  long colreadct=0;
      
  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-in") {
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
    } else if(string(argv[i]) == "-matchrad" || string(argv[i]) == "-match_rad" || string(argv[i]) == "-mrad" || string(argv[i]) == "-matching_radius" || string(argv[i]) == "--matching_radius" || string(argv[i]) == "-searchrad") {
      if(i+1 < argc) {
	//There is still something to read;
	matchrad=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Matching radius keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--output_file") {
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
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }
  cout.precision(17);  

  cout << "Input file is called " << infile << "\n";
  if(colformatfile_set==0) {
    cout << "WARNING: No column formatting file supplied\n";
  } else cout << "Column formatting file for is " << colformatfile << "\n";
 
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
  cout << "Successfully opened " << infile << "\n";
  reachedeof = lct = 0;
  // We assume there is no header
  MJDvec = Decvec = RAvec = {};
  linevec = IDvec = {};
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
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
      } 
      startpoint = endpoint+1;
    }
    if(reachedeof == 0 && lnfromfile.size()>=15) {
      if(MJD==0.0) {
	cerr << "ERROR: MJD not read from line " << lct << " of input file " << infile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(RA==-999.9) {
	cerr << "ERROR: RA not read from line " << lct << " of input file " << infile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(Dec==-999.9) {
	cerr << "ERROR: Dec not read from line " << lct << " of input file " << infile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(idstring.size()<=0) {
	cerr << "ERROR: Object ID not read from line " << lct << " of input file " << infile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      MJDvec.push_back(MJD);
      RAvec.push_back(RA);
      Decvec.push_back(Dec);
      IDvec.push_back(idstring);
      linevec.push_back(lnfromfile);
      linesizevec.push_back(long(lnfromfile.size()));
    }
  }
  instream1.close();
  
  detnum = MJDvec.size();
  if(detnum!=long(RAvec.size()) || detnum!=long(Decvec.size()) || detnum!=long(IDvec.size()) || detnum!=long(linevec.size()) || detnum!=long(linesizevec.size())) {
    cerr << "ERROR: vector mismatch after reading file " << infile << " : " << MJDvec.size() << " " << RAvec.size() << " " << Decvec.size() << " " << IDvec.size() << " " << linevec.size() << " " << linesizevec.size() << "\n";
    return(2);
  }
  cout << "Successfully read " << detnum << " lines from input file " << infile << "\n";
  outstream1.open(outfile);
  for(detct=0;detct<detnum;detct++) {
    long bestdup=detct;
    long bestsize = linesizevec[detct];
    //cout << "Working on line " << linevec[detct] << "\n";
    if(bestsize>0) {
      // Line number detct has not already been marked as a duplicate.
      for(i=detct+1;i<detnum;i++) {
	dist = 3600.0*distradec01(RAvec[detct],Decvec[detct],RAvec[i],Decvec[i]);
	timediff = 3600.0*timescale*fabs(MJDvec[detct]-MJDvec[i]);
	dist = sqrt(dist*dist+timediff*timediff);
	if(dist<matchrad) {
	  // Line i is a duplicate of line detct.
	  if(linesizevec[i]>bestsize) {
	    // Prefer i over the current best duplicate because it seems to have more information.
	    // Wipe the current best duplicate
	    linesizevec[bestdup]=0;
	    // Load the new duplicate as the best.
	    bestsize = linesizevec[i];
	    bestdup = i;
	  } else {
	    // Line i is not as good as line detct. Wipe line i.
	    linesizevec[i]=0;
	  }
	}
      }
      // Write out the best line identified to the output file.
      // Note that if the line detct is unique (no duplicates), detct
      // itself will be this best line.
      //cout << "Writing line " << linevec[bestdup] << "\n";
      outstream1 << linevec[bestdup] << "\n";
      // Wipe this line so it does not get written again later.
      linesizevec[bestdup]=0;
    } //else cout << "This line has already been marked as a duplicate\n";
  }
  outstream1.close();
  return(0);
}
