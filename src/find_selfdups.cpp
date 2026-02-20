// January 31, 2026: find_selfdups
// Reads one file, and try to merge overlapping (duplicated)
// linkages within it.

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
  cerr << "Usage: find_selfdups -infile infile -colformat column_format_file -matchrad match_in_arcsec -timescale day_to_deg_conversion -mindupnum mindupnum -uniquefile unique output file -dupfile duplicated output file\n";
}

int main(int argc, char *argv[])
{
  ofstream outstream1;
  ofstream outstream2;
  ifstream instream1;
  string colformatfile,infile,stest,idstring;
  vector <string> detfiles;
  string detfile, uniquefile,dupfile,lnfromfile;
  vector <double> MJDvec;
  vector <double> RAvec;
  vector <double> Decvec;
  vector <long> detnumvec;
  long detnum;
  vector <string> linevec;
  vector <string> IDvec;
  double MJD,RA,Dec;
  double timescale = DAY_TO_DEG_CONV;
  double matchrad = 1.0;
  long i,j,bandlen,status;
  i = j = bandlen = status = 0;
  long kdroot=0;
  long splitpoint=0;
  long index=0;
  long colreadct,lct;
  int reachedeof = 0;
  int startpoint,endpoint;
  point4d_index onepoint = point4d_index(0,0,0,0,0);
  point4d_index querypoint = point4d_index(0,0,0,0,0);
  vector <point4d_index> poolvec;
  KD_point4d_index kdpoint = KD_point4d_index(onepoint,-1,-1,1,0);
  vector <KD_point4d_index> kdvec;
  vector <long> indexvec;
  double mjdref = 0.0;
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int colformatfile_set=0;
  vector <string> linklines;
  vector <long> linkdets;
  string thisID;
  long isdupnum=0;
  string dupID;
  vector <string> svec;
  vector <string> svec2;
  vector <string> dupIDvec;
  long detct=0;
  long mindupnum=1;
    
  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
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
    } else if(string(argv[i]) == "-mindupnum") {
      if(i+1 < argc) {
	//There is still something to read;
	mindupnum=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Timescale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-uniquefile" || string(argv[i]) == "-outunique" || string(argv[i]) == "--outunique_file") {
      if(i+1 < argc) {
	//There is still something to read;
	uniquefile=argv[++i];
	i++;
      }
      else {
	cerr << "Output unique file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-dupfile" || string(argv[i]) == "-outdup" || string(argv[i]) == "--outdup_file") {
      if(i+1 < argc) {
	//There is still something to read;
	dupfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output unique file keyword supplied with no corresponding argument\n";
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
  cout << "Matching radius will be " << matchrad << " arcseconds\n";
  cout << "Timescale for converting 1 day of time to equivalent degrees will be " << timescale << "\n";
  cout << "Output unique file will be called " << uniquefile << "\n";
  cout << "Output duplicate file will be called " << dupfile << "\n";
  
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
  cout << "Successfully opened file " << infile << "\n";
  reachedeof = lct = 0;
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  MJDvec = Decvec = RAvec = {};
  linevec = IDvec = {};
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    cout << "Just read line " << lnfromfile << "\n";
    cout << "reachedeof = " << reachedeof << " " << instream1.eof() << " " << instream1.fail() << " " << instream1.bad() << " " << lct << "\n";
	
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
      MJDvec.push_back(MJD);
      RAvec.push_back(RA);
      Decvec.push_back(Dec);
      IDvec.push_back(idstring);
      linevec.push_back(lnfromfile);
    }
  }
  instream1.close();
  detnum = MJDvec.size();
  if(detnum!=long(RAvec.size()) || detnum!=long(Decvec.size()) || detnum!=long(IDvec.size()) || detnum!=long(linevec.size())) {
    cerr << "ERROR: vector mismatch after reading file " << detfile << " : " << MJDvec.size() << " " << RAvec.size() << " " << Decvec.size() << " " << IDvec.size() << " " << linevec.size() << "\n";
    return(2);
  }
  cout << "Successfully read " << detnum << " lines from input file " << infile << "\n";
  // Find the median MJD
  mjdref = dmedian(MJDvec);
  
  dupIDvec={};
  
  // Create k-d tree
  poolvec = {};
  kdvec = {};
  for(detct=0;detct<detnum;detct++) {
    onepoint = point4d_index((MJDvec[detct]-mjdref)*timescale/DEGPRAD,cos(RAvec[detct]/DEGPRAD)*cos(Decvec[detct]/DEGPRAD),sin(RAvec[detct]/DEGPRAD)*cos(Decvec[detct]/DEGPRAD), sin(Decvec[detct]/DEGPRAD), detct);
    poolvec.push_back(onepoint);
  }
  cout << "Finished loading pool of " << detnum << " detections for file " << infile << "\n";
  // Form KDtree
  kdvec={};
  kdroot = splitpoint = 0;
  splitpoint=medind_4d_index(poolvec,1);
  kdpoint = KD_point4d_index(poolvec[splitpoint],-1,-1,1,0);
  kdvec.push_back(kdpoint);
  kdtree_4d_index(poolvec,1,splitpoint,kdroot,kdvec);
  cout << "Finished constructing k-d tree for file " << infile << " with " << kdvec.size() << " items\n";
  
  // Loop over linkages, identifying any duplicates
  outstream1.open(uniquefile);
  outstream2.open(dupfile);
  linklines = {};
  linkdets = {};
  thisID = IDvec[0];
  isdupnum=0;
  detct=0;
  while(detct<detnum) {
    while(detct<detnum && thisID==IDvec[detct]) {
      if(RAvec[detct] != -999.9) {
	cout << "Detection " << detct << " appears valid, linklines.size() = " << linklines.size() << "\n";
	// This point has not already been ouput as a duplicate of something.
	linklines.push_back(linevec[detct]);
	linkdets.push_back(detct);
	querypoint = point4d_index((MJDvec[detct]-mjdref)*timescale/DEGPRAD,cos(RAvec[detct]/DEGPRAD)*cos(Decvec[detct]/DEGPRAD),sin(RAvec[detct]/DEGPRAD)*cos(Decvec[detct]/DEGPRAD), sin(Decvec[detct]/DEGPRAD), detct);
	indexvec = {};
	status = kdrange_4d_index(kdvec, querypoint, matchrad/ASECPRAD, indexvec);
	for(i=0;i<long(indexvec.size());i++) {
	  index = kdvec[indexvec[i]].point.index;
	  if(IDvec[index] != thisID) {
	    // A duplicate was found
	    isdupnum+=1;
	    // Load the stringID of the duplicate to dupIDvec
	    dupIDvec.push_back(IDvec[index]);
	  }
	}
      } else cout << "Detection " << detct << " has previously been marked as a duplicate\n";
      detct++;
    }
    if(linklines.size()>0) {
      // We have loaded a linkage.
      if(isdupnum<mindupnum) {
	// It is a linkage with no duplicates. Write it to the unique output file
	for(i=0;i<long(linklines.size());i++) outstream1 << linklines[i] << "," << IDvec[linkdets[i]] << "\n";
      } else {
	// The linkage has duplicates.
	// Construct the combined ID for their union
	dupID = thisID;
	svec =  dupIDvec;
	if(svec.size()>0) {
	  sort(svec.begin(),svec.end());
	  svec2={};
	  svec2.push_back(svec[0]);
	  for(i=1;i<long(svec.size());i++) {
	    if(svec[i]!=svec2[svec2.size()-1]) svec2.push_back(svec[i]);
	  }
	  for(i=0;i<long(svec2.size());i++) {
	    dupID = dupID + "_" + svec2[i];
	  }
	}
	// Now we have constructed the combined dupID, which is probably quite a monstronsity.
	// Anyway, write out all of the primary linkage lines, flagged with this dupID
	for(i=0;i<long(linklines.size());i++) outstream2 << linklines[i] << "," << dupID << "\n";
	// Now write out all of the duplication lines
	for(i=0;i<long(svec2.size());i++) {
	  for(j=0;j<detnum;j++) {
	    if(svec2[i]==IDvec[j]) {
	      // This line is a duplicate. Write it out to the duplicate file
	      outstream2 << linevec[j] << "," << dupID << "\n";
	      // And mark it as an already-identified duplicate
	      RAvec[j] = -999.9;
	    }
	  }
	}
	// End case where we found bona fide duplicates.
      }
    }
    // Set up for the next linkage
    linklines = {};
    linkdets = {};
    if(detct<detnum) thisID = IDvec[detct];
    isdupnum=0;
    dupIDvec={};
  }
  outstream1.close();
  outstream2.close();

  return(0);
}
