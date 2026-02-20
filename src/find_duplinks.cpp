// January 29, 2026: find_duplinks
// Reads two files, an 'unlabled' file and a
// 'labeling' file. Constructs a 4-dimensional k-d tree of the
// 'unlabeled' file, which is expected to be larger than the
// 'labeling' file. Performs a range-query of the
// 'unlabeled' k-d tree around each point in the 'labeling'
// file. If there is a match within a given tolerance, the
// nearest neighbor is than labeled with string ID from the
// 'labeling' file, which supercedes the string ID originally
// read from the 'unlabeled' file, if any. Finally, the newly
// labeled 'unlabeled' file is written out. This output file
// is guaranteed to have the same number of lines as the
// input 'unlabeled' file, and the same string ID's except
// where matches were found.

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
  cerr << "Usage: find_duplinks -listfile listfile -colformat column_format_file -matchrad match_in_arcsec -timescale day_to_deg_conversion -uniquefile unique output file -dupfile duplicated output file\n";
}

int main(int argc, char *argv[])
{
  ofstream outstream1;
  ofstream outstream2;
  ifstream instream1;
  ifstream instream2;
  string colformatfile,listfile,stest,idstring;
  vector <string> detfiles;
  string detfile, uniquefile,dupfile,lnfromfile;
  vector <vector <double>> fileMJDvecs;
  vector <vector <double>> fileRAvecs;
  vector <vector <double>> fileDecvecs;
  vector <double> MJDvec;
  vector <double> RAvec;
  vector <double> Decvec;
  vector <long> detnumvec;
  long detnum;
  vector <string> linevec;
  vector <vector <string>> file_linevecs;
  vector <string> IDvec;
  vector <vector <string>> file_IDvecs;
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
  vector <vector <KD_point4d_index>> file_kdtrees;
  vector <long> indexvec;
  vector <double> median_MJDvec;
  double median_MJD;
  double mjdref = 0.0;
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int colformatfile_set=0;
  int filenum=0;
  int filect=0;
  int compfilect=0;
  vector <string> linklines;
  vector <long> linkdets;
  string thisID;
  int isdup;
  string dupID;
  vector <string> svec;
  vector <string> svec2;
  vector <vector <string>> dupIDvec;
  long detct=0;
    
  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-listfile" || string(argv[i]) == "-list" || string(argv[i]) == "-filelist" || string(argv[i]) == "-file_list" || string(argv[i]) == "--listfile" || string(argv[i]) == "--file_list") {
      if(i+1 < argc) {
	//There is still something to read;
	listfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input list file keyword supplied with no corresponding argument\n";
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

  cout << "Input list file is called " << listfile << "\n";
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
  
  // Read input list file
  instream1.open(listfile);
  filenum=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    detfile="";
    instream1 >> detfile;
    if(detfile.size()>0) {
      instream2.open(detfile);
      if(!instream2) {
	cerr << "can't open input detection file " << detfile << "\n";
	return(1);
      }
      detfiles.push_back(detfile);
      filenum++;
      cout << "Successfully opened file " << filenum << ", called " << detfiles[filenum-1] << "\n";
      reachedeof = lct = 0;
      // Skip one-line header
      getline(instream2,lnfromfile);
      lct++;
      MJDvec = Decvec = RAvec = {};
      linevec = IDvec = {};
      while(reachedeof==0) {
	getline(instream2,lnfromfile);
	lct++;
	if(!instream2.eof() && !instream2.fail() && !instream2.bad()) ; // Read on.
	else if(instream2.eof()) reachedeof=1; //End of file, fine.
	else if(instream2.fail()) reachedeof=-1; //Something wrong, warn
	else if(instream2.bad()) reachedeof=-2; //Worse problem, warn
	cout << "Just read line " << lnfromfile << "\n";
	cout << "reachedeof = " << reachedeof << " " << instream2.eof() << " " << instream2.fail() << " " << instream2.bad() << " " << lct << "\n";
	
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
      instream2.close();
      detnum = MJDvec.size();
      if(detnum!=long(RAvec.size()) || detnum!=long(Decvec.size()) || detnum!=long(IDvec.size()) || detnum!=long(linevec.size())) {
	cerr << "ERROR: vector mismatch after reading file " << detfile << " : " << MJDvec.size() << " " << RAvec.size() << " " << Decvec.size() << " " << IDvec.size() << " " << linevec.size() << "\n";
	return(2);
      }
      cout << "Successfully read " << detnum << " lines from input file " << detfile << "\n";
      fileMJDvecs.push_back(MJDvec);
      fileRAvecs.push_back(RAvec);
      fileDecvecs.push_back(Decvec);
      file_linevecs.push_back(linevec);
      file_IDvecs.push_back(IDvec);
      detnumvec.push_back(detnum);
      // Find the median MJD
      median_MJD = dmedian(MJDvec);
      median_MJDvec.push_back(median_MJD);
    }
  }
  instream1.close();
  cout << "Successfully read " << filenum << " input files\n";
  
  // Load dupIDvec with filenum empty string vectors.
  dupIDvec={};
  svec={};
  for(filect=0;filect<filenum;filect++) dupIDvec.push_back(svec);
  
  // Find the master median MJD
  mjdref = dmedian(median_MJDvec);

  // Create k-d trees
  for(filect=0;filect<filenum;filect++) {
    poolvec = {};
    kdvec = {};
    for(detct=0;detct<detnumvec[filect];detct++) {
      onepoint = point4d_index((fileMJDvecs[filect][detct]-mjdref)*timescale/DEGPRAD,cos(fileRAvecs[filect][detct]/DEGPRAD)*cos(fileDecvecs[filect][detct]/DEGPRAD),sin(fileRAvecs[filect][detct]/DEGPRAD)*cos(fileDecvecs[filect][detct]/DEGPRAD), sin(fileDecvecs[filect][detct]/DEGPRAD), detct);
      poolvec.push_back(onepoint);
    }
    cout << "Finished loading pool of " << detnumvec[filect] << " detections for file number " << filect << "\n";
    // Form KDtree
    kdvec={};
    kdroot = splitpoint = 0;
    splitpoint=medind_4d_index(poolvec,1);
    kdpoint = KD_point4d_index(poolvec[splitpoint],-1,-1,1,0);
    kdvec.push_back(kdpoint);
    kdtree_4d_index(poolvec,1,splitpoint,kdroot,kdvec);
    cout << "Finished constructing k-d tree for file " << filect << " with " << kdvec.size() << " items\n";
    file_kdtrees.push_back(kdvec);
  }
  cout << "Finished making all k-d trees\n";
  
  // Loop over files, identifying any duplicates
  outstream1.open(uniquefile);
  outstream2.open(dupfile);
  for(filect=0;filect<filenum;filect++) {
    cout << "Searching for duplicates in file " << filect << " with " << detnumvec[filect] << " detections\n";
    linklines = {};
    linkdets = {};
    thisID = file_IDvecs[filect][0];
    isdup=0;
    for(i=0;i<filenum;i++) dupIDvec[i]={};
    detct=0;
    while(detct<detnumvec[filect]) {
      while(detct<detnumvec[filect] && thisID==file_IDvecs[filect][detct]) {
	if(fileRAvecs[filect][detct] != -999.9) {
	  cout << "Detection " << detct << " appears valid, linklines.size() = " << linklines.size() << "\n";
	  // This point has not already been ouput as a duplicate of something.
	  linklines.push_back(file_linevecs[filect][detct]);
	  linkdets.push_back(detct);
	  querypoint = point4d_index((fileMJDvecs[filect][detct]-mjdref)*timescale/DEGPRAD,cos(fileRAvecs[filect][detct]/DEGPRAD)*cos(fileDecvecs[filect][detct]/DEGPRAD),sin(fileRAvecs[filect][detct]/DEGPRAD)*cos(fileDecvecs[filect][detct]/DEGPRAD), sin(fileDecvecs[filect][detct]/DEGPRAD), detct);
	  for(compfilect=filect+1;compfilect<filenum;compfilect++) {
	    indexvec = {};
	    status = kdrange_4d_index(file_kdtrees[compfilect], querypoint, matchrad/ASECPRAD, indexvec);
	    if(indexvec.size()>0) {
	      // A duplicate was found
	      isdup=1;
	      // Load the stringID of the duplicate to dupIDvec[compfilect]
	      for(i=0;i<long(indexvec.size());i++) {
		// Label the match with an index that maps all the way back to fileMJDvecs, etc
		index = file_kdtrees[compfilect][indexvec[i]].point.index;
		// Load the corresponding string ID.
		dupIDvec[compfilect].push_back(file_IDvecs[compfilect][index]);
	      }
	    }
	  }
	} else cout << "Detection " << detct << " has previously been marked as a duplicate\n";
	detct++;
      }
      if(linklines.size()>0) {
	// We have loaded a linkage.
	if(isdup==0) {
	  // It is a linkage with no duplicates. Write it to the unique output file
	  for(i=0;i<long(linklines.size());i++) outstream1 << linklines[i] << "," << detfiles[filect] << "_" << file_IDvecs[filect][linkdets[i]] << "\n";
	} else {
	  // The linkage has duplicates.
	  // Construct the combined ID for their union
	  dupID = detfiles[filect] + "_" + thisID;
	  for(compfilect=filect+1;compfilect<filenum;compfilect++) {
	    svec =  dupIDvec[compfilect];
	    if(svec.size()>0) {
	      sort(svec.begin(),svec.end());
	      svec2={};
	      svec2.push_back(svec[0]);
	      for(i=1;i<long(svec.size());i++) {
		if(svec[i]!=svec2[svec2.size()-1]) svec2.push_back(svec[i]);
	      }
	      for(i=0;i<long(svec2.size());i++) {
		dupID = dupID + "_" + detfiles[compfilect] + "_" + svec2[i];
	      }
	    }
	  }
	  // Now we have constructed the combined dupID, which is probably quite a monstronsity.
	  // Anyway, write out all of the primary linkage lines, flagged with this dupID
	  for(i=0;i<long(linklines.size());i++) outstream2 << linklines[i] << "," << dupID << "\n";
	  // Now write out all of the duplication lines from each file
	  for(compfilect=filect+1;compfilect<filenum;compfilect++) {
	    svec =  dupIDvec[compfilect];
	    if(svec.size()>0) {
	      sort(svec.begin(),svec.end());
	      svec2={};
	      svec2.push_back(svec[0]);
	      for(i=1;i<long(svec.size());i++) {
		if(svec[i]!=svec2[svec2.size()-1]) svec2.push_back(svec[i]);
	      }
	      for(i=0;i<long(svec2.size());i++) {
		for(j=0;j<detnumvec[compfilect];j++) {
		  if(svec2[i]==file_IDvecs[compfilect][j]) {
		    // This line in file compfilect is a duplicate. Write it out to the duplicate file
		    outstream2 << file_linevecs[compfilect][j] << "," << dupID << "\n";
		    // And mark it as an already-identified duplicate
		    fileRAvecs[compfilect][j] = -999.9;
		  }
		}
	      }
	    }
	  }
	  // End case where we found bona fide duplicates.
	}
      }
      // Set up for the next linkage
      linklines = {};
      linkdets = {};
      if(detct<detnumvec[filect]) thisID = file_IDvecs[filect][detct];
      isdup=0;
      for(i=0;i<filenum;i++) dupIDvec[i]={};
    }
  }
  outstream1.close();
  outstream2.close();

  return(0);
}
