// September 24, 2026: split_csv_on_column.cpp
// Given an input CSV file and a column number to split on,
// divide it into subfiles based on changes in the split column.
// For example, if the first ten entries all have 1 in the split
// column and then next thirteen entries all have 87 in the split
// column, an output file called outroot.01 will be written with
// the first ten lines and one called outroot.87 will be written
// with the next thirteen.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: split_csv_on_column -infile infile -splitcol splitcol -outroot outroot\n";
}

int main(int argc, char *argv[])
{
  long i=0;
  string infile, outfile;
  string outroot = "splitout";
  string stest,idstring,oldid;
  long splitcol=16;
  int colct,startpoint,endpoint;
  ifstream instream1;
  ofstream outstream1;
  string lnfromfile;
  int max_id_length=0;
  long filect,filenum;
  int verbose=0;
  
  if(argc<0 && argc<1) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-in" || string(argv[i]) == "-csv") {
      if(i+1 < argc) {
	//There is still something to read;
	infile=argv[++i];
	i++;
      }
      else {
	cerr << "Input CSV file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-splitcol" || string(argv[i]) == "-split" || string(argv[i]) == "-col") {
      if(i+1 < argc) {
	//There is still something to read;
	splitcol=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Split column keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outroot") {
      if(i+1 < argc) {
	//There is still something to read;
	outroot=argv[++i];
	i++;
      }
      else {
	cerr << "Output file root name keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERB" || string(argv[i]) == "-VERBOSE") {
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

  // Check for required arguments
  if(infile.size()<=0) {
    cerr << "ERROR: you have not specified an input file\n";
    show_usage();
    return(1);
  }
 
  // Read the input file. Note that we assume there is no header line.
  // First read is just to find the maximum number of characters in the split column
  instream1.open(infile);
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    getline(instream1,lnfromfile);
    if(lnfromfile.size()>0) {
      startpoint=0;
      endpoint=0;
      colct=0;
      while(startpoint<int(lnfromfile.size())) {
	colct++;
	stest="";
	endpoint = get_csv_string01(lnfromfile, stest, startpoint);
	startpoint=endpoint+1;
	if(colct==splitcol) {
	  if(int(stest.size()) > max_id_length) max_id_length = stest.size();
	}
      }
    }
  }
  instream1.close();
  cout << "Longest identifier in splitcol has " << max_id_length << " characters\n";
  // Second read for actual splitting
  instream1.open(infile);
  idstring=oldid="";
  filect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    getline(instream1,lnfromfile);
    if(lnfromfile.size()>0) {
      startpoint=0;
      endpoint=0;
      colct=0;
      while(startpoint<int(lnfromfile.size())) {
	colct++;
	stest="";
	endpoint = get_csv_string01(lnfromfile, stest, startpoint);
	startpoint=endpoint+1;
	if(colct==splitcol) idstring=stest;
      }
      if(idstring.size()>0 && idstring!=oldid) {
	if(filect>0) {
	  // Close the old file we've been writing.
	  outstream1.close();
	}
	// Construct the name of a new file.
	filect++;
	oldid = idstring;
	outfile = outroot;
	for(i=idstring.size();i<max_id_length;i++) outfile = outfile + "0";
	outfile = outfile + idstring + ".csv";
	outstream1.open(outfile);
	if(verbose>0) cout << "Writing output file called " << outfile << "\n";
	// Write the current line to the new file
	outstream1 << lnfromfile << "\n";
      } else if(idstring.size()>0) {
	// Write another line to the file we already have open.
	outstream1 << lnfromfile << "\n";
      }
    }
  }
  outstream1.close();
  instream1.close();

  filenum = filect;
  cout << "Finished writing a total of " << filenum << " output files.\n";
  
  return(0);
}
