// February 26, 2026: undumbify01a.cpp
// Remove b'' nonsense from a file.
// This is a python-related thing where the quantity string
// gets rendered b'string'. In a csv, one can look for all
// the three-character sequences matching  ,b'  and
// all the two-character sequences matching ',
// In both cases the sequence should be replaced with just a comma.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: undumbify01a infile outfile\n";
}

int main(int argc, char *argv[])
{
  string infile,outfile,lnfromfile,newline,stest;
  long i,j;
  i=j=0;
  
  if(argc!=3) {
    show_usage();
    return(1);
  }
  infile=argv[1];
  outfile=argv[2];
  
  ifstream instream1 {infile};
  ofstream outstream1 {outfile};
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    getline(instream1,lnfromfile);
    newline="";
    for(i=0;i<long(lnfromfile.size());i++) {
      if(lnfromfile[i]!='b' && lnfromfile[i]!='\'') {
	// It's safe to add this character: go ahead.
	newline.push_back(lnfromfile[i]);
      } else if(i+1<long(lnfromfile.size()) && lnfromfile[i]=='b' && lnfromfile[i+1]=='\'') {
	// This is one of the problem-sequences: skip forward by 1
	i++;
      } else if(lnfromfile[i]=='\'') {
	// This is the other problem-sequence: do nothing
      } else {
	// This is not one of the problem-sequences: go ahead
	newline.push_back(lnfromfile[i]);
      }
    }
    if(newline.size()>0) outstream1 << newline << "\n";
  }
  instream1.close();
  outstream1.close();
  return(0);
}

