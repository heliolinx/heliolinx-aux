// December 11, 2025: hldet_trkan.cpp
// Given an input file in hldet format, assumed to correspond
// to observations of a single object, divides it into tracklets and
// analyzes them.
// Assumes there is no header.

#define NORMPLOTSCALE 900.0
#define INSETPLOTSCALE 50.0

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: hldet2mpc.cpp hldet_file mpcfile\n";
}

int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  string hldetfile;
  string outfile;
  long i=0;
  int status=0;
  long obsnum=0;
  int verbose=0;
  vector <double> mjdvec;
  vector <double> timeintervals;
  vector <double> timediffs;
  double timediff;
  double maxmag;
  double totaltime,rescaletime,mininterval;
  
  if(argc<3) {
    show_usage();
    return(1);
  }
  hldetfile=argv[1];
  outfile=argv[2];

  cout.precision(17);  
  if(verbose>0) cout << "input observation file " << hldetfile << "\n";
  if(verbose>0) cout << "output file " << outfile << "\n";
  
  // Read input observation file.
  detvec={};
  status=read_hldet_file(hldetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << hldetfile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  obsnum = detvec.size();
  if(verbose>0) cout << "Read " << obsnum << " data lines from observation file " << hldetfile << "\n";
  for(i=0;i<obsnum;i++) mjdvec.push_back(detvec[i].MJD);
  sort(mjdvec.begin(),mjdvec.end());
  timediff = mjdvec[1] - mjdvec[0];
  timediffs.push_back(timediff);
  for(i=1;i<obsnum-1;i++) {
    timediff = mjdvec[i] - mjdvec[i-1];
    timeintervals.push_back(timediff);
    if((mjdvec[i+1]-mjdvec[i]) < timediff) timediff = mjdvec[i+1]-mjdvec[i];
    timediffs.push_back(timediff);
  }
  i=obsnum-1;  
  timediff = mjdvec[i] - mjdvec[i-1];
  timediffs.push_back(timediff);
  timeintervals.push_back(timediff);
  maxmag = -LARGERR3;
  for(i=0;i<obsnum;i++) {
    if(detvec[i].mag>maxmag) maxmag = detvec[i].mag;
  }
  sort(timediffs.begin(),timediffs.end());
  cout << maxmag << " " << timediffs[obsnum-1] << "\n";

  mininterval = LARGERR3;
  for(i=0;i<obsnum-1;i++) {
    if(timeintervals[i]<mininterval) mininterval=timeintervals[i];
    cout << "int, min " << timeintervals[i] << " " << mininterval << "\n";
  }
  totaltime=0.0;
  for(i=0;i<obsnum-1;i++) {
    timeintervals[i] = sqrt(timeintervals[i]);
    totaltime+=timeintervals[i];
    cout << "int, total " << timeintervals[i] << " " << totaltime << "\n";
  }
  rescaletime = INSETPLOTSCALE;
  cout << rescaletime << "\n";
  for(i=0;i<obsnum-1;i++) {
    rescaletime += timeintervals[i]*NORMPLOTSCALE/totaltime;
    cout << rescaletime << "\n";
  }
    
  ofstream outstream1 {outfile};
  outstream1 << maxmag << " " << timediffs[obsnum-1] << "\n";
  outstream1.close();
  
  return(0);
}

