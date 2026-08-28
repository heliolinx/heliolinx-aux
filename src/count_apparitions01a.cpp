// August 25, 2026: count_apparitions01a
// Given an input linkage in hldet format, count the number of
// apparitions and the number of tracklets per apparition,
// and print out the numbers of tracklets and the year of
// the first tracklet in each apparation. For example, if
// an object was seen on two tracklets in 2014, one each
// in 2015 and 2022, and three tracklets in 2026, the output
// will be: 4 app: 2+1+1+3 2014,2015,2022,2026.
//
// Note that the input file is NOT expected to have a header line.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define DEBUG 0
#define HNUM 4

#define APPNIGHTSTEP 0.6
#define APPSTEP 240.0

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain observations in heliolinc's hldet format.
static void show_usage()
{
  cerr << "Usage: count_apparitions01a infile.\n";
}

int main(int argc, char *argv[])
{
  ifstream instream1;
  string infile;
  vector <hldet> detvec = {};
  vector <double> obsMJD;
  long obsnum,obsct,ptct,trackletct,appnum,appct,alltrknum;
  obsnum = obsct = ptct = appnum = appct = alltrknum = 0;
  vector <long> trackletnum;
  vector <double> appyear;
  double timediff,day;
  timediff = day = 0.0;
  int verbose, year, month, num_singletons, status;
  verbose = year = month = num_singletons = status = 0;
  
  
  if(argc!=2) {
    show_usage();
    return(1);
  }

  infile=argv[1];
  // Read input observation file.
  detvec={};
  status=read_hldet_file(infile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << infile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  obsnum = detvec.size();
  // cout << "Read " << obsnum << " data lines from observation file " << infile << "\n";

  // Load the MJD, RA, and Dec vectors.
  obsMJD = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD.push_back(detvec[obsct].MJD);
  }
  sort(obsMJD.begin(),obsMJD.end());

  // Perform tracklet/apparition analysis
  trackletnum = {};
  appyear = {};
  // First point is guaranteed to start a new apparation
  obsct=0;
  mjd2mpcdate(obsMJD[obsct], year, month, day);
  appyear.push_back(year);
  ptct=trackletct=alltrknum=1;
  for(obsct=1;obsct<obsnum;obsct++) {
    timediff = obsMJD[obsct] - obsMJD[obsct-1];
    if(timediff>0.0 && timediff<=APPNIGHTSTEP) {
      // This is an additional point in the same tracklet
      ptct++;
    } else if(timediff>APPNIGHTSTEP && timediff<=APPSTEP) {
      // This is the start of a new tracklet in the same apparition.
      trackletct++;
      alltrknum++;
      if(ptct<=1) num_singletons+=1;
      ptct=1;
    } else if(timediff>APPSTEP) {
      // This is the start of a new tracklet in a new apparition.
      trackletnum.push_back(trackletct);
      if(ptct<=1) num_singletons+=1;
      trackletct = ptct = 1;
      alltrknum++;
      mjd2mpcdate(obsMJD[obsct], year, month, day);
      appyear.push_back(year);
    } else {
      cerr << "ERROR: illogical time difference " << timediff << "\n";
      return(1);
    }
  }
  // Wrap up the final tracklet and apparition.
  if(ptct<=1) num_singletons+=1;
  trackletnum.push_back(trackletct);

  appnum = appyear.size();
  if(long(trackletnum.size())!=appnum) {
    cerr << "ERROR: apparition year and tracklet count vector sized do not match!\n";
    cerr << appyear.size() << " vs. " << trackletnum.size() << "\n";
    return(2);
  }
  
  cout << infile << " " << appnum << " app: " << trackletnum[0];
  for(appct=1;appct<appnum;appct++) cout << "+" << trackletnum[appct];
  cout << " " << appyear[0];
  for(appct=1;appct<appnum;appct++) cout << "," << appyear[appct];
  cout << "\n";
  cout << num_singletons << " of " << alltrknum << " total tracklets were singletons\n";
  
  return(0);
}
