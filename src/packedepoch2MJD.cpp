// January 07, 2026: packedepoch2MJD.cpp:
// Given an epoch in packed form from MPCORB, convert it to MJD,
// using the prescription from
// http://www.minorplanetcenter.org/iau/info/PackedDates.html

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: packedepoch2MJD epochstring\n";
}

int main(int argc, char *argv[])
{
  string epochstring;
  double mjdout=0.0;
  
  
  if(argc!=2) {
    show_usage();
    return(1);
  }
  epochstring=argv[1];
  if(epochstring.size()<5) {
    cerr << "Error: epoch string " << epochstring << " is too short to convert to a valid MJD\n";
    return(1);
  }
  mjdout=0;
  int status = packepoch2MJD(epochstring, mjdout);
  if(status!=0) {
    cerr << "ERROR: pack2MJD failed with error status " << status << "\n";
    return(status);
  }
  cout << fixed << setprecision(8) << mjdout << "\n";
  return(0);
}
    
