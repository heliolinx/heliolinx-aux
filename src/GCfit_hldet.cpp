// February 15, 2025: hldet2mpc.cpp
// Converts a file in hldet format to MPC80 format.
// Assumes there is no header.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: hldet2mpc.cpp hldet_file\n";
}

int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  string hldetfile;
  long i=0;
  long status=0;
  int verbose=0;
  long obsnum=0;
  
  if(argc!=2) {
    show_usage();
    return(1);
  }
  hldetfile=argv[1];

  cout.precision(17);  
  cout << "input observation file " << hldetfile << "\n";
  
  // Read input observation file.
  detvec={};
  status=read_hldet_file(hldetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << hldetfile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  obsnum = detvec.size();
  cout << "Read " << obsnum << " data lines from observation file " << hldetfile << "\n";
  
  double poleRA,poleDec,angvel,pa,crosstrack,alongtrack;
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> residuals;
  
  status = greatcircresid(detvec,poleRA,poleDec,angvel,pa,crosstrack,alongtrack,fitRA,fitDec,residuals,verbose);
  double GCRnew=0.0;
  for(i=0;i<long(residuals.size());i++) {
    GCRnew += residuals[i]*residuals[i];
  }
  GCRnew = sqrt(GCRnew/double(residuals.size()));

  cout << "\nMJD RA Dec fitRA fitDec resid\n";
  for(long i=0;i<long(residuals.size());i++) {
    cout << fixed << setprecision(9) << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " " << fitRA[i] << " " << fitDec[i] << " ";
    cout << fixed << setprecision(6) << residuals[i] << "\n";
  }
  cout << "pointnum, angvel(deg/day), angvel(arcsec/hr), PA(degrees)\n";
  cout << fixed << setprecision(6) <<  residuals.size() << " " << angvel << " " << angvel*150.0 << " " << pa << "\n";
  cout << "across(arcsec), along, GCR, GCRnew, arcdist, pole coords (degrees)\n";
  cout << fixed << setprecision(6) << crosstrack << " " << alongtrack << " " << sqrt(crosstrack*crosstrack+alongtrack*alongtrack) << " " << GCRnew << " ";
  cout << fixed << setprecision(9) << poleRA << " " << poleDec << "\n";
  
  return(0);
}

