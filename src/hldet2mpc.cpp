// February 15, 2025: hldet2mpc.cpp
// Converts a file in hldet format to MPC80 format.
// Assumes there is no header.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: hldet2mpc hldet_file mpcfile\n";
}

int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  string hldetfile,mpcfile;
  long i=0;
  long j=0;
  long bandlen=0;
  long status=0;
  long obsnum=0;
  int verbose=0;
  int year,month;
  double day;
  int rahr,ramin;
  int decdeg,decmin;
  double Dec,rasec,decsec;
  string signstring;

  
  if(argc<3) {
    show_usage();
    return(1);
  }
  hldetfile=argv[1];
  mpcfile=argv[2];

  cout.precision(17);  
  cout << "input observation file " << hldetfile << "\n";
  cout << "output file " << mpcfile << "\n";
  
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
  
  ofstream outstream1 {mpcfile};

  for(i=0; i<obsnum; i++) {
    // The temporary ID will just be hldet00
    outstream1 << "     hldet00" << "  C";
    mjd2mpcdate(detvec[i].MJD,year,month,day);
    outstream1 << year << " ";
    if(month<10) outstream1 << "0";
    outstream1 << month << " ";
    day = round(day*1000000.0l);
    day /= 1000000.0l;
    if(day<10.0l) outstream1  << fixed << setprecision(6) << "0";
    outstream1  << fixed << setprecision(6) << day;
    // Convert RA, Dec from decimal degrees to sexagesimal format.
    rahr = int(detvec[i].RA/15.0l);
    ramin = int(detvec[i].RA*4.0l - double(rahr)*60.0l);
    rasec = detvec[i].RA*240.0l - double(rahr)*3600.0l - double(ramin)*60.0l;
    rasec = round(rasec*1000.0l);
    rasec /= 1000.0l;
    if(detvec[i].Dec>=0) {
      signstring="+";
      Dec = detvec[i].Dec;
    } else {
      signstring="-";
      Dec = -detvec[i].Dec;
    }
    decdeg = int(Dec);
    decmin = int(Dec*60.0l - double(decdeg)*60.0l);
    decsec = Dec*3600.0l - double(decdeg)*3600.0l - double(decmin)*60.0l;
    decsec = round(decsec*100.0l);
    decsec /= 100.0l;
    // Write out RA and Dec with appropriate zero-padding.
    if(rahr<10) outstream1 << "0";
    outstream1 << rahr << " ";
    if(ramin<10) outstream1 << "0";
    outstream1 << ramin << " ";
    if(rasec<10.0l) outstream1 << "0";
    outstream1 << fixed << setprecision(3) << rasec << signstring;
    if(decdeg<10) outstream1 << "0";
    outstream1 << decdeg << " ";
    if(decmin<10) outstream1 << "0";
    outstream1 << decmin << " ";
    if(decsec<10.0l) outstream1 << "0";
    outstream1 << fixed << setprecision(2) << decsec << "         ";
    // Write out magnitude and band.
    if(detvec[i].mag!=0.0) outstream1 << fixed << setprecision(1) << detvec[i].mag << " ";
    else outstream1 << "99.9 ";
    if(detvec[i].band[0]!='\0' && detvec[i].band[0]!=' ') {
      outstream1 << detvec[i].band;
      // Determine number of spaces to add after band.
      bandlen = j = 0;
      while(j<MINSTRINGLEN && detvec[i].band[j]!='\0') {
	bandlen++;
	j++;
      }
    } else {
      // Write dummy band x, set band length to 1
      outstream1 << "x";
      bandlen=1;
    }
    for(j=0;j<7-bandlen;j++) outstream1 << " ";
    // Write out obscode
    outstream1 << detvec[i].obscode << "\n";
  }
  outstream1.close();
  
  return(0);
}

