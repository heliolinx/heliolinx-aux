// October 16, 2025: mpc80_to_hldet
// Translate a file that is in MPC 80-column format into a file of hldet format,
// that is, column-delimited file with 16 columns:
// MJD, RA (deg), Dec (deg), mag
// Then the following, which are dummy: trail_len, trail_PA, sigmag, sig_across, sig_along, image
// Object_name, band, obscode
// Then the following, which are all dummy: known_obj, det_qual, index.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: mpc80_to_hldet -infile mpc_file -outfile outfile\n";
}

int main(int argc, char *argv[])
{
  int i=0;
  string mpcfile,outfile;
  string stest;
  ifstream instream1;
  ofstream outstream1;
  string lnfromfile;
  string object, band, obscode;
  double MJD, RA, Dec, mag;
  int status=0;
    
  if(argc<0 && argc<1) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-mpcfile" || string(argv[i]) == "-mpcin" || string(argv[i]) == "-mpcobs" || string(argv[i]) == "--mpcobs" || string(argv[i]) == "-mpcinfile" || string(argv[i]) == "--mpcfile" || string(argv[i]) == "--mpcin" || string(argv[i]) == "--mpcinfile") {
      if(i+1 < argc) {
	//There is still something to read;
	mpcfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input mpc file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--out" || string(argv[i]) == "-o") {
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
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  // Check for required arguments
  if(mpcfile.size()<=0) {
    cerr << "ERROR: you have not specified an input mpc obs file\n";
    show_usage();
    return(1);
  }
  if(outfile.size()<=0) {
    cerr << "ERROR: you have not specified an output file\n";
    show_usage();
    return(1);
  }
  
  // Read file of mpc observations
  instream1.open(mpcfile,ios_base::in);
  outstream1.open(outfile);
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the mpc file
    getline(instream1,lnfromfile);
    if(!instream1.bad() && !instream1.fail() && !instream1.eof() && lnfromfile.size()>=80) {
      status = mpc80_parseline(lnfromfile, object, &MJD, &RA, &Dec, &mag, band, obscode);
      if(status==0) {
	if(band==" ") band="?";
	outstream1 << fixed << setprecision(7) << MJD << "," << RA << "," << Dec << ",";
	outstream1 << fixed << setprecision(4) << mag;
	outstream1 << ",0.0,90.0,0.3,0.5,0.5,0,";
	outstream1 << object << "," << band << "," << obscode << ",0,999,0\n"; 
      }
    }
  }
  instream1.close();
  outstream1.close();
  return(0);
}

