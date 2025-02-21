// June 21, 2023: statevec_to_ELvec.cpp
// Read a file with a one-line header, followed by asteroid
// orbits, one per line, with the following fields:
// (1) stringID, (2) MJD at epoch, (3-5) heliocentric X, Y, Z in km,
// and (6-8) heliocentric VX, VY, VZ in km. Write a new file
// of the form required for ingest by, e.g., 4D_KDtree_rangequery01,
// with string ID, energy, and then the three components of the
// vector angular momentuk.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: statevec_to_ELvec -infile input state vector file -outfile output file\n\n";
  cerr << "Note the required format for the input file:\n";
  cerr << "A one-line header, followed by asteroid orbits,\n";
  cerr << "one per line, with the following fields: (1) stringID,\n";
  cerr << "(2) MJD at epoch, (3-5) heliocentric X, Y, Z in km,\n";
  cerr << "and (6-8) heliocentric VX, VY, VZ in km.\n";
}
  
int main(int argc, char *argv[])
{
  string infile,outfile;
  string lnfromfile;
  string desig;
  ifstream instream1;
  ofstream outstream1;
  string epoch;
  int i=0;
  point3d outpos = point3d(0l,0l,0l);
  point3d outvel = point3d(0l,0l,0l);
  point3d Lvec = point3d(0l,0l,0l);
  double E = 0l;
  double heliodist = 0l;
  double Lscalar = 0l;
  
  if(argc<5) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile") {
      if(i+1 < argc) {
	//There is still something to read;
	infile=argv[++i];
	i++;
      }
      else {
	cerr << "Input state vector file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-outfile") {
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

  // Open input file
  instream1.open(infile);
  outstream1.open(outfile);
  // Skip the one-line header
  getline(instream1,lnfromfile);
  // Read the body of the file
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    outpos = point3d(0l,0l,0l);
    outvel = point3d(0l,0l,0l);
    Lvec = point3d(0l,0l,0l);
    E = 0l;
    heliodist = 0l;
    instream1 >> desig >> epoch >> outpos.x >> outpos.y >> outpos.z >> outvel.x >> outvel.y >> outvel.z;
    heliodist = vecabs3d(outpos);
    Lvec = crossprod3d(outpos,outvel);
    Lscalar = vecabs3d(Lvec);
    if(!instream1.fail() && !instream1.bad() && heliodist>0l && Lscalar>0l) {
      E = 0.5l*dotprod3d(outvel, outvel) - GMSUN_KM3_SEC2/heliodist;
      cout << "Energy: " << E << " Lvec: " << Lvec.x << " " << Lvec.y << " " << Lvec.z << " " << Lscalar << "\n";
      outstream1 << fixed << setprecision(6) << desig << " " << E << " ";
      outstream1 << fixed << setprecision(1) << "" << Lvec.x << " " << Lvec.y << " " << Lvec.z << "\n";      
    }
  }
  instream1.close();
  outstream1.close();
  
  return(0);
}
