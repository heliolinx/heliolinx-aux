// April 16, 2026: Kepler_testharness01.cpp:
// Speed-test Kepler functions Kepler_fg_func_int() and
// Kepler_fg_func_int_Ryder_Halley() against one another.
// If the parameter runkep is 0, the program will only generate
// a bunch of random state vectors. If runkep is 1, it will
// read a file of random state vectors and integrate them
// using the older function Kepler_fg_func_int. For runkep 2,
// it will read a file of random state vectors and integrate them
// using the new function Kepler_fg_func_int_Ryder_Halley.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MJDSTART 60650.0
#define MJDEND 60660.0
#define TESTNUM 5000000
#define MAXDIST_AU 3.0

static void show_usage()
{
  cerr << "Usage: Kepler_testharness01 -ranseed random_number_seed -infile infile -outfile outfile -runkep runkep\n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int runkep=0;
  int status=0;
  ifstream instream1;
  string stest;
  string infile;
  string outfile;
  double mjdstart = MJDSTART;
  double mjdend = MJDEND;
  double maxdist = MAXDIST_AU;
  long testnum = TESTNUM;
  ofstream outstream1;
  string seedstring;
  double sundist,vesc,vtot,vresid;
  point3d startpos = point3d(0.0,0.0,0.0);
  point3d startvel = point3d(0.0,0.0,0.0);
  point3d endpos, endvel;
  double x,y,z,vx,vy,vz;
  double MGsun = GMSUN_KM3_SEC2;
    
  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    if(string(argv[i]) == "-rs" || string(argv[i]) == "-ranseed" || string(argv[i]) == "-rseed" || string(argv[i]) == "-seed" || string(argv[i]) == "--ranseed" || string(argv[i]) == "--randomseed" || string(argv[i]) == "--randomnumberseed" || string(argv[i]) == "--random_number_seed") {
      if(i+1 < argc) {
	// There is still something to read;
	seedstring=argv[++i];
	i++;
      } else {
	cerr << "Random number seed keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-infile") {
      if(i+1 < argc) {
	// There is still something to read;
	infile=argv[++i];
	i++;
      } else {
	cerr << "Input file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile") {
      if(i+1 < argc) {
	// There is still something to read;
	outfile=argv[++i];
	i++;
      } else {
	cerr << "output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-runkep") {
      if(i+1 < argc) {
	//There is still something to read;
	runkep=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "runkep keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  if(runkep==0) {
    // Simulate a bunch of bound state vectors, and write them to the output file
    outstream1.open(outfile);
    for(i=0;i<testnum;i++) {
      startpos = point3d(0.0,0.0,0.0);
      startvel = point3d(0.0,0.0,0.0);
      startpos.x = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      startpos.y = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      startpos.z = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      sundist = vecabs3d(startpos);
      vesc = sqrt(2.0*MGsun/sundist);
      startvel.x = vesc - 2.0*vesc*unitvar(generator);
      vtot = vecabs3d(startvel);
      vresid = sqrt(vesc*vesc - vtot*vtot);
      startvel.y = vresid - 2.0*vresid*unitvar(generator);
      vtot = vecabs3d(startvel);
      vresid = sqrt(vesc*vesc - vtot*vtot);
      startvel.z = vresid - 2.0*vresid*unitvar(generator);
      vtot = vecabs3d(startvel);
      while(vtot>=vesc) {
	startvel.x *= 0.9;
	vtot = vecabs3d(startvel);
      }
      outstream1 << fixed << setprecision(6) << startpos.x << " " << startpos.y << " " << startpos.z << " ";
      outstream1 << fixed << setprecision(10) << startvel.x << " " << startvel.y << " " << startvel.z << "\n";
    }
    outstream1.close();
  } else if(runkep==5) {
    // Simulate a bunch of bound state vectors, and write them to the output file
    outstream1.open(outfile);
    for(i=0;i<testnum;i++) {
      startpos = point3d(0.0,0.0,0.0);
      startvel = point3d(0.0,0.0,0.0);
      startpos.x = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      startpos.y = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      startpos.z = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      sundist = vecabs3d(startpos);
      vesc = 0.9*sqrt(2.0*MGsun/sundist);
      startvel.x = vesc - 2.0*vesc*unitvar(generator);
      vtot = vecabs3d(startvel);
      vresid = sqrt(vesc*vesc - vtot*vtot);
      startvel.y = vresid - 2.0*vresid*unitvar(generator);
      vtot = vecabs3d(startvel);
      vresid = sqrt(vesc*vesc - vtot*vtot);
      startvel.z = vresid - 2.0*vresid*unitvar(generator);
      vtot = vecabs3d(startvel);
      while(vtot>=vesc) {
	startvel.x *= 0.9;
	vtot = vecabs3d(startvel);
      }
      outstream1 << fixed << setprecision(6) << startpos.x << " " << startpos.y << " " << startpos.z << " ";
      outstream1 << fixed << setprecision(10) << startvel.x << " " << startvel.y << " " << startvel.z << "\n";
    }
    outstream1.close();
  } else if (runkep==1 || runkep==2) {
    instream1.open(infile);
    outstream1.open(outfile);
    // Read input file, line by line
    while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      instream1 >> x >> y >> z >> vx >> vy >> vz;
      startpos = point3d(x,y,z);
      startvel = point3d(vx,vy,vz);
      if(runkep==1) Kepler_fg_func_int(MGsun, mjdstart, startpos, startvel, mjdend, endpos, endvel);
      else Kepler_fg_func_int_Ryder_Halley(MGsun, mjdstart, startpos, startvel, mjdend, endpos, endvel);
      outstream1 << fixed << setprecision(5) << endpos.x << " " << endpos.y << " " << endpos.z << " ";
      outstream1 << fixed << setprecision(10) << endvel.x << " " << endvel.y << " " << endvel.z << "\n";
    }
  } else if(runkep==3 || runkep==4) {
    for(i=0;i<testnum;i++) {
      startpos = point3d(0.0,0.0,0.0);
      startvel = point3d(0.0,0.0,0.0);
      startpos.x = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      startpos.y = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      startpos.z = (maxdist - 2.0*maxdist*unitvar(generator))*AU_KM;
      sundist = vecabs3d(startpos);
      vesc = sqrt(2.0*MGsun/sundist);
      startvel.x = vesc - 2.0*vesc*unitvar(generator);
      vtot = vecabs3d(startvel);
      vresid = sqrt(vesc*vesc - vtot*vtot);
      startvel.y = vresid - 2.0*vresid*unitvar(generator);
      vtot = vecabs3d(startvel);
      vresid = sqrt(vesc*vesc - vtot*vtot);
      startvel.z = vresid - 2.0*vresid*unitvar(generator);
      vtot = vecabs3d(startvel);
      while(vtot>=vesc) {
	startvel.x *= 0.9;
	vtot = vecabs3d(startvel);
      }
      if(runkep==3) Kepler_fg_func_int(MGsun, mjdstart, startpos, startvel, mjdend, endpos, endvel);
      else Kepler_fg_func_int_Ryder_Halley(MGsun, mjdstart, startpos, startvel, mjdend, endpos, endvel);
    }
  }
  return(0);
}
