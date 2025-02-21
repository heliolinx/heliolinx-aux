// April 15, 2024: parse_itf01.cpp.
// Read an MPC-80 column formatted file, perform a Great Circle fit to
// each tracklet, and print out the central MJD, RA, Dec, angular velocities,
// number of points, etc.

#include "solarsyst_dyn_geo01.h"
#include "cmath"


static void show_usage()
{
  cerr << "Usage: parse_itf01 infile outfile\n";
}
int main(int argc, char *argv[])
{
  string infile,outfile,lnfromfile,object,band,obscode,oldobject;
  int reachedeof=0;
  long lct=0;
  long pnum=0;
  int status=0;
  ifstream instream1;
  ofstream outstream1;
  double MJD,RA,Dec,mag;
  MJD = RA = Dec = mag = 0.0l;
  vector <double> MJDvec;
  vector <double> RAvec;
  vector <double> Decvec;
  vector <double> magvec;
  double poleRA, poleDec, angvel, pa, crosstrack, alongtrack;
  poleRA = poleDec = angvel = pa = crosstrack = alongtrack = 0.0l;
  double RAvel,Decvel,avgmag, norm, GCR;
  RAvel = Decvel = avgmag = norm = GCR = 0.0;

  if(argc!=3) {
    show_usage();
    return(1);
  }

  infile = argv[1];
  outfile = argv[2];

  instream1.open(infile);
  outstream1.open(outfile);
  outstream1 << "idstring MJD RA Dec angvel PA RAvel Decvel avgmag pnum GCR alongtrack crosstrack\n";
  oldobject="";
  while(reachedeof==0) {
    object="";
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    lct++;
    MJD = RA = Dec = 0.0l;
    status=mpc80_parseline(lnfromfile, object, &MJD, &RA, &Dec, &mag, band, obscode);
    //cout << "Line read with status " << status << ": " << oldobject << " " << object << " " << MJD << " " << RA << " " << obscode << "\n"; 
    if(object!=oldobject && MJDvec.size()>1) {
      // We have just finished reading a new valid tracklet.
      pnum = MJDvec.size();
      // Perform Great Circle fit.
      //cout << "Fitting tracklet for object " << oldobject << "\n";
      greatcircfit(MJDvec, RAvec, Decvec, poleRA, poleDec, angvel, pa, crosstrack, alongtrack);
      // Calculate RA and Dec components of the angular velocity
      RAvel = angvel*sin(pa/DEGPRAD);
      Decvel = angvel*cos(pa/DEGPRAD);
      // Calculate the Great Circle residual
      GCR = sqrt(crosstrack*crosstrack + alongtrack*alongtrack);
      // Calculate the average magnitude.
      avgmag = norm = 0.0l;
      for(int i=0; i<long(magvec.size()); i++) {
	if(magvec[i]>0.0) {
	  avgmag += magvec[i];
	  norm += 1.0l;
	}
      }
      if(norm>0.0l) avgmag/=norm;
      // Write data to the output file.
      outstream1 << oldobject << " ";
      outstream1 << fixed << setprecision(6) << MJDvec[0] << " " << RAvec[0] << " " << Decvec[0] << " ";
      outstream1 << fixed << setprecision(3) << angvel << " " << pa << " " << RAvel << " " << Decvel << " ";
      outstream1 << fixed << setprecision(2) << avgmag << " " << pnum << " ";
      outstream1 << fixed << setprecision(3) << GCR << " " << alongtrack << " " << crosstrack << "\n";
      cout << oldobject << " ";
      cout << fixed << setprecision(6) << MJDvec[0] << " " << RAvec[0] << " " << Decvec[0] << " ";
      cout << fixed << setprecision(3) << angvel << " " << pa << " " << RAvel << " " << Decvel << " ";
      cout << fixed << setprecision(2) << avgmag << " " << pnum << " ";
      cout << fixed << setprecision(3) << GCR << " " << alongtrack << " " << crosstrack << "\n";
      // Wipe the vectors to accommodate the next tracklet
      MJDvec = RAvec = Decvec = magvec = {};
      if(status==0) {
	// Load the first point for the next tracklet
	oldobject = object;
	MJDvec.push_back(MJD);
	RAvec.push_back(RA);
	Decvec.push_back(Dec);
	magvec.push_back(mag);
      }
    } else if((object==oldobject || MJDvec.size()==0) && status==0) {
      // Add a new point to the current tracklet
      oldobject=object;
      MJDvec.push_back(MJD);
      RAvec.push_back(RA);
      Decvec.push_back(Dec);
      magvec.push_back(mag);
    } else if(object!=oldobject && MJDvec.size()==1) {
      cout << "Warning: possible dud tracklet identified with string id " << oldobject << "\n";
      // Wipe the vectors to accommodate the next tracklet
      MJDvec = RAvec = Decvec = magvec = {};
      if(status==0) {
	// Load the first point for the next tracklet
	oldobject = object;
	MJDvec.push_back(MJD);
	RAvec.push_back(RA);
	Decvec.push_back(Dec);
	magvec.push_back(mag);
      }
    }
  }
  instream1.close();
  outstream1.close();
  
  return(0);
}
