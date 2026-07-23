// April 27, 2026: galeclipelong01a.cpp
// Given an input file with MJD, RA, and Dec,
// calculate the galactic longitude and latitude,
// the ecliptic longitude and latitude, and the
// solar elongation.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: galeclipelong01a infile\n";
}
    
int main(int argc, char *argv[])
{
  string infile;
  double MJD,RA,Dec;
  double ecliplon,ecliplat,gal_lon,gal_lat,sunra,sundec,sunelong;
  long linect=0;
  
  
  if(argc!=2) {
    show_usage();
    return(1);
  }
  infile = argv[1];
  // Read input detection file
  ifstream instream1 {infile};
  cout << "MJD RA Dec gal_lon gal_lat ecliplon ecliplat sunelong\n";
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the mpc file
    MJD=RA=Dec=-999.9;
    instream1 >> MJD >> RA >> Dec;
    linect++;
    if(MJD!=-999.9 && RA!=-999.9 && Dec!=-999.9) {
      // Convert RA and Dec to ecliptic coordinates.
      poleswitch02(RA, Dec, NEPRA, NEPDEC, 90.0l, ecliplon, ecliplat);
      // Convert RA and Dec to Galactic coordinates
      poleswitch02(RA, Dec, NGPRA, NGPDEC, NCPGAL_LON, gal_lon, gal_lat);
      sunradec_approx01(MJD, sunra, sundec);
      sunelong = distradec01(RA, Dec, sunra, sundec);
      cout << fixed << setprecision(6) << MJD << " " << RA << " " << Dec << " " << gal_lon << " " << gal_lat << " " << ecliplon << " " << ecliplat << " " << sunelong << "\n";
    }
    else if(!instream1.eof()) cout << "Warning: line " << linect << " not successfully read\n";
  }
  instream1.close();

  return(0);
}
