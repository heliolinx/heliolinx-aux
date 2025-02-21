// June 20, 2023
// Read a file with lines in the format of the Minor Planet Center's MPCORB.DAT,
// and do something interesting with them.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: read_mpcorb01 -mpcorb input MPCORB file -outfile output file\n";
}
  
int main(int argc, char *argv[])
{
  string mpcorb,outfile;
  string lnfromfile;
  string desig;
  ifstream instream1;
  ofstream outstream1;
  string epoch;
  int ondata=0;
  double semimaj_axis = 0l; // in AU
  double eccentricity = 0l; // unitless
  double inclination = 0l;  // in degrees
  double long_ascend_node = 0l; // Longitude of the ascending node, in degrees
  double arg_perihelion = 0l;   // Argument of perihelion, in degrees
  double mean_anom = 0l;        // Mean anomaly at the epoch, in degrees
  double mjd_epoch = 0l;        // Epoch for the orbit in MJD
  double mean_daily_motion = 0l; // in degrees/day
  double H = 0l; 
  double G = 0l;
  int i=0;
  int status=0;
  point3d outpos = point3d(0l,0l,0l);
  point3d outvel = point3d(0l,0l,0l);
  point3d Lvec = point3d(0l,0l,0l);
  double E = 0l;
  double heliodist = 0l;
  double Lscalar = 0l;
  double lonLvec = 0l;
  double latLvec = 0l;
  
  asteroid_orbit oneorb = asteroid_orbit(desig,semimaj_axis, eccentricity, inclination, long_ascend_node, arg_perihelion, mean_anom, mjd_epoch, mean_daily_motion,H,G);

  vector <asteroid_orbit> orbvec;
  
  if(argc<5) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-mpcorb") {
      if(i+1 < argc) {
	//There is still something to read;
	mpcorb=argv[++i];
	i++;
      }
      else {
	cerr << "Input MPCORB file keyword supplied with no corresponding argument\n";
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

  // Open MPCORB file
  instream1.open(mpcorb);
  outstream1.open(outfile);
  // Find start of data
  ondata=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad() && !ondata) {
    getline(instream1,lnfromfile);
    ondata=1;
    for(i=0;i<4;i++) if(lnfromfile[i]!='-') ondata=0;
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    status=1;
    if(!instream1.fail() && !instream1.bad() && lnfromfile.size()>103) status = read_orbline(lnfromfile, oneorb);
    if(status==0) {
      cout << fixed << setprecision(6) << oneorb.desig << " " << oneorb.semimaj_axis << " " << oneorb.eccentricity << " " << oneorb.inclination << " " << oneorb.long_ascend_node << " " << oneorb.arg_perihelion << " " << oneorb.mean_anom << " " << oneorb.mjd_epoch << " " << oneorb.mean_daily_motion << "\n";
      Kepler2dyn(oneorb.mjd_epoch, oneorb, outpos,  outvel);
      outpos.x *= AU_KM;
      outpos.y *= AU_KM;
      outpos.z *= AU_KM;
      outvel.x *= AU_KM/SOLARDAY;
      outvel.y *= AU_KM/SOLARDAY;
      outvel.z *= AU_KM/SOLARDAY;
      heliodist = vecabs3d(outpos);
      Lvec = crossprod3d(outpos,outvel);
      Lscalar = vecabs3d(Lvec);
      E = 0.5l*dotprod3d(outvel, outvel) - GMSUN_KM3_SEC2/heliodist;
      //cout << fixed << setprecision(6) << "State vectors: " << outpos.x << " " << outpos.y << " " << outpos.z << " " << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
      //cout << "Energy: " << E << " Lvec: " << Lvec.x << " " << Lvec.y << " " << Lvec.z << " " << Lscalar << "\n";
      E = -GMSUN_KM3_SEC2/2.0l/AU_KM/oneorb.semimaj_axis;
      Lscalar = sqrt((1.0l - DSQUARE(oneorb.eccentricity))*oneorb.semimaj_axis*AU_KM*GMSUN_KM3_SEC2);
      // You subtract 90 degrees from the longitude of the ascending node to get the longitute of the
      // angular momentum vector.
      lonLvec = oneorb.long_ascend_node - 90.0l;
      latLvec = 90.0l - oneorb.inclination;
      Lvec.x = Lscalar*cos(lonLvec/DEGPRAD)*cos(latLvec/DEGPRAD);
      Lvec.y = Lscalar*sin(lonLvec/DEGPRAD)*cos(latLvec/DEGPRAD);
      Lvec.z = Lscalar*sin(latLvec/DEGPRAD);
      cout << "Energy: " << E << " Lvec: " << Lvec.x << " " << Lvec.y << " " << Lvec.z << " " << Lscalar << "\n";
      outstream1 << fixed << setprecision(6) << oneorb.desig << " " << E << " ";
      outstream1 << fixed << setprecision(1) << "" << Lvec.x << " " << Lvec.y << " " << Lvec.z << "\n";      
    }
  }
  instream1.close();
  
  outstream1.close();
  return(0);
}
