// January 14, 2026
// Given an input file in hldet format, divides it into tracklets and
// calculates the sun-relative ecliptic coordinates and ecliptic angular
// velocities of each tracklet. Produces an output file with the same
// format as mpc80_trkangvel01.cpp

#define NORMPLOTSCALE 900.0
#define INSETPLOTSCALE 50.0

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: hldet_ecliptrk hldet_file outfile\n";
}

int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  string hldetfile;
  string outfile;
  string obscode,oldcode;
  long clusterindex,oldindex;
  long obsct=0;
  int status=0;
  long obsnum=0;
  int verbose=1;
  vector <double> MJDvec;
  vector <double> RAvec;
  vector <double> Decvec;
  double poleRA,poleDec,angvel,pa,crosstrack,alongtrack,arc,timespan,GCR;
  double ecliplon,ecliplat,lonmed,lonmean,lonrms,latmean,MJDmean,oldmjd;
  double sunra,sundec,sunlon,sunlat,ecliplon_relative_to_sun,ecliplon_relative_to_opposition;

  if(argc<3) {
    show_usage();
    return(1);
  }
  hldetfile=argv[1];
  outfile=argv[2];

  cout.precision(17);  
  if(verbose>0) cout << "input observation file " << hldetfile << "\n";
  if(verbose>0) cout << "output file " << outfile << "\n";
  
  // Read input observation file.
  detvec={};
  status=read_hldet_file(hldetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << hldetfile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  obsnum = detvec.size();
  if(verbose>0) cout << "Read " << obsnum << " data lines from observation file " << hldetfile << "\n";
  if(obsnum<=0) {
    cerr << "ERROR: no lines read from input file\n";
    return(1);
  }
  ofstream outstream1 {outfile};
  obsct=0;
  MJDvec = RAvec = Decvec = {};
  clusterindex = detvec[0].index;
  oldindex = clusterindex;
  oldmjd = detvec[0].MJD;
  obscode = string(detvec[0].obscode);
  oldcode = obscode;
  while(obsct<obsnum) {
    clusterindex = detvec[obsct].index;
    obscode = string(detvec[obsct].obscode);
    // Convert RA and Dec to ecliptic coordinates.
    poleswitch02(detvec[obsct].RA, detvec[obsct].Dec, NEPRA, NEPDEC, 90.0l, ecliplon, ecliplat);
    if(verbose>0) cout << "RA, Dec, ecliplon,ecliplat = " << detvec[obsct].RA << " " << detvec[obsct].Dec << " " << ecliplon << " " << ecliplat << "\n";
    if(verbose>0) cout << "clusterindex, oldindex, obscode, oldcode: " << clusterindex << " " << oldindex << " " << obscode << " " << oldcode << " " << detvec[obsct].MJD << " " << oldmjd << "\n";
    if(clusterindex==oldindex && obscode==oldcode && fabs(detvec[obsct].MJD-oldmjd)<NIGHTSTEP) {
      // This line appears to continue a previously started tracklet.
      // Add the new point to the tracklet vectors.
      MJDvec.push_back(detvec[obsct].MJD);
      RAvec.push_back(ecliplon);
      Decvec.push_back(ecliplat);
    } else {
      // We appear to have transitioned to a new tracklet.
      if(verbose>0) cout << "Finishing out a tracklet with " << MJDvec.size() << " points\n";
      if(long(MJDvec.size())>=2) {
	// Analyze the newly completed tracklet.
	arc = timespan = GCR = 0.0;
	status = greatcircfit(MJDvec, RAvec, Decvec, poleRA, poleDec, angvel, pa, crosstrack, alongtrack);
	if(status==0) {
	  // Great Circle fit succeeded. Calculate angular arc, timespan, and GCR.
	  timespan = MJDvec[MJDvec.size()-1] - MJDvec[0];
	  arc = angvel*timespan*3600.0;
	  GCR = sqrt(crosstrack*crosstrack + alongtrack*alongtrack);
	  if(verbose>0) cout << "angvel, pa, timespan, arc, GCR = " << angvel << " " << pa << " " << timespan << " " << arc << " " << GCR << "\n";
	  // Calculate the mean ecliptic coordinates
	  lonmed = lonmean = lonrms = latmean = MJDmean = 0.0;
	  anglevec_meanrms(RAvec, 360.0, &lonmed, &lonmean, &lonrms);
	  latmean = dmean01(Decvec);
	  // Calculate the mean MJD and the solar ecliptic coordinates at that time.
	  MJDmean = dmean01(MJDvec);
	  sunradec_approx01(MJDmean, sunra, sundec);
	  poleswitch02(sunra, sundec, NEPRA, NEPDEC, 90.0l, sunlon, sunlat);
	  // Calculate the object's distance from the sun in ecliptic longitude
	  ecliplon_relative_to_sun = lonmean - sunlon;
	  while(ecliplon_relative_to_sun <= -180.0) ecliplon_relative_to_sun += 360.0;
	  while(ecliplon_relative_to_sun > 180.0) ecliplon_relative_to_sun -= 360.0;
	  // Calculate the object's distance from opposition in ecliptic longitude
	  ecliplon_relative_to_opposition = lonmean - sunlon-180.0;
	  while(ecliplon_relative_to_opposition <= -180.0) ecliplon_relative_to_opposition += 360.0;
	  while(ecliplon_relative_to_opposition > 180.0) ecliplon_relative_to_opposition -= 360.0;

	  if(verbose>0) cout << "MJDmean, sunra, sundec, sunlon, sunlat, sunelong: " << MJDmean << " " << sunra << " " << sundec << " " << sunlon << " " << sunlat << " " << ecliplon_relative_to_sun << "\n";
	  // Write sun-relative coordinates and ecliptic angular velocities to output file
	  outstream1 << fixed << setprecision(6) << ecliplon_relative_to_sun << " " << ecliplon_relative_to_opposition << " " << latmean << " " << angvel*sin(pa/DEGPRAD) << " " << angvel*cos(pa/DEGPRAD) << " " << MJDmean << " " << oldindex << " " << oldcode << "\n";
	}
      }
      // Set up to load the new tracklet.
      MJDvec={};
      RAvec={};
      Decvec={};
      oldindex = clusterindex;
      oldcode = obscode;
      oldmjd = detvec[obsct].MJD;
      // Load the current point as the first point of the new tracklet.
      MJDvec.push_back(detvec[obsct].MJD);
      RAvec.push_back(ecliplon);
      Decvec.push_back(ecliplat);
    }
    obsct++;
  }
  // Handle a final tracklet
  if(long(MJDvec.size())>=2) {
    // Analyze the newly completed tracklet.
    arc = timespan = GCR = 0.0;
    status = greatcircfit(MJDvec, RAvec, Decvec, poleRA, poleDec, angvel, pa, crosstrack, alongtrack);
    if(status==0) {
      // Great Circle fit succeeded. Calculate angular arc, timespan, and GCR.
      timespan = MJDvec[MJDvec.size()-1] - MJDvec[0];
      arc = angvel*timespan*3600.0;
      GCR = sqrt(crosstrack*crosstrack + alongtrack*alongtrack);
      if(verbose>0) cout << "angvel, pa, timespan, arc, GCR = " << angvel << " " << pa << " " << timespan << " " << arc << " " << GCR << "\n";
      // Calculate the mean ecliptic coordinates
      lonmed = lonmean = lonrms = latmean = MJDmean = 0.0;
      anglevec_meanrms(RAvec, 360.0, &lonmed, &lonmean, &lonrms);
      latmean = dmean01(Decvec);
      // Calculate the mean MJD and the solar ecliptic coordinates at that time.
      MJDmean = dmean01(MJDvec);
      sunradec_approx01(MJDmean, sunra, sundec);
      poleswitch02(sunra, sundec, NEPRA, NEPDEC, 90.0l, sunlon, sunlat);
      // Calculate the object's distance from the sun in ecliptic longitude
      ecliplon_relative_to_sun = lonmean - sunlon;
      while(ecliplon_relative_to_sun <= -180.0) ecliplon_relative_to_sun += 360.0;
      while(ecliplon_relative_to_sun > 180.0) ecliplon_relative_to_sun -= 360.0;
      // Calculate the object's distance from opposition in ecliptic longitude
      ecliplon_relative_to_opposition = lonmean - sunlon-180.0;
      while(ecliplon_relative_to_opposition <= -180.0) ecliplon_relative_to_opposition += 360.0;
      while(ecliplon_relative_to_opposition > 180.0) ecliplon_relative_to_opposition -= 360.0;

      if(verbose>0) cout << "MJDmean, sunra, sundec, sunlon, sunlat, sunelong: " << MJDmean << " " << sunra << " " << sundec << " " << sunlon << " " << sunlat << " " << ecliplon_relative_to_sun << "\n";
      // Write sun-relative coordinates and ecliptic angular velocities to output file
      outstream1 << fixed << setprecision(6) << ecliplon_relative_to_sun << " " << ecliplon_relative_to_opposition << " " << latmean << " " << angvel*sin(pa/DEGPRAD) << " " << angvel*cos(pa/DEGPRAD) << " " << MJDmean << " " << oldindex << " " << oldcode << "\n";
    }
  }
  outstream1.close();
  return(0);
}

