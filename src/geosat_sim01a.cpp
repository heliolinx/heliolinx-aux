// November 26, 2024: geosat_sim01a.cpp
// Simulates geocentric satellites, to find the range of angular
// velocity and GCR that they might exhibit.
//
// Steps:
// 1. Randomly select a geocentric distance from a log-uniform distribution that goes from 6500 km (200 km altitude) out to 1e6 km.
// 2. Randomly select a time, and a geocentric RA and Dec. Reject if not visible from Vera Rubin Observatory.
// 3. Randomly select a velocity that is less than the geocentric escape velocity.
// 4. Create a simulated tracklet.
// 5. Calculate the GCR and angular velocity of this tracklet.
//
// Problems to be overcome:
// 1. The GCR will depend on the number of points in the tracklet. 3-point tracklets could have low GCR in cases where tracklets with more points would not, because the middle point could happen to fall on an exact line between the first and last points.
//     Solution: for the first trial, simulate many points so we get a good approximation of the GCR of a continuous trajectory (i.e., the limit of infinitely many points).
// 2. The prescription above will produce many retrograde orbits, which are highly improbable.
//     Solution: don't worry about it, just simulate so many that every part of the parameter space is well-probed.
// 3. What integrator should be used? Kepler_fg_func_vec()

#define MAXNORM 1.0
#define MINNORM 0.1
#define MINALT 20.0

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: geosat_sim01a -ranseed random_number_seed -numobs observations-per-tracklet -obslength tracklet-time-span(hr) -mjdstart mjdstart -mjdend mjdend -obslon obslon -plxcos plxcos -plxsin plxsin -mindist mindist -maxdist maxdist -simnum simnum -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  long numobs=0;
  double obslength=0.5;
  int status=0;
  ifstream instream1;
  string stest;
  string outfile;
  double mjdstart = 0.0;
  double mjdend = 0.0;
  double mjdnow = 0.0;
  long simnum=0;
  long simct=0;
  vector <point3d> observerpos;
  vector <point3d> targpos;
  vector <point3d> targvel;
  vector <double> trackletRA;
  vector <double> trackletDec;
  vector <double> trackletMJD;
  point3d temp3vec = point3d(0,0,0);
  point3d obs_to_targ = point3d(0,0,0);
  point3d startpos = point3d(0,0,0);
  point3d startvel = point3d(0,0,0);
  point3d outpos = point3d(0,0,0);
  point3d earthpole = point3d(0,0,0);
  string lnfromfile;
  vector <string> linestringvec;
  int above_horiz=0;
  double obslon = 289.25058l; // Changed to match X05 obscode
  double plxcos = 0.864981l;
  double plxsin = -0.500958l;
  ofstream outstream1;
  string seedstring;
  double mindist =  EARTHEQUATRAD;
  double maxdist =  1.0e6;
  double dist,vel,vesc,logdist,tempnorm,alt,minalt,RA,Dec;
  double inclination=0.0;
  double a,e,E,L;

  if(argc<15) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
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
    } else if(string(argv[i]) == "-numobs") {
      if(i+1 < argc) {
	//There is still something to read;
        numobs=stol(argv[++i]);
	i++;
	if(numobs<3) {
	  cerr << "Error: at least 3 observations required for a meaningful tracklet\n";
	  return(2);
	}
      }
      else {
	cerr << "Observations-per-tracklet keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obslength") {
      if(i+1 < argc) {
	//There is still something to read;
        obslength=stod(argv[++i]);
	i++;
	if(!isnormal(obslength) || obslength<=0.0) {
	  cerr << "Error: invalid tracklet time span (" << obslength << " hr) supplied.\n";
	  cerr << "Timespan must be strictly positive!\n";
	  return(2);
	}
      }
      else {
	cerr << "Observed tracklet timespan keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdstart" || string(argv[i]) == "-ms" || string(argv[i]) == "-MJDstart" || string(argv[i]) == "--mjdstart" || string(argv[i]) == "--MJDstart" || string(argv[i]) == "--ModifiedJulianDaystart" || string(argv[i]) == "--modifiedjuliandaystart") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdstart=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD start keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdend" || string(argv[i]) == "-me" || string(argv[i]) == "-MJDend" || string(argv[i]) == "--mjdend" || string(argv[i]) == "--MJDend" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdend=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obslon") {
      if(i+1 < argc) {
	//There is still something to read;
	obslon=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Observatory longitude keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-plxcos") {
      if(i+1 < argc) {
	//There is still something to read;
	plxcos=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Observatory plxcos keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-plxsin") {
      if(i+1 < argc) {
	//There is still something to read;
	plxsin=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Observatory plxsin keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mindist") {
      if(i+1 < argc) {
	//There is still something to read;
	mindist=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum distance keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxdist") {
      if(i+1 < argc) {
	//There is still something to read;
	maxdist=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Maximum distance keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-simnum" || string(argv[i]) == "-sn" || string(argv[i]) == "-snum" || string(argv[i]) == "--simnum") {
      if(i+1 < argc) {
        //There is still something to read;
        simnum=stol(argv[++i]);
        i++;
      }
      else {
        cerr << "Simulation number keyword supplied with too few corresponding arguments\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outorb" || string(argv[i]) == "--outorbits") {
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

  cout.precision(17);  
  cout << "randum number seed string " << seedstring << "\n";
  cout << "MJD range " << mjdstart << "--" << mjdend << "\n";
  cout << "Observatory parameters: " << obslon << " " << plxcos << " " << plxsin << "\n";
  cout << "geocentric distance range " << mindist << "--" << maxdist << " km\n";
  cout << "number of encounters to simulate: " << simnum << "\n";
  cout << "output file " << outfile << "\n";

  // Calculate the ecliptic vector for Earth's pole
  status = celestial_to_stateunit(0.0, 90.0, earthpole);
  if(status!=0) return(status);
  cout << "Earth's polar vector: " << earthpole.x << " " << earthpole.y << " " << earthpole.z << "\n";
  
  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "simnum dist angvel pa GCR crosstrack alongtrack minalt vescfrac startRA startDec semimaj_axis eccentricity inclination\n";

  // Loop on simulated objects.
  for(simct=0; simct<simnum; simct++) {
    // Randomly assign a time.
    mjdnow = mjdstart + (mjdend-mjdstart)*unitvar(generator);
    // Load the trackletMJD vector
    trackletMJD={};
    observerpos={};
    for(i=0;i<numobs;i++) {
      trackletMJD.push_back(mjdnow+double(i)*obslength/24.0/double(numobs-1));
      status = observer_geocoords01(trackletMJD[i], obslon, plxcos, plxsin, outpos);
      if(status!=0) {
	cerr << "ERROR: observer_geocoords01 failed with status " << status << "\n";
	return(status);
      }
      observerpos.push_back(outpos);
    }
    // Randomly assign a distance
    logdist = log(mindist) + (log(maxdist) - log(mindist))*unitvar(generator);
    dist = exp(logdist);
    // Randomly assign a distance unit vector, making sure it
    // corresponds to a location that is above the horizons
    // at trackletMJD[0].
    above_horiz=0;
    while(above_horiz<=0) {
      tempnorm=0.0;
      while(tempnorm<MINNORM && tempnorm<MAXNORM) {
	temp3vec.x = 2.0*unitvar(generator) - 1.0;
	temp3vec.y = 2.0*unitvar(generator) - 1.0;
	temp3vec.z = 2.0*unitvar(generator) - 1.0;
	tempnorm = vecabs3d(temp3vec);
      }
      startpos.x = temp3vec.x*dist/tempnorm;
      startpos.y = temp3vec.y*dist/tempnorm;
      startpos.z = temp3vec.z*dist/tempnorm;
      obs_to_targ.x = startpos.x - observerpos[0].x;
      obs_to_targ.y = startpos.y - observerpos[0].y;
      obs_to_targ.z = startpos.z - observerpos[0].z;
      // Use the fact that observerpos[0] points to the local
      // zenith to figure out if an object at startpos would
      // be above the horizon (and above the minimum altitude MINALT).
      tempnorm = dotprod3d(observerpos[0], obs_to_targ)/vecabs3d(obs_to_targ)/vecabs3d(observerpos[0]);
      if(tempnorm > cos((90.0-MINALT)/DEGPRAD)) above_horiz=1;
    }
    // cout << "tempnorm = " << tempnorm << ", alt = " << 90.0 - DEGPRAD*acos(tempnorm) << "\n";
    // Randomly assign a velocity
    vesc = sqrt(2.0*GMEARTH_KM3_SEC2/dist);
    vel = vesc*unitvar(generator);
    // Randomly assign a unit vector for the velocity
    tempnorm=0.0;
    while(tempnorm<MINNORM && tempnorm<MAXNORM) {
      temp3vec.x = 2.0*unitvar(generator) - 1.0;
      temp3vec.y = 2.0*unitvar(generator) - 1.0;
      temp3vec.z = 2.0*unitvar(generator) - 1.0;
      tempnorm = vecabs3d(temp3vec);
    }
    startvel.x = temp3vec.x*vel/tempnorm;
    startvel.y = temp3vec.y*vel/tempnorm;
    startvel.z = temp3vec.z*vel/tempnorm;
    temp3vec = crossprod3d(startpos, startvel);
    L=vecabs3d(temp3vec);
    vecnorm3d(temp3vec);
    tempnorm = dotprod3d(temp3vec,earthpole);
    inclination = acos(tempnorm)*DEGPRAD;
    E = 0.5*vel*vel - GMEARTH_KM3_SEC2/dist;
    a = -GMEARTH_KM3_SEC2/2.0/E;
    e = sqrt(1.0 + 2*E*L*L/GMEARTH_KM3_SEC2/GMEARTH_KM3_SEC2);
    
    // Integrate the Keplerian orbit to find the object's position at
    // each time in trackletMJD.
    Kepler_fg_func_vec(GMEARTH_KM3_SEC2, trackletMJD[0], startpos, startvel, trackletMJD, targpos, targvel);
    // Calculate the RA and Dec at each point
    trackletRA = trackletDec = {};
    minalt = 100.0;
    for(i=0;i<numobs;i++) {
      obs_to_targ.x = targpos[i].x - observerpos[i].x;
      obs_to_targ.y = targpos[i].y - observerpos[i].y;
      obs_to_targ.z = targpos[i].z - observerpos[i].z;
      vecnorm3d(obs_to_targ);
      // Calculate the altitude
      tempnorm = dotprod3d(observerpos[i], obs_to_targ)/vecabs3d(observerpos[i]);
      alt = 90.0 - acos(tempnorm)*DEGPRAD;
      if(alt<minalt) minalt=alt;
      // Calculate RA and Dec
      stateunit_to_celestial(obs_to_targ, RA, Dec);
      trackletRA.push_back(RA);
      trackletDec.push_back(Dec);
    }
    double poleRA, poleDec, angvel, pa, crosstrack, alongtrack, GCR;
    status = greatcircfit(trackletMJD, trackletRA, trackletDec, poleRA, poleDec, angvel, pa, crosstrack, alongtrack);
    GCR = sqrt(crosstrack*crosstrack + alongtrack*alongtrack);
    // dist angvel pa GCR crosstrack alongtrack minalt\n";
    outstream1 << simct << " " << dist << " " << angvel << " " << pa << " " << GCR << " " << crosstrack << " " << alongtrack << " " << minalt << " " << vel/vesc << " " << trackletRA[0] << " " << trackletDec[0] << " " << a << " " << e << " " << inclination << "\n";
  }
  outstream1.close();
  
  return(0);
}
