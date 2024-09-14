// September 10, 2024: calc_heliovanemat.cpp
// Calculates a hypothesis matrix for heliovane, by interpolating
// an input file giving sampling in each dimension,
// as a function of the difference in heliocentric
// ecliptic longitude between the asteroid and Earth.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define STEPMAX 10000000

static void show_usage()
{
  cerr << "Usage: calc_heliovanemat -interp interpolation file -mind minimum ecliptic longitude difference (deg) -maxd maximum ecliptic longitude difference (deg) -minv minimum velocity in ecliptic longitude -maxv maximum velocity in ecliptic longitude -minsunelong minimum solar elongation -unbound 0 -vinf max_v_inf -distfac factor by which to scale distance sampling -velfac factor by which to scale velocity sampling -accfac factor by which to scale acceleration sampling -accrange fraction of full acceleration range to use -dsign Negative, Positive, or Both -middle usemiddle -outfile output file\n";

}

int main(int argc, char *argv[])
{
  double mindist = 5.0;
  double maxdist = 90.0l;
  double minvel = 0.5;
  double maxvel = LARGERR;
  double diststep = 1.0l;
  double velstep = 1.0e-2;
  double distfac,velfac,accfac;
  distfac = velfac = accfac = 1.0l;
  int unbound = 0;
  long i=0;
  long stepct=0;
  int verbose=0;
  ifstream instream1;
  ofstream outstream1;
  double accelstep,maxacc,vtan;
  accelstep = maxacc = vtan = 0.0l;
  long accelnum,accelct,velnum,velct;
  string interpfile,outfile;
  int usemiddle=0;
  double minsunelong = 45.0;
  double ds,vs,as;
  vector <double> distvec,distsampvec,velsampvec,accsampvec;
  long ilo,ihi;
  string lnfromfile;
  double max_v_inf=0.0;
  double delta_lambda = 0.0;
  double lambda_dot,lambda_ddot;
  lambda_dot = lambda_ddot = 0.0;
  double tanvel = 0.0;
  double rhodot = 0.0;
  double accrange = 1.0;
  string dsign="Both";
 
  if(argc<5) {
    show_usage();
    return(1);
  }
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-interp" || string(argv[i]) == "-interpfile" || string(argv[i]) == "-infile" || string(argv[i]) == "-interpolation" || string(argv[i]) == "--interp" || string(argv[i]) == "--interpfile" || string(argv[i]) == "--interpolation" || string(argv[i]) == "--interpolation_file") {
      if(i+1 < argc) {
	//There is still something to read;
	interpfile = argv[++i];
	i++;
      }
      else {
	cerr << "Interpolation file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mindist" || string(argv[i]) == "-mind" || string(argv[i]) == "-mindistance" || string(argv[i]) == "-minimumdist" || string(argv[i]) == "--minimum_distance" || string(argv[i]) == "--mindist" || string(argv[i]) == "--mind" || string(argv[i]) == "--mindist") {
      if(i+1 < argc) {
	//There is still something to read;
	mindist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input minimum heliocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxdist" || string(argv[i]) == "-maxd" || string(argv[i]) == "-maxdistance" || string(argv[i]) == "-maximumdist" || string(argv[i]) == "--maximum_distance" || string(argv[i]) == "--maxdist" || string(argv[i]) == "--maxd" || string(argv[i]) == "--maxdist") {
      if(i+1 < argc) {
	//There is still something to read;
	maxdist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input maximum heliocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minv" || string(argv[i]) == "-minvel") {
      if(i+1 < argc) {
	//There is still something to read;
	minvel = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input minimum velocity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxv" || string(argv[i]) == "-maxvel") {
      if(i+1 < argc) {
	//There is still something to read;
	maxvel = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input maximum velocity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minsunelong" || string(argv[i]) == "-minelong" || string(argv[i]) == "-minse") {
      if(i+1 < argc) {
	//There is still something to read;
	minsunelong=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum solar elongation keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-unbound" || string(argv[i]) == "--unbound") {
      if(i+1 < argc) {
	// There is still something to read;
	unbound=stoi(argv[++i]);
	i++;
      } else {
	cerr << "Unbound keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-vinf" || string(argv[i]) == "-maxvinf" || string(argv[i]) == "-max_v_inf") {
      if(i+1 < argc) {
	//There is still something to read;
	max_v_inf=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum impact parameter keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-distfac" || string(argv[i]) == "--distfac" ) {
      if(i+1 < argc) {
	//There is still something to read;
	distfac = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Scaling factor for distance sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velfac" || string(argv[i]) == "--velfac") {
      if(i+1 < argc) {
	//There is still something to read;
	velfac=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Scaling factor for velocity sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-accfac" || string(argv[i]) == "-acfac" || string(argv[i]) == "--acfac" || string(argv[i]) == "--accfac" ) {
      if(i+1 < argc) {
	//There is still something to read;
	accfac = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Scaling factor for acceleration sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-accrange") {
      if(i+1 < argc) {
	//There is still something to read;
	accrange = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Acceleration range keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-middle" || string(argv[i]) == "-usemiddle" || string(argv[i]) == "-hitmiddle") {
      if(i+1 < argc) {
	//There is still something to read;
	usemiddle = stoi(argv[++i]); 
	i++;
      }
      else {
	cerr << "Middle keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-dsign") {
      if(i+1 < argc) {
	//There is still something to read;
	dsign=argv[++i];
	i++;
      }
      else {
	cerr << "Delta_lambda sign keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--out" || string(argv[i]) == "--outfile") {
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
    }  else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  cout << "Ecliptic longitudes will be probed from " << mindist << " to " << maxdist << "deg, relative to Earth\n";
  cout << "Output file to be written is called " << outfile << "\n";

  // Catch required parameters if missing
  if(interpfile.size()<=0) {
    cout << "\nERROR: input interpolation file is required\n";
    show_usage();
    return(1);
  } else if(maxdist<=mindist) {
    cout << "\nERROR: maximum distnace must be strictly greater than minimum distance\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cout << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  // Read interpolation file
  instream1.open(interpfile);
  getline(instream1,lnfromfile); // Skips one-line header
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    instream1 >> delta_lambda >> ds >> vs >> as;
    distvec.push_back(delta_lambda);
    distsampvec.push_back(ds);
    velsampvec.push_back(vs);
    accsampvec.push_back(as);
  }
  instream1.close();
  long interpnum = distvec.size();
  if(long(distsampvec.size())!=interpnum || long(velsampvec.size())!=interpnum || long(accsampvec.size())!=interpnum) {
    cerr << "ERROR: vectors read from input file do not all have\n";
    cerr << "the same length: " << interpnum << " " << distsampvec.size() << " " << velsampvec.size() << " " << accsampvec.size() << "\n";
  }
  cout << "Read " << interpnum << " lines from file " << interpfile << "\n";
  for(i=0;i<interpnum;i++) {
    cout << "Interp file: " << distvec[i] << " " << distsampvec[i] << " " << velsampvec[i] << " " << accsampvec[i] << "\n";
  }
  
  outstream1.open(outfile);
  outstream1 << "#delta_lambda(deg) lambda_vel(deg/day) lambda_vel(deg/day^2)\n";
  delta_lambda = mindist;
  while(delta_lambda<=maxdist && stepct<=STEPMAX) {
    if(delta_lambda<distvec[0]) {
      // We're off the lower end of the interpolation file.
      // Just set all the step sizes to the ones corresponding
      // to the smallest heliocentric distance covered by the file.
      diststep = distsampvec[0];
      velstep = velsampvec[0];
      accelstep = accsampvec[0];
    } else if(delta_lambda>distvec[interpnum-1]) {
      // We're off the upper end of the interpolation file.
      // Just set all the step sizes to the ones corresponding
      // to the largest heliocentric distance covered by the file.
      diststep = distsampvec[interpnum-1];
      velstep = velsampvec[interpnum-1];
      accelstep = accsampvec[interpnum-1];
    } else {
      // Find which two points we are interpolating between
      ihi=0;
      while(ihi < interpnum && distvec[ihi]<=delta_lambda) ihi++;
      if(ihi<=0) {
	cerr << "Logically impossible case at delta_lambda = " << delta_lambda << ", distvec[0] = " << distvec[0] << "\n";
	return(1);
      }
      ilo=ihi-1;
      diststep = exp(log(distsampvec[ilo]) + (log(distsampvec[ihi])-log(distsampvec[ilo]))*(log(delta_lambda)-log(distvec[ilo]))/(log(distvec[ihi])-log(distvec[ilo])));
      velstep = exp(log(velsampvec[ilo]) + (log(velsampvec[ihi])-log(velsampvec[ilo]))*(log(delta_lambda)-log(distvec[ilo]))/(log(distvec[ihi])-log(distvec[ilo])));
      accelstep = exp(log(accsampvec[ilo]) + (log(accsampvec[ihi])-log(accsampvec[ilo]))*(log(delta_lambda)-log(distvec[ilo]))/(log(distvec[ihi])-log(distvec[ilo])));
    }
    diststep *= distfac;
    velstep *= velfac;
    accelstep *= accfac;

    // Now dist is the difference in heliocentric ecliptic longitude
    // between Earth and the asteroid. Calculate the corresponding
    // heliocentric distance, assuming the Earth-asteroid-sun phase
    // angle is 90 degrees -- or, if the difference in ecliptic longitude
    // is too large for that to be possible, assuming the asteroid is at
    // the specified minimum solar elongation.
    
    double rdelta_lambda = delta_lambda/DEGPRAD; // Converts from degrees to radians
    double rminsunelong = minsunelong/DEGPRAD; // Converts from degrees to radians
    double rho = 0.0l; // heliocentric distance projected on ecliptic plane
    double obsdist_ecl = 0.0l; // geocentric distance projected on ecliptic plane
 
    // Difference in ecliptic longitude between Earth and the asteroid, in degrees.
    if(M_PI/2.0l-fabs(rdelta_lambda) > rminsunelong) {
      // Do calculation for a sun-target-observer phase angle of 90 degrees,
      // since this is possible.
      rho = AU_KM*cos(rdelta_lambda);
      obsdist_ecl = AU_KM*sin(fabs(rdelta_lambda));
    } else {
      // A phase angle of 90 degrees would imply the object is too close to the sun.
      // Assume instead that it is at the minimum observable solar elongation.
      rho = AU_KM*sin(rminsunelong)/sin(rminsunelong+fabs(rdelta_lambda));
      obsdist_ecl = AU_KM*sin(fabs(rdelta_lambda))/sin(rminsunelong+fabs(rdelta_lambda));
    }

    // Calculate the circular velocity at this distance from the Sun
    double circvel = sqrt(GMSUN_KM3_SEC2/rho);
    double vesc = sqrt(2.0*GMSUN_KM3_SEC2/rho);
    double vmax = vesc;
    if(unbound==1) vmax = sqrt(vesc*vesc + max_v_inf*max_v_inf);
    double maxangvel = vmax/rho*DEGPRAD*SOLARDAY; // velocity in degrees per day
    double tminvel=minvel;
    double tmaxvel=maxvel;
    // Make sure tminvel and tmaxvel are consistent with the limiting
    // physical velocity vmax.
    if(tminvel < -maxangvel) tminvel = -maxangvel;
    if(tmaxvel > maxangvel) tmaxvel = maxangvel;

    if(verbose>=1) cout << "Delta_lambda is " << delta_lambda << " deg, rho = " << rho/AU_KM << " AU, vmax = " << maxangvel << " deg/day, velstep = " << velstep << "deg/day, accelstep = " << accelstep << " deg/day^2\n";
    velnum = ceil((tmaxvel-tminvel)/velstep);
    if(verbose>=1) cout << "tminvel, tmaxvel, velstep, velnum: " << tminvel << " " << tmaxvel << " " << velstep << " " << velnum << "\n";
    if(velnum<=1) {
      // Velocity step is so large, we'll only probe the central velocity
      lambda_dot = 0.5*(tminvel+tmaxvel);
      tanvel = rho*fabs(lambda_dot)/DEGPRAD/SOLARDAY; // km/sec
      if(tanvel<vmax) {
	rhodot = sqrt(vmax*vmax - tanvel*tanvel); // km/sec
	rhodot *= SOLARDAY; // km/day
	maxacc = 2.0*accrange*rhodot*fabs(lambda_dot)/rho; // deg/day^2
	accelnum = ceil(2.0*maxacc/accelstep);
	if(usemiddle==1 && accelnum>1 && accelnum%2==0) {
	  // We are currently set to probe an even number of acceleration
	  // points. This means we won't probe the center of the range
	  // of possible accelerations, which is often the most probable.
	  // Add one additional point, to ensure that the midpoint of
	  // the possible acceleration range gets probed.
	  accelnum+=1;
	}
	if(accelnum<=1) {
	  // Acceleration step is so large, we'll only use the central value,
	  // which for heliovane is always zero.
	  if(dsign!="Negative") outstream1 << delta_lambda << " " << lambda_dot << " " << 0.0 << "\n";
	  if(dsign!="Positive") outstream1 << -delta_lambda << " " << lambda_dot << " " << 0.0 << "\n";
	  if(verbose>=1) cout << "Writing point " << stepct << ": " << delta_lambda << " " << lambda_dot << " " << 0.0 << "\n";
	  stepct++;
	} else {
	  // Customize the acceleration sampling to have approximately the nominal
	  // step size, but optimally sample the full range with the specified number of steps
	  for(accelct=0;accelct<accelnum;accelct++) {
	    lambda_ddot = -maxacc + (double(accelct)+0.5l)*2.0*maxacc/double(accelnum);
	    if(dsign!="Negative") outstream1 << delta_lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	    if(dsign!="Positive") outstream1 << -delta_lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	    if(verbose>=1) cout << "Writing point " << stepct << ": " << delta_lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	    stepct++;
	  }
	}
      } else {
	cerr << "Warning: lambda = " << delta_lambda << " deg, rho = " << rho/AU_KM << " AU, lambda_dot = " << lambda_dot << " deg/day produces min tangential velocity " << tanvel << " km/sec, larger than vmax = " << vmax << "\n";
      }
    } else {
      // We have more than one velocity step. Customize the velocity sampling to have
      // approximately the nominal step size, but to optimally sample the full range
      // with the specified number of steps.
      // Now for the loop, probing everything except the endpoints
      for(velct=0;velct<velnum;velct++) {
	lambda_dot = tminvel + (double(velct)+0.5)*(tmaxvel-tminvel)/double(velnum);
	tanvel = rho*fabs(lambda_dot)/DEGPRAD/SOLARDAY; // km/sec
	if(tanvel<vmax) {
	  rhodot = sqrt(vmax*vmax - tanvel*tanvel); // km/sec
	  rhodot *= SOLARDAY; // km/day
	  maxacc = 2.0*accrange*rhodot*fabs(lambda_dot)/rho; // deg/day^2
	  accelnum = ceil(2.0*maxacc/accelstep);
	  if(usemiddle==1 && accelnum>1 && accelnum%2==0) {
	    // We are currently set to probe an even number of acceleration
	    // points. This means we won't probe the center of the range
	    // of possible accelerations, which is often the most probable.
	    // Add one additional point, to ensure that the midpoint of
	    // the possible acceleration range gets probed.
	    accelnum+=1;
	  }
	  if(accelnum<=1) {
	    // Acceleration step is so large, we'll only use the central value,
	    // which for heliovane is always zero.
	    if(dsign!="Negative") outstream1 << delta_lambda << " " << lambda_dot << " " << 0.0 << "\n";
	    if(dsign!="Positive") outstream1 << -delta_lambda << " " << lambda_dot << " " << 0.0 << "\n";
	    if(verbose>=1) cout << "Writing point " << stepct << ": " << delta_lambda << " " << lambda_dot << " " << 0.0 << "\n";
	    stepct++;
	  } else {
	    // Customize the acceleration sampling to have approximately the nominal
	    // step size, but optimally sample the full range with the specified number of steps
	    for(accelct=0;accelct<accelnum;accelct++) {
	      lambda_ddot = -maxacc + (double(accelct)+0.5l)*2.0*maxacc/double(accelnum);
	      if(dsign!="Negative") outstream1 << delta_lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	      if(dsign!="Positive") outstream1 << -delta_lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	      if(verbose>=1) cout << "Writing point " << stepct << ": " << delta_lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	      stepct++;
	    }
	  }
	} else {
	  cerr << "Warning: lambda = " << delta_lambda << " deg, rho = " << rho/AU_KM << " AU, lambda_dot = " << lambda_dot << " deg/day produces min tangential velocity " << tanvel << " km/sec, larger than vmax = " << vmax << "\n";
	}
	// Close loop over velocity steps
      }
      // Close else-statement checking that we have more than one velocity step
    }
    delta_lambda += diststep;
  }
	
  return(0);
}
