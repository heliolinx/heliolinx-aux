// February 05, 2025: arctrace01d.cpp
// Given an input file giving observer positions and measured RA and Dec
// for an object at various times, given as MJD, fit the object's orbit.
// Start by identifying a short sequence of observations that can be
// reasonably fit with a Keplerian orbit, and performing this fit using
// the Method of Herget. Then, using the approximate state vectors output
// from the Keplerian fit to perform a full, 6-dimensional orbit fit
// including planetary perturbations. This version uses code analogous
// to 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define MINCOLS 9

// Note: configfile contains the masses and ephemerides for all of the planets.
// The observation file must contain
// MJD, observer X, Y, Z, VX, VY, VZ, RA, and Dec. The observer coordinates and velocities must
// be in km and km/sec, relative to the Sun. The RA and Dec must be in decimal degrees.
static void show_usage()
{
  cerr << "Usage: arctrace01b -cfg configfile -observations obsfile -kepspan time_span_for_Keplerian_fit(day) -minchi min_chi_change -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  ifstream instream1;
  string stest;
  string configfile;
  string obsfile;
  string planetfile;
  string outfile;
  double kepspan=15.0;
  long double x,y,z,vx,vy,vz;
  x=y=z=vx=vy=vz=0.0L;
  long double ldval=0.0L;
  double dval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> planetpos_reverse;
  vector <point3LD> Earthpos;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  vector <point3LD> observerpos;
  vector <point3LD> observervel;
  vector <point3LD> orbit05pos;
  vector <point3LD> orbit05vel;
  vector <long double> orbit05MJD;
  vector <point3d> Kepobserverpos;
  vector <long double> obsMJD;
  vector <double> obsRA;
  vector <double> obsDec;
  vector <double> KepMJD;
  vector <double> KepRA;
  vector <double> KepDec;
  vector <double> sigastrom;
  string lnfromfile;
  vector <string> linestringvec;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  point3LD tppos = point3LD(0,0,0);
  point3LD tpvel = point3LD(0,0,0);
  point3d dpos = point3d(0,0,0);
  long double outRA,outDec;
  int fieldnum=0;
  vector <long double> scalestate;
  make_LDvec(60,scalestate);
  vector <long double> newstate;
  make_LDvec(60,newstate);
  vector <long double> chivec;
  make_LDvec(70,chivec);
  vector <long double> beststatevec;
  long double chisq=0L;
  long double newchi=0L;
  long double bestchi=0L;
  vector <vector <long double>> simplex;
  make_LDmat(70,60,simplex);
  vector <double> fitRA;
  vector <double> fitDec;
  vector <double> fitresid;
  vector <double> bestRA;
  vector <double> bestDec;
  vector <double> bestresid;
  long double pvtimescale=PV_TIMESCALE;
  double geodist1 = 1.0;
  double geodist2 = 1.1;
  double ftol = FTOL_HERGET_SIMPLEX;
  double simplex_scale = SIMPLEX_SCALEFAC;
  int obsnum,obsct,i,j,status,configread,polyorder,planetnum,planetct,pctSun;
  int pctEarth,badread,reachedeof,worstpoint,bestpoint,countworse,point1,verbose;
  int simptype,point2,kepnum,kepmax;
  obsnum=obsct=i=j=status=configread=polyorder=planetnum=planetct=pctSun=0;
  pctEarth=badread=reachedeof=worstpoint=bestpoint=countworse=point1=verbose=0;
  simptype=point2=kepnum=kepmax=1;
  vector <double> orbit;
  vector <long double> gradient;
  vector <long double> gradunit;
  long double pos_step = POS_STEP;
  long double beststep = 0.0l;
  long double poschange = 0.0l;
  int itct,chifailct,shortstepct;
  itct=chifailct=shortstepct=0;
  long double minchichange = 0.001l;
  long double chichange = 0.001l;
  vector <vector <long double>> H;
  vector <vector <long double>> OP1;
  vector <vector <long double>> OP2;
  vector <vector <long double>> OP3;
  make_LDmat(6, 6, H);
  vector <long double> oldgrad;
  vector <long double> svec;
  make_LDvec(6,svec);
  vector <long double> oldsvec;
  vector <long double> changevec;
  vector <long double> statevecdiff;
  vector <long double> gradientdiff;
  vector <long double> Hgrad;
  vector <long double> u;
  vector <long double> u2;
  vector <long double> thisway;
  int goodstep=0;
  long double dotprod_gradstate=0.0l;
  long double dotprod_gradHgrad=0.0l;
  int planetfile_startpoint=0;
  int planetfile_refpoint=0;
  int planetfile_endpoint=0;
  
  if(argc<9) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-c" || string(argv[i]) == "-cfg" || string(argv[i]) == "-config" || string(argv[i]) == "--config" || string(argv[i]) == "--configfile" || string(argv[i]) == "--configuration" || string(argv[i]) == "--ConfigFile") {
      if(i+1 < argc) {
	//There is still something to read;
	configfile=argv[++i];
	i++;
	// Read the configuration file. This must happen here, so that
	// default values supplied in this file can be overwritten later
	// if the user desires.
	// Read configuration file.
	ifstream instream1 {configfile};
	if(!instream1) {
	  cerr << "ERROR: can't open input config file " << configfile << "\n";
	  return(1);
	}
	// Read polyorder
	status=readconfigint(instream1,&polyorder);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&polyorder);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Polynomial order for orbit integration read as " << polyorder << "\n";
	// Read the number of planets.
	status=readconfigint(instream1,&planetnum);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&planetnum);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Number of planets read as " << planetnum << "\n";
	// Read the index of the Sun within the planet vectors
	status=readconfigint(instream1,&pctSun);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctSun);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "The Sun is planet number " << pctSun << "\n";
	// Read the index of Earth within the planet vectors
	status=readconfigint(instream1,&pctEarth);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctEarth);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Earth is planet number " << pctEarth << "\n";
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the mass for planet number planetct
	  status=readconfigLD(instream1,&ldval);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigLD(instream1,&ldval);
	  }
	  planetmasses.push_back(ldval);
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
	}
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the ephemeris file for planet number planetct
	  status=readconfigstring(instream1,planetfile);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigstring(instream1,planetfile);
	  }
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "Ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
	  mjdtest={};
	  temppos={};
	  tempvel={};
	  read_horizons_fileLD(planetfile,mjdtest,temppos,tempvel);
	  if(planetct==0) planetmjd=mjdtest;
	  else {
	    for(j=0;j<long(planetmjd.size());j++) {
	      if(mjdtest[j]!=planetmjd[j]) {
		cout << "ERROR: time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(1);
	      }
	    }
	  }
	  for(j=0;j<long(temppos.size());j++) {
	    planetpos.push_back(temppos[j]);
	  }
	  if(planetct == pctEarth) Earthpos = temppos;
	  if(planetct == pctSun) {
	    Sunpos = temppos;
	    Sunvel = tempvel;
	  }
	  cout << "Finished reading ephemeris file " << planetfile << "\n";
	}
	// Close input stream that was reading the config file.
	instream1.close();
	configread=1;
      } else {
	cerr << "Configuration file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }
    // Stop execution if we haven't read a config file successfully by now.
    if(configread!=1) {
      cerr << "ERROR: configuration file must be supplied and successfully";
      cerr << "read before any of the other parameters. This is because\n";
      cerr << "it contains defaults that can optionally be overridden by\n";
      cerr << "user-supplied values in later arguments.\n";
      return(2);
    } else cout << "Configuration file read successfully\n";
    if(string(argv[i]) == "-o" || string(argv[i]) == "-obs" || string(argv[i]) == "-obsfile" || string(argv[i]) == "-observations" || string(argv[i]) == "--obsfile" || string(argv[i]) == "--observation" || string(argv[i]) == "--observationfile" || string(argv[i]) == "--observations") {
      if(i+1 < argc) {
	// There is still something to read;
	obsfile=argv[++i];
	i++;
      } else {
	cerr << "Observation file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-kepspan") {
      if(i+1 < argc) {
	//There is still something to read;
	kepspan=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Keplerian time-span keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minchi") {
      if(i+1 < argc) {
	//There is still something to read;
	minchichange=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum chi change keyword supplied with no corresponding argument\n";
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
  cout << "input configuration file " << configfile << "\n";
  cout << "input observation file " << obsfile << "\n";
  cout << "input timespan for Keplerian fit " << kepspan << "\n";
  cout << "output file " << outfile << "\n";

  // Read input observation file.
  instream1.open(obsfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << obsfile << "\n";
    return(1);
  }
  badread=0;
  while(reachedeof==0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(lnfromfile.size()>20 && (reachedeof==0 || reachedeof==1)) {
      fieldnum = get_col_vector01(lnfromfile, linestringvec);
      if(fieldnum>=MINCOLS) {
	ldval = stold(linestringvec[0]);
	obsMJD.push_back(ldval);
	x = stold(linestringvec[1]);
	y = stold(linestringvec[2]);
	z = stold(linestringvec[3]);
	vx = stold(linestringvec[4]);
	vy = stold(linestringvec[5]);
	vz = stold(linestringvec[6]);
	outpos = point3LD(x,y,z);
	observerpos.push_back(outpos);
	outpos = point3LD(vx,vy,vz);
	observervel.push_back(outpos);
	dval = stod(linestringvec[7]);
	obsRA.push_back(dval);
	dval = stod(linestringvec[8]);
	obsDec.push_back(dval);
      } else badread=1;
    }
    cout << obsMJD.size() << " " << obsMJD[obsMJD.size()-1] << " " << obsRA[obsMJD.size()-1] << " " << obsDec[obsMJD.size()-1] << "\n";
  }
  obsnum=obsMJD.size();
  if(obsnum<2 || long(obsRA.size())!=obsnum || long(obsDec.size())!=obsnum || long(observerpos.size())!=obsnum || long(observervel.size())!=obsnum) {
    cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << obsnum << " " << obsRA.size() << " " << obsDec.size() << " " << observerpos.size() << " " << observervel.size() << "\n";
    return(1);
  }
  if(badread>=1){
    cerr << "Error reading " << obsfile << ": apparent short line\n";
    return(1);
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();

  // Calculate a timescale for converting velocity to position units.
  pvtimescale = SOLARDAY*(obsMJD[obsnum-1]-obsMJD[0])/2.0l;

  // Find the time period of length kepspan
  // with the largest number of observations
  bestpoint=kepnum=kepmax=0;
  for(i=0;i<obsnum;i++) {
    j=i;
    kepnum=0;
    while(j<obsnum && obsMJD[j]-obsMJD[i] < kepspan) {
      kepnum++;
      j++;
    }
    if(kepnum>kepmax) {
      kepmax=kepnum;
      bestpoint=i;
    }
  }
  kepnum=kepmax;
  cout << "Best sequence contains " << kepnum << " points, and runs from obsMJD[" << bestpoint << "]=" << obsMJD[bestpoint] << " to obsMJD[" << bestpoint+kepnum-1 << "]=" << obsMJD[bestpoint+kepnum-1] << "\n";
  KepMJD = KepRA = KepDec = sigastrom = {};
  Kepobserverpos = {};
  for(i=bestpoint;i<bestpoint+kepnum;i++) {
    KepMJD.push_back(double(obsMJD[i]));
    dpos = point3d(observerpos[i].x,observerpos[i].y,observerpos[i].z);
    Kepobserverpos.push_back(dpos);
    KepRA.push_back(obsRA[i]);
    KepDec.push_back(obsDec[i]);
    sigastrom.push_back(1.0);
  }
  fitDec = fitRA = fitresid = orbit = {};
  chisq = Hergetfit_vstar(geodist1, geodist2, simplex_scale, simptype, ftol, 1, kepnum, Kepobserverpos, KepMJD, KepRA, KepDec, sigastrom, fitRA, fitDec, fitresid, orbit, verbose);
  cout << "Keplerian fit produced chisq = " << chisq << "\n";

  planetfile_startpoint = planetfile_refpoint = planetfile_endpoint = -99;
  j=0;
  while(j<long(planetmjd.size()) && planetmjd[j]<obsMJD[0]) j++;
  planetfile_startpoint = j-1; // This is the last point in planetmjd that is before the start of the observations
  j=planetfile_startpoint;
  while(j<long(planetmjd.size()) && planetmjd[j]<0.5*(KepMJD[0]+KepMJD[kepnum-1])) j++;
  planetfile_refpoint = j-1;
  j=planetfile_startpoint;
  while(j<long(planetmjd.size()) && planetmjd[j]<obsMJD[obsnum-1]) j++;
  planetfile_endpoint = j; // This is the first point in planetmjd that is after the end of the observations

  planetfile_startpoint -= polyorder+1; // Pad starting and ending points so interpolation works.
  planetfile_endpoint += polyorder+1;
  
  if(planetfile_startpoint<0 || planetfile_endpoint<0 || planetfile_refpoint<0) {
    cerr << "ERROR: observation start " << obsMJD[0] << " and/or end " << obsMJD[obsnum-1] << "and or reference time" << 0.5*(KepMJD[0]+KepMJD[kepnum-1]) << " could not be matched\nto any timestep in the planet files\n";
    return(1);
  } else {
    cout << "observation start " << obsMJD[0] << " corresponds to timestep " << planetfile_startpoint << " in the planet files\n";
    cout << "reference time " << 0.5*(KepMJD[0]+KepMJD[kepnum-1]) << " corresponds to timestep " << planetfile_refpoint << " in the planet files\n";
    cout << "observation end " << obsMJD[obsnum-1] << " corresponds to timestep " << planetfile_endpoint << " in the planet files\n";
  }
  

  // Integrate the Keplerian state vectors from their current reference MJD
  // to the MJD corresponding to planetfile_refpoint
  point3d dstartpos=point3d(orbit[3],orbit[4],orbit[5]);
  point3d dstartvel=point3d(orbit[6],orbit[7],orbit[8]);
  point3d dendpos=point3d(0,0,0);
  point3d dendvel=point3d(0,0,0);
  
  Kepler_univ_int(GMSUN_KM3_SEC2, orbit[2], dstartpos, dstartvel, planetmjd[planetfile_refpoint], dendpos, dendvel);
  
  // Convert Keplerian state vectors from heliocentric to barycentric coords.
  startpos.x = dendpos.x + Sunpos[planetfile_refpoint].x;
  startpos.y = dendpos.y + Sunpos[planetfile_refpoint].y;
  startpos.z = dendpos.z + Sunpos[planetfile_refpoint].z;
  startvel.x = dendvel.x + Sunvel[planetfile_refpoint].x;
  startvel.y = dendvel.y + Sunvel[planetfile_refpoint].y;
  startvel.z = dendvel.z + Sunvel[planetfile_refpoint].z;

  // Convert observer positions from heliocentric to barycentric coords.
  sigastrom = {};
  for(i=0;i<obsnum;i++) {
    planetposvel01LD(obsMJD[i], polyorder, planetmjd, Sunpos, Sunvel, outpos, outvel);
    observerpos[i].x += outpos.x;
    observerpos[i].y += outpos.y;
    observerpos[i].z += outpos.z;
    observervel[i].x += outvel.x;
    observervel[i].y += outvel.y;
    observervel[i].z += outvel.z;
    sigastrom.push_back(1.0);
  }

  // Integrate the orbit to the planetfile timestamps
  integrate_orbit05LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, startpos, startvel, planetfile_startpoint, planetfile_refpoint, planetfile_endpoint, orbit05MJD, orbit05pos, orbit05vel);
  for(i=0;i<long(orbit05MJD.size());i++) {
    cout << orbit05MJD[i] << " " << planetmjd[planetfile_refpoint] << " " << orbit05pos[i].x << " " << orbit05pos[i].y << " " << orbit05pos[i].z << "\n";
  }
  
  // Interpolate to get positions at the actual times of the observations
  obspos={};
  for(obsct=0;obsct<obsnum;obsct++) {
    planetpos01LD(obsMJD[obsct], polyorder, orbit05MJD, orbit05pos, outpos);
    obspos.push_back(outpos);
  }
  
  // Find the corresponding topocentric celestial coords, and compare the
  // predicted positions to those that were observed.
  chisq = 0.0l;
  for(obsct=0;obsct<obsnum;obsct++) {
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    stateunitLD_to_celestial(outpos, outRA, outDec);
    cout  << fixed << setprecision(6) << "Input MJD " << obsMJD[obsct] << ": " << obsRA[obsct] << " "  << obsDec[obsct] << " "  << " Output: " << outRA << ": " << outDec <<  "\n";
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    chisq += LDSQUARE(dval/sigastrom[obsct]);
  }
  cout << "Chi-square value for input state vector is " << chisq << "\n";

  fitDec = fitRA = fitresid = {};
  chisq = tortoisechi03(polyorder,planetnum,planetmjd,planetmasses,planetpos,observerpos,obsMJD,obsRA,obsDec,sigastrom,startpos,startvel,planetfile_refpoint,fitRA,fitDec,fitresid);
  cout << "Chi-square value for input state vector is " << chisq << "\n";
  bestchi = chisq;
  bestRA = fitRA;
  bestDec = fitDec;
  bestresid = fitresid;
  outpos = startpos;
  outvel = startvel;
  // Load a scaled version of the starting state vectors
  svec[0] = startpos.x;
  svec[1] = startpos.y;
  svec[2] = startpos.z;
  svec[3] = startvel.x*pvtimescale;
  svec[4] = startvel.y*pvtimescale;
  svec[5] = startvel.z*pvtimescale;

  // Intialize the approximate Hessian matrix H to the identity matrix
  for(i=0;i<6;i++) {
    for(j=0;j<6;j++) {
      if(i==j) H[i][j]=1.0l;
      else H[i][j]=0.0l;
    }
  }
  
  itct=chifailct=goodstep=0;
  chichange = LARGERR;
  while(chisq>1.0 && itct<ITMAX && pos_step>MINPOS_STEP && chifailct<=CHIFAILMAX) {
    itct+=1;
    // Calculate the gradient of the chi-squared surface
    orbgrad01b(chisq, pos_step, pvtimescale, polyorder, planetnum, planetmjd, planetmasses, planetpos, observerpos, obsMJD, obsRA, obsDec, sigastrom, startpos, startvel, planetfile_refpoint, gradient);
    for(i=0;i<6;i++) gradient[i] *= -1.0l; // Now it really is the gradient.
    // Calculate unit vector corresponding to the gradient
    gradunit = gradient;
    ldval = nvecnorm(gradunit);
    // Take the dot product of the gradient with its own unit vector
    ldval = nvecdotprod(gradunit,gradient);
    // This should be the expected change in the chi-square value,
    // given a one-unit change in the position in the ideal direction.
    poschange = chisq*CHI_REDFAC/ldval;
    // This should be the distance we have to go in the ideal direction to
    // reduce the chi-square value by the factor CHI_REDFAC

    if(goodstep==1) {
      // BFGS step: add a correction to the Hessian matrix H
      // Calculate state vector difference
      statevecdiff={};
      for(i=0;i<6;i++) statevecdiff.push_back(svec[i] - oldsvec[i]);
      // Calculate gradient difference
      gradientdiff={};
      for(i=0;i<6;i++) gradientdiff.push_back(gradient[i] - oldgrad[i]);
      // Normalize both difference vectors
      nvecnorm(statevecdiff);
      nvecnorm(gradientdiff);
      cout << "State vector diff: ";
      for(i=0;i<6;i++) cout << statevecdiff[i] << " ";
      cout << "\n";
      cout << "Gradient diff: ";
      for(i=0;i<6;i++) cout << gradientdiff[i] << " ";
      cout << "\n";
      
      // FIRST TERM OF CORRECTION TO THE HESSIAN MATRIX
      status = vector_outerprod(statevecdiff,statevecdiff,OP1);
      if(status!=0) {
	cerr << "Failure in vector_outerprod(statevecdiff,statevecdiff,OP1);\n";
	return(status);
      }
      // Take the dot product of the state vector difference and the gradient difference
      dotprod_gradstate = nvecdotprod(statevecdiff,gradientdiff);
      // Divide the matrix OP1 by this factor
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) OP1[i][j]/=dotprod_gradstate;
      }
      cout << "OP1:\n";
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) {
	  cout << OP1[i][j] << " ";
	}
	cout << "\n";
      }
      
      // SECOND TERM OF CORRECTION TO THE HESSIAN MATRIX
      // Transform the gradient difference by the old versio of the Hessian
      status = matXvec(H, gradientdiff, Hgrad);
      if(status!=0) {
	cerr << "Failure in matXvec(H, gradientdiff, Hgrad);\n";
	return(status);
      }
      status = vector_outerprod(Hgrad,Hgrad,OP2);
      if(status!=0) {
	cerr << "Failure in vector_outerprod(Hgrad,Hgrad,OP2);\n";
	return(status);
      }
      // Take the dot-product of the gradient difference with the new,
      // scaled version of it that is vec1
      dotprod_gradHgrad = nvecdotprod(Hgrad,gradientdiff);
      
      // Divide the matrix OP2 by this factor
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) OP2[i][j]/=dotprod_gradHgrad;
      }
      cout << "OP2:\n";
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) {
	  cout << OP2[i][j] << " ";
	}
	cout << "\n";
      }
       
      // THIRD TERM OF CORRECTION TO THE HESSIAN MATRIX
      // construct first part of u vector
      u = statevecdiff;
      for(i=0;i<6;i++) u[i]/=dotprod_gradstate;
      // Apply correction to produce final u vector
      for(i=0;i<6;i++) u[i] -= Hgrad[i]/dotprod_gradHgrad;
      // Calculate outer product
      status = vector_outerprod(u,u,OP3);
      if(status!=0) {
	cerr << "Failure in vector_outerprod(u,u,OP3);\n";
	return(status);
      }
      // Scale the outer product
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) OP3[i][j]*=dotprod_gradHgrad;
      }
      cout << "OP3:\n";
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) {
	  cout << OP3[i][j] << " ";
	}
	cout << "\n";
      }
 
      // Add all the terms of the Hessian together
      for(i=0;i<6;i++) {
	for(j=0;j<6;j++) H[i][j] += OP1[i][j] - OP2[i][j] + OP3[i][j];
      }
    }
    cout << "Hessian:\n";
    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) {
	cout << H[i][j] << " ";
      }
      cout << "\n";
    }
       
    // BFGS step: multiply the gradient unit vector by the Hessian matrix to get the new direction
    matXvec(H, gradient, statevecdiff);
    thisway = statevecdiff;
    ldval = nvecnorm(thisway); // Renormalize
    for(i=0;i<6;i++) thisway[i] *= -1.0l; // Sign-flip: we want to go DOWN, not up.
    cout << "Search vector: " << fixed << setprecision(6) << thisway[0] << " "   << thisway[1] << " "   << thisway[2] << " "   << thisway[3] << " "   << thisway[4] << " "   << thisway[5] << "\n";

    // Perform one-dimensional optimization along the direction defined by
    // the Hessian-scaled gradient vector thisway
    newchi = orb1Dmin01b(chisq, poschange, pvtimescale, minchichange, thisway, polyorder, planetnum, planetmjd, planetmasses, planetpos, observerpos, obsMJD, obsRA, obsDec, sigastrom, startpos, startvel, planetfile_refpoint, fitRA, fitDec, fitresid, &beststep, tppos, tpvel);
    
    cout << "Iteration " << itct << ": chi-square value for input state vector is " << fixed << setprecision(3) << newchi << " step: " << fixed << setprecision(3) << pos_step  << " motion: " << fixed << setprecision(3) << beststep << ", chichange = " << chichange << "\n";

    if(newchi>0.0l && newchi<chisq) {
      // The 1-D search succeeded.
      chichange = chisq-newchi;
      if(chichange<minchichange) chifailct++;
      if(beststep<pos_step) {
	shortstepct++;
	if(shortstepct>SHORTSTEPMAX) {
	  shortstepct=0;
	  pos_step/=2.0l;
	}
      }
      chisq = newchi;
      startpos = tppos;
      startvel = tpvel;
      bestchi = newchi;
      bestRA = fitRA;
      bestDec = fitDec;
      bestresid = fitresid;
      outpos = startpos;
      outvel = startvel;
      goodstep=1;
      // BFGS stuff starts here
      oldsvec = svec;
      oldgrad = gradient;
      svec[0] = startpos.x;
      svec[1] = startpos.y;
      svec[2] = startpos.z;
      svec[3] = startvel.x*pvtimescale;
      svec[4] = startvel.y*pvtimescale;
      svec[5] = startvel.z*pvtimescale;
    } else {
      // The 1-D search failed. The gradient must have been bad.
      // Shorten the scale used for calculating the derivative
      pos_step/=2.0l;
      goodstep=0;
    }
  }
  long double astromrms=0.0l;
  for(obsct=0;obsct<obsnum;obsct++) astromrms += bestresid[obsct]*bestresid[obsct];
  astromrms = sqrt(astromrms/double(obsnum));
  
  cout << "Best chi-squared value was " << bestchi << ", astrometric RMS = " << astromrms << "\n";
  cout << fixed << setprecision(10) << "Best state vectors at MJD " << planetmjd[planetfile_refpoint] << " : " << fixed << setprecision(3) << outpos.x << " " << outpos.y << " " << outpos.z << " "  << fixed << setprecision(10) << outvel.x << " " << outvel.y << " " << outvel.z << "\n";

  ofstream outstream1 {outfile};
  for(obsct=0;obsct<obsnum;obsct++) {
    outstream1 << fixed << setprecision(8) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << bestRA[obsct] << " " << bestDec[obsct] << fixed << setprecision(8) << " " << (obsRA[obsct]-bestRA[obsct])*cos(obsDec[obsct]/DEGPRAD)*3600.0l << " " << (obsDec[obsct]-bestDec[obsct])*3600.0l << " " << bestresid[obsct] << "\n";
    }
  outstream1.close();
  
  return(0);
}

