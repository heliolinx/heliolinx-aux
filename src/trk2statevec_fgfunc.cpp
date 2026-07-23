#define MAXTANVELCUT 50.0 // Maximum value that can be placed on the minimum tangential velocity in km/sec
                          // for a valid tracklet in heliolinc. Note that setting it anywhere near
                          // this would reject almost all tracklets.

int trk2statevec_fgfunc(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar, double max_v_inf, int NotKepler)
{
  allstatevecs={};
  long imnum = image_log.size();
  long imct=0;
  long pairnum = tracklets.size();
  long pairct=0;
  int badpoint=0;
  int status1=0;
  int status2=0;
  int num_dist_solutions=0;
  int solnct=0;
  double mjdavg=0l;
  vector <double> heliodistvec;
  double delta1 = 0.0l;
  double RA,Dec;
  long i1,i2;
  i1=i2=0;
  point6dx2 statevec1 = point6dx2(0l,0l,0l,0l,0l,0l,0,0);
  point6ix2 stateveci = point6ix2(0,0,0,0,0,0,0,0);
  point3d observerpos1 = point3d(0l,0l,0l);
  point3d observerpos2 = point3d(0l,0l,0l);
  point3d targpos1 = point3d(0l,0l,0l);
  point3d targpos2 = point3d(0l,0l,0l);
  point3d targvel1 = point3d(0l,0l,0l);
  point3d targvel2 = point3d(0l,0l,0l);
  point3d unitbary = point3d(0l,0l,0l);
  vector <point3d> targposvec1;
  vector <point3d> targposvec2;
  int glob_warning=0;
  vector <double> deltavec1;
  vector <double> deltavec2;
  double absvelocity=0l;
  double impactpar=0l;
  double timediff=0l;
  double E = 0.0l;
  double v_inf = 0.0l;
  double obstanvel = MAXTANVELCUT;
  double meanobsdist,trackletarc,trackletangvel;
  meanobsdist = trackletarc = trackletangvel = 0.0;
 
  // Calculate approximate heliocentric distances from the
  // input quadratic approximation.
  // Note: heliodist is in km, heliovel is in km/day, and helioacc is in km/day^2
  heliodistvec={};
  if(NotKepler) {
    // Use the old-style, three-term Taylor Series approximation of heliocentric distance
    // as a function of time.
    for(imct=0;imct<imnum;imct++) {
      delta1 = image_log[imct].MJD - mjdref; // Units are days
      heliodistvec.push_back(heliodist + heliovel*delta1 + 0.5*helioacc*delta1*delta1);
      if(heliodistvec[imct]<=0.0l) {
	badpoint=1;
	return(1);
      }
    }
  } else {
    // Use Matt Holman's idea of getting heliocentric distance r(t) from an
    // actual Keplerian orbit. It turns out that the same parameters (r, r-dot, and r-double-dot)
    // that specify the old three-term Taylor Series also uniquely determine the Keplerian
    // r(t) -- that is, they specify enough of the orbit that only one Keplerian solution
    // for r(t) exists, even though other aspects of the orbit (i.e., orientation) are not
    // specified.
    // SOLVE FOR THE SQUARE OF THE TANGENTIAL VELOCITY
    double localg = GMSUN_KM3_SEC2/DSQUARE(heliodist); // Units are km/sec^2
    double physacc = helioacc/DSQUARE(SOLARDAY); // Units are km/sec^2
    double vesc = 2.0*GMSUN_KM3_SEC2/heliodist; // this is the square of the escape velocity in km/sec
    double tanvel = heliodist*(physacc+localg); // this is the square of the tangential velocity in km/sec
    // CHECK FOR UNPHYSICAL AND UNBOUND CASES
    if(tanvel<0.0l) {
      cerr << fixed << setprecision(6) << "ERROR: hypothesis point " << heliodist/AU_KM << ", " << heliovel/AU_KM << ", " << -helioacc/DSQUARE(SOLARDAY)/localg << " is not possible for any trajectory\n";
      return(1);
    }
    if(vesc < tanvel + LDSQUARE(heliovel/SOLARDAY)) {
      cerr << fixed << setprecision(6) << "ERROR: hypothesis point " << heliodist/AU_KM << ", " << heliovel/AU_KM << ", " << -helioacc/DSQUARE(SOLARDAY)/localg << " is not possible for a bound orbit\n";
      cerr << fixed << setprecision(6) << "vesc = " << sqrt(vesc) << ", tanvel = " << sqrt(tanvel) << ", heliovel = " << heliovel/SOLARDAY << ", totvel = " << sqrt(tanvel + LDSQUARE(heliovel/SOLARDAY)) << "\n";
      return(1);
    }
    // If we get here, a sensible bound solution exists. Solve for the true tangential velocity
    tanvel = sqrt(tanvel);
    // Construct state vectors producing the required orbit (in the x-y plane, for simplicity).
    point3d startpos = point3d(heliodist,0l,0l);
    point3d startvel = point3d(heliovel/SOLARDAY,tanvel,0l);
    point3d endpos = point3d(0l,0l,0l);
    point3d endvel = point3d(0l,0l,0l);
    for(imct=0;imct<imnum;imct++) {
      // Integrate the orbit to find the heliocentric distance as a function of time.
      status1 = Kepler_fg_func_int(GMSUN_KM3_SEC2, mjdref, startpos, startvel, image_log[imct].MJD, endpos, endvel);
      if(status1!=0) {
	cerr << "ERROR: Keplerian integration failed for r(t) hypothesis point " << heliodist/AU_KM << ", " << heliovel/AU_KM << ", " << -helioacc/DSQUARE(SOLARDAY)/localg << ", at MJD = " << image_log[imct].MJD << "\n";
	return(status1);
      }
      heliodistvec.push_back(vecabs3d(endpos));
    }
  }
  if(badpoint==0 && long(heliodistvec.size())!=imnum) {
    cerr << "ERROR: number of heliocentric distance values does\n";
    cerr << "not match the number of input images!\n";
    return(2);
  }
  for(pairct=0; pairct<pairnum; pairct++) {
    badpoint=0;
    // Obtain indices to the image_log and heliocentric distance vectors.
    i1=tracklets[pairct].Img1;
    i2=tracklets[pairct].Img2;
    // Project the first point
    RA = tracklets[pairct].RA1;
    Dec = tracklets[pairct].Dec1;
    celestial_to_stateunit(RA,Dec,unitbary);
    observerpos1 = point3d(image_log[i1].X,image_log[i1].Y,image_log[i1].Z);
    targposvec1={};
    deltavec1={};
    status1 = helioproj02(unitbary,observerpos1, heliodistvec[i1], deltavec1, targposvec1);
    RA = tracklets[pairct].RA2;
    Dec = tracklets[pairct].Dec2;
    celestial_to_stateunit(RA,Dec,unitbary);
    observerpos2 = point3d(image_log[i2].X,image_log[i2].Y,image_log[i2].Z);
    targposvec2={};
    deltavec2={};
    status2 = helioproj02(unitbary, observerpos2, heliodistvec[i2], deltavec2, targposvec2);
    // Calculate the mean distance from the observer over the two points in the tracklet.
    // If the projection equations had two solutions (meaning the hypothesis distance was
    // interior to the Earth), use the larger solutions.
    meanobsdist = (deltavec1[0] + deltavec2[0])/2.0;
    if(status1==2 && status2==2 &&  deltavec1[1]>deltavec1[0] && deltavec2[1]>deltavec2[0]) meanobsdist = (deltavec1[1] + deltavec2[1])/2.0;
    if(minimpactpar > 0.0 && minimpactpar < MAXTANVELCUT && meanobsdist < mingeoobs*AU_KM && status1 > 0 && status2 > 0 && badpoint==0) {
      // New check added May 13, 2026: Calculate observer-centric tangential velocity,
      // and reject the point if that is too low. The value of minimpactpar is here
      // interpreted as a tangential velocity in km/sec, since it is nonzero but too small to
      // make sense as an impact parameter in km.
      trackletarc = distradec01(tracklets[pairct].RA1, tracklets[pairct].Dec1, tracklets[pairct].RA2, tracklets[pairct].Dec2)/DEGPRAD; // Tracklet arc in radians
      timediff = (image_log[i2].MJD - image_log[i1].MJD)*SOLARDAY; // Time difference in seconds
      trackletangvel = trackletarc/timediff; // radians per second
      obstanvel = trackletangvel*meanobsdist; // km/sec
      if(obstanvel<minimpactpar) {
	// Tangential velocity relative to the observer is too low
	badpoint=1;
      }
    }
    if(status1 > 0 && status2 > 0 && badpoint==0) {
      // Calculate time difference between the observations
      timediff = (image_log[i2].MJD - image_log[i1].MJD)*SOLARDAY;
      // Did helioproj find two solutions in both cases, or only one?
      num_dist_solutions = status1;
      if(num_dist_solutions > status2) num_dist_solutions = status2;
      // Loop over solutions (num_dist_solutions can only be 1 or 2).
      for(solnct=0; solnct<num_dist_solutions; solnct++) {
	// Calculate the object's v_inf relative to the sun.
	targpos1 = targposvec1[solnct];
	targpos2 = targposvec2[solnct];
	  
	targvel1.x = (targpos2.x - targpos1.x)/timediff;
	targvel1.y = (targpos2.y - targpos1.y)/timediff;
	targvel1.z = (targpos2.z - targpos1.z)/timediff;

	targpos1.x = 0.5L*targpos2.x + 0.5L*targpos1.x;
	targpos1.y = 0.5L*targpos2.y + 0.5L*targpos1.y;
	targpos1.z = 0.5L*targpos2.z + 0.5L*targpos1.z;

	E = 0.5l*dotprod3d(targvel1,targvel1) - GMSUN_KM3_SEC2/vecabs3d(targpos1);
	if(E>0.0l) v_inf = sqrt(2.0l*E);
	else if(!isnormal(E)) v_inf=0.0l;
	else v_inf = -sqrt(-2.0l*E); // This is a bit weird, but we allow the user to
	                             // set a negative v_inf, if desired, to rule out
	                             // objects that are barely bound to the sun.
	if(v_inf>max_v_inf) continue; // Skip further calculation if v_inf is too high.

	// Begin new stuff added to eliminate 'globs'
	// These are spurious linkages of unreasonably large numbers (typically tens of thousands)
	// of detections that arise when the hypothetical heliocentric distance at a time when
	// many observations are acquired is extremely close to, but slightly greater than,
	// the heliocentric distance of the observer. Then detections over a large area of sky
	// wind up with projected 3-D positions in an extremely small volume -- and furthermore,
	// they all have similar velocities because the very small geocentric distance causes
	// the inferred velocities to be dominated by the observer's motion and the heliocentric
	// hypothesis, with only a negligible contribution from the on-sky angular velocity.
	glob_warning=0;
	if(deltavec1[solnct]<mingeoobs*AU_KM && deltavec2[solnct]<mingeoobs*AU_KM) {
	  // New-start
	  // Load target positions
	  targpos1 = targposvec1[solnct];
	  targpos2 = targposvec2[solnct];
	  // Calculate positions relative to observer
	  targpos1.x -= observerpos1.x;
	  targpos1.y -= observerpos1.y;
	  targpos1.z -= observerpos1.z;
	    
	  targpos2.x -= observerpos2.x;
	  targpos2.y -= observerpos2.y;
	  targpos2.z -= observerpos2.z;
	    
	  // Calculate velocity relative to observer
	  targvel1.x = (targpos2.x - targpos1.x)/timediff;
	  targvel1.y = (targpos2.y - targpos1.y)/timediff;
	  targvel1.z = (targpos2.z - targpos1.z)/timediff;
   
	  // Calculate impact parameter (past or future).
	  absvelocity = vecabs3d(targvel1);
	  impactpar = dotprod3d(targpos1,targvel1)/absvelocity;
	  // Effectively, we've projected targpos1 onto the velocity
	  // vector, and impactpar temporarily holds the magnitude of this projection.
	  // Subtract off the projection of the distance onto the velocity unit vector
	  targpos1.x -= impactpar*targvel1.x/absvelocity;
	  targpos1.y -= impactpar*targvel1.y/absvelocity;
	  targpos1.z -= impactpar*targvel1.z/absvelocity;
	  // Now targpos1 is the impact parameter vector at projected closest approach.
	  impactpar  = vecabs3d(targpos1); // Now impactpar is really the impact parameter
	  if(impactpar<=minimpactpar) {
	    // The hypothesis implies the object already passed within minimpactpar km of the Earth
	    // in the likely case that minimpactpar has been set to imply an actual impact,
	    // it's not our problem anymore.
	    glob_warning=1;
	  }
	}
	if(!glob_warning) {
	  targpos1 = targposvec1[solnct];
	  targpos2 = targposvec2[solnct];
	  
	  targvel1.x = (targpos2.x - targpos1.x)/timediff;
	  targvel1.y = (targpos2.y - targpos1.y)/timediff;
	  targvel1.z = (targpos2.z - targpos1.z)/timediff;

	  targpos1.x = 0.5L*targpos2.x + 0.5L*targpos1.x;
	  targpos1.y = 0.5L*targpos2.y + 0.5L*targpos1.y;
	  targpos1.z = 0.5L*targpos2.z + 0.5L*targpos1.z;
      
	  // Integrate orbit to the reference time.
	  mjdavg = 0.5l*image_log[i1].MJD + 0.5l*image_log[i2].MJD;
	  status1 = Kepler_fg_func_int(GMSUN_KM3_SEC2,mjdavg,targpos1,targvel1,mjdref,targpos2,targvel2);
	  if(status1 == 0 && badpoint==0) {
	    statevec1 = point6dx2(targpos2.x,targpos2.y,targpos2.z,chartimescale*targvel2.x,chartimescale*targvel2.y,chartimescale*targvel2.z,pairct,0);
	    // Note that the multiplication by chartimescale converts velocities in km/sec
	    // to units of km, for apples-to-apples comparison with the positions.
	    stateveci = conv_6d_to_6i(statevec1,INTEGERIZING_SCALEFAC);
	    allstatevecs.push_back(stateveci);
	  } else {
	    // Kepler integration encountered unphysical situation.
	    continue;
	  }
	}
      }
    } else {
      badpoint=1;
      // Heliocentric projection found no physical solution.
      continue;
    }
  }
  return(0);
}
