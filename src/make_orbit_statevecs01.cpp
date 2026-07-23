// June 22, 2026: make_orbit_statevecs01.cpp
// Given an input x, y, z position and an allowed maximum eccentricity
// and inclination, plus sampling intervals in velocity (in km/sec),
// eccentricty (unitless) and angle (in degrees), create a set of state
// vectors that spans the space of possible Keplerian orbits passing
// through the specified x, y, z position and having eccentricity
// and inclination within the specified ranges.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

void vector_crossprod(const vector <double> &vec1, const vector <double> &vec2, vector <double> &outvec);

static void show_usage()
{
  cerr << "Usage: make_orbit_statevecs01 -x xpos -y ypos -z zpos -maxecc max_eccentricity -maxinc max_inclination -velsamp velocity_sampling -eccsamp eccentricity_sampling -angsamp angle_sampling -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  string outfile;
  ofstream outstream1;
  double x,y,z,maxecc,maxinc,velsamp,eccsamp,angsamp;
  vector <double> posvec;
  vector <double> statevec;
  double E,L,r0,vmin,vmax,vel;
  double MGsun = GMSUN_KM3_SEC2;
  make_dvec(6,statevec);
  long i,vnum,vct,eccen_num,eccen_ct;
  int status=0;
  double eccen;
  double vpar,vperp,norm,theta,sine_rvtheta;
  double a, e, incl;
  vector <double> perpvec1;
  vector <double> perpvec2;
  vector <double> unitvec;
  vector <double> tempvec;
  vector <double> velvec;
  make_dvec(3,posvec);
  make_dvec(3,perpvec1);
  make_dvec(3,perpvec2);
  make_dvec(3,unitvec);
  make_dvec(3,velvec);
  velsamp = 1.0;
  angsamp = 5.0;
  eccsamp = 0.01;
  maxinc = 50.0;
  maxecc = 0.50;
  x=y=z=0.0;

  if(argc<17) {
    show_usage();
    return(1);
  }
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile") {
      if(i+1 < argc) {
	// There is still something to read;
	outfile=argv[++i];
	i++;
      } else {
	cerr << "Output file file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-x"  || string(argv[i]) == "-xpos") {
      if(i+1 < argc) {
	//There is still something to read;
	x=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "x position keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-y"  || string(argv[i]) == "-ypos") {
      if(i+1 < argc) {
	//There is still something to read;
	y=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "y position keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-z"  || string(argv[i]) == "-zpos") {
      if(i+1 < argc) {
	//There is still something to read;
	z=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "z position keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxecc"  || string(argv[i]) == "-max_eccen") {
      if(i+1 < argc) {
	//There is still something to read;
	maxecc=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "max eccentricity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxinc"  || string(argv[i]) == "-max_inc") {
      if(i+1 < argc) {
	//There is still something to read;
	maxinc=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "max inclination keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velsamp"  || string(argv[i]) == "-vel_samp") {
      if(i+1 < argc) {
	//There is still something to read;
	velsamp=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "velocity sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-eccsamp"  || string(argv[i]) == "-ecc_samp") {
      if(i+1 < argc) {
	//There is still something to read;
	eccsamp=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "eccentricity sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }else if(string(argv[i]) == "-angsamp"  || string(argv[i]) == "-ang_samp") {
      if(i+1 < argc) {
	//There is still something to read;
	angsamp=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "angle sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
  r0 = sqrt(x*x + y*y + z*z);
  statevec[0] = posvec[0] = x;
  statevec[1] = posvec[1] = y;
  statevec[2] = posvec[2] = z;
  outstream1.open(outfile);
  eccen_num = maxecc/eccsamp + 1;
  for(eccen_ct=0; eccen_ct<eccen_num; eccen_ct++) {
    eccen = maxecc*(0.5+double(eccen_ct))/double(eccen_num);
    vmin = sqrt(MGsun * (1.0 - eccen)/r0);
    vmax = sqrt(MGsun * (1.0 + eccen)/r0);
    vnum = (vmax-vmin)/velsamp + 1;
    for(vct=0;vct<vnum;vct++) {
      vel = vmin + (vmax-vmin)*(0.5+double(vct))/double(vnum);
      E = 0.5*vel*vel - MGsun/r0;
      L = sqrt((eccen*eccen-1.0)*MGsun*MGsun/2.0/E);
      sine_rvtheta = L/r0/vel;
      if(sine_rvtheta>=1) {
	cerr << "WARNING: sine rvtheta = " << sine_rvtheta << "\n";
	sine_rvtheta = 1.0;
      }
      vpar = vel*sqrt(1.0-sine_rvtheta*sine_rvtheta);
      vperp = vel*sine_rvtheta;
      // Find two unit vectors that are perpendicular to the position vector
      norm=0.0;      
      for(i=0;i<3;i++) {
	unitvec[0] = unitvec[1] = unitvec[2] = 0.0;
	unitvec[i] = 1.0;
	vector_crossprod(statevec, unitvec, tempvec);
	if(nvecabs(tempvec) > norm) {
	  norm = nvecabs(tempvec);
	  perpvec1 = tempvec;
	}
      }
      // Now perpvec1 is the longest vector that is perpendicular both to
      // one of the original basis vectors and the position vector. The
      // remaining vector we need is perpendicular to both perpvec1 and the
      // position vector.
      vector_crossprod(statevec, perpvec1, perpvec2);
      nvecnorm(perpvec1);
      nvecnorm(perpvec2);
      // Loop over all the possible ways to rotate the velocity vector
      // around the axis defined by the position vector
      theta = 0.0;
      while(theta<360.0) {
	velvec[0] = vpar*posvec[0]/r0 + vperp*(perpvec1[0]*cos(theta/DEGPRAD) + perpvec2[0]*sin(theta/DEGPRAD));
	velvec[1] = vpar*posvec[1]/r0 + vperp*(perpvec1[1]*cos(theta/DEGPRAD) + perpvec2[1]*sin(theta/DEGPRAD));
	velvec[2] = vpar*posvec[2]/r0 + vperp*(perpvec1[2]*cos(theta/DEGPRAD) + perpvec2[2]*sin(theta/DEGPRAD));
	// Calculate the inclination.
	statevec[3] = velvec[0];
	statevec[4] = velvec[1];
	statevec[5] = velvec[2];
	a = e = incl = 0;
	status = statevec2kep_easy(MGsun, statevec, a, e, incl);
	cout << "a = " << a/AU_KM << ", e = " << e << " should be " << eccen << ", incl = " << incl << "\n";
	// If inclination is within range, write to the output file
	if(incl<=maxinc) {
	  outstream1 << fixed << setprecision(4) << statevec[0] << " " << statevec[1] << " " << statevec[2] << " ";
	  outstream1 << fixed << setprecision(8) << statevec[3] << " " << statevec[4] << " " << statevec[5] << " ";
	  outstream1 << fixed << setprecision(4) << a << " " << e << " " << incl << "\n";
	}
	theta+=angsamp;
      }
    }
  }
  outstream1.close();
  return(0);
}

void vector_crossprod(const vector <double> &vec1, const vector <double> &vec2, vector <double> &outvec)
{
  outvec={};
  make_dvec(3,outvec);
  outvec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  outvec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  outvec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}
