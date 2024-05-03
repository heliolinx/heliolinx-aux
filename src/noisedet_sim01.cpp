// May 02, 2024: noisedet_sim01.cpp:
// Given an input image file of the form written by make_tracklets_new
// and read by heliolinc_kd, simulate a bunch of false detections. They
// are intended either to be added on top of existing simulations of real
// data, or else used alone as a test to see how frequently heliolinc will
// generate linkages of purely spurious points.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: noisedet_sim01 -imgs imfile -ranseed random_number_seed -imrad image radius(deg) -ndet bad detections per image -mag dummy magnitude -idroot id_rootstring -band band -out outputfile\n";
  cerr << "\nor, at minimum:\n\n";
  cerr << "noisedet_sim01 -imgs imfile -ranseed random_number_seed\n";
}
    
int main(int argc, char *argv[])
{
  vector <hlimage> image_log;
  string imfile;
  string outfile = "noisedets01a.csv";
  string band = "V";
  string idroot = "noisedet";
  double mag = 99.9;
  double imagerad = 2.0;
  long ndet = 100;
  int default_outfile, default_band, default_imagerad;
  int default_mag, default_ndet, default_idroot;
  default_outfile = default_band = default_imagerad = default_mag = default_ndet = default_idroot = 1;
  ofstream outstream1;
  string seedstring;
  long i;
  long imct;
  int status=0;
  double x, y, rad, theta, RA, Dec;
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	//There is still something to read;
	imfile=argv[++i];
	i++;
      }
      else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rs" || string(argv[i]) == "-ranseed" || string(argv[i]) == "-rseed" || string(argv[i]) == "-seed" || string(argv[i]) == "--ranseed" || string(argv[i]) == "--randomseed" || string(argv[i]) == "--randomnumberseed" || string(argv[i]) == "--random_number_seed") {
      if(i+1 < argc) {
	// There is still something to read;
	seedstring=argv[++i];
	i++;
      } else {
	cerr << "Random number seed keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-imagerad" || string(argv[i]) == "-imagerad" ) {
      if(i+1 < argc) {
	//There is still something to read;
        imagerad=stod(argv[++i]);
	default_imagerad = 0;
	i++;
	if(!isnormal(imagerad) || imagerad<=0.0) {
	  cerr << "Error: invalid image radius (" << imagerad << " deg) supplied.\n";
	  cerr << "Image radius must be strictly positive!\n";
	  return(2);
	}
      }
      else {
	cerr << "Image radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-ndet" || string(argv[i]) == "-detnum" ) {
      if(i+1 < argc) {
	//There is still something to read;
        ndet=stol(argv[++i]);
	default_ndet = 0;
	i++;
	if(ndet<=0) {
	  cerr << "Error: invalid number of detections per image (" << ndet << ") supplied.\n";
	  cerr << "Number of detections per image must be strictly positive!\n";
	  return(2);
	}
      }
      else {
	cerr << "Number of detection per image keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mag" || string(argv[i]) == "-magnitude" ) {
      if(i+1 < argc) {
	//There is still something to read;
        mag=stod(argv[++i]);
	default_mag = 0;
	i++;
      }
      else {
	cerr << "Magnitude keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-idroot" || string(argv[i]) == "-IDroot" || string(argv[i]) == "--root" || string(argv[i]) == "--idroot") {
      if(i+1 < argc) {
	//There is still something to read;
	idroot=argv[++i];
	default_idroot = 0;
	i++;
      }
      else {
	cerr << "ID root keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-band" || string(argv[i]) == "-bnd" || string(argv[i]) == "--band" || string(argv[i]) == "--bnd") {
      if(i+1 < argc) {
	//There is still something to read;
	band=argv[++i];
	default_band = 0;
	i++;
      }
      else {
	cerr << "Band keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--out" || string(argv[i]) == "--outfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	default_outfile = 0;
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
  cout << "Seedstring is " << seedstring << "\n";
  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine
  cout << "Random: " << unitvar(generator) << " " << unitvar(generator) << " " << unitvar(generator) << " " << unitvar(generator) << " " << unitvar(generator) << "\n";

  
  if(argc<5)
    {
      cerr << "Too few arguments even for minimalist invocation:\n";
      show_usage();
      return(1);
    }
  
  cout.precision(17);  
  cout << "input image file " << imfile << "\n";
  
  if(default_imagerad==1) cout << "Defaulting to image radius = " << imagerad << "deg\n";
  else cout << "input image radius " << imagerad << "deg\n";
  if(default_ndet==1) cout << "Defaulting to " << ndet << " detections per image\n";
  else cout << "Using user-selected value of " << ndet << " detections per image\n";
  if(default_ndet==1) cout << "Defaulting to dummy magnitude = " << mag << "\n";
  else cout << "Dummy magnitude will be " << mag << "\n";
  if(default_idroot==1) cout << "Defaulting to ID root name " << idroot << "\n";
  else cout << "ID root name will be " << idroot << "\n";
  if(default_band==1) cout << "Defaulting to photometric band " << band << "\n";
  else cout << "Photomeric band will be " << band << "\n";
  if(default_outfile==1) cout << "Defaulting to output file name " << outfile << "\n";
  else cout << "Output file will be named " << outfile << "\n";
  
  image_log={};
  status=read_image_file2(imfile, image_log);
  if(status!=0) {
    cerr << "ERROR: could not successfully read image file " << imfile << "\n";
    cerr << "read_image_file2 returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << image_log.size() << " data lines from image file " << imfile << "\n";
  

  outstream1.open(outfile);
  outstream1 << "#IDstring,MJD,RA,Dec,mag,band,ObsCode\n";
  for(imct=0; imct<long(image_log.size()); imct++) {
    i=0;
    while(i<ndet) {
      x = 2.0l*imagerad*unitvar(generator) - imagerad;
      y = 2.0l*imagerad*unitvar(generator) - imagerad;
      rad = sqrt(x*x + y*y);
      if(rad<=imagerad) {
	if(x==0l && y>=0l) theta = 0l;
	else if(x==0l && y<0l) theta=180.0l;
	else if(x>0) theta = 90.0l - DEGPRAD*atan(y/x);
	else if(x<0) theta = 270.0l - DEGPRAD*atan(y/x);
	else {
	  cerr << "Illogical trigonmetry case with x = " << x << " and y = " << y << "\n";
	  return(1);
	}
	i++;
	arc2cel01(image_log[imct].RA,image_log[imct].Dec,rad,theta,RA,Dec);
	outstream1 << idroot << imct << ",";
	outstream1 << fixed << setprecision(8) << image_log[imct].MJD << ",";
	outstream1 << fixed << setprecision(6) << RA << "," << Dec << ",";
	outstream1 << fixed << setprecision(3) << mag << "," << band << "," << image_log[imct].obscode << "\n";
      }
    }
  }
  outstream1.close();
  return(0);
}
