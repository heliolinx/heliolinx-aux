// October 10, 2025: Horizons_resample.cpp
// Resample a JPL Horizons output file to align with
// the Everhart spacing.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define HNUM 4
#define TIMEMATCHTOL 1.0e-10l

static void show_usage()
{
  cerr << "Usage: Horizons_resample -infile infile -timestep timestep_in_days -hnum hnum -polyorder polyorder -outfile outfile_forward outfile_reverse\n";
}

int main(int argc, char *argv[])
{
  int status=0;
  ifstream instream1;
  ofstream outstream1;
  string stest;
  string planetfile;
  string outfile1;
  string outfile2;
  long double timestep = 5.0;
  int polyorder=5;
  long hnum = HNUM;
  vector <long double> hspace;
  vector <long double> htimes;
  vector <long double> outstate;
  vector <vector <long double>> statevecs;
  vector <long double> planetmjd;
  long itimestep;
  long iMJD;
  long mjdnum,mjdct;
  long i=0;
  long double mjd0=0;
  long double mjdnow=0;

  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-in") {
      if(i+1 < argc) {
	//There is still something to read;
	planetfile=argv[++i];
	i++;
      }
      else {
	cerr << "infile keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timestep" || string(argv[i]) == "-ts" || string(argv[i]) == "-step" || string(argv[i]) == "--timestep" || string(argv[i]) == "--time_step_days" || string(argv[i]) == "--ModifiedJulianDayref" || string(argv[i]) == "--timestep") {
      if(i+1 < argc) {
	//There is still something to read;
	timestep=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Timestep keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-polyorder") {
      if(i+1 < argc) {
	//There is still something to read;
	polyorder=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "polyorder keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-hnum") {
      if(i+1 < argc) {
	//There is still something to read;
	hnum=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "hnum keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outorb" || string(argv[i]) == "--outorbits") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile1=argv[++i];
	cout << "Read outfile1 as " << outfile1 << "\n";
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	outfile2=argv[++i];
	cout << "Read outfile2 as " << outfile2 << "\n";
	i++;
      }
      else {
	cerr << "Output file keyword supplied with only one out of two required arguments\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  cout.precision(17);  
  cout << "output file " << outfile1 << "\n";

  // Load hspace vector based on selected value of hnum
  if(hnum<3 || hnum>8) {
    cerr << "Warning: hnum set to out-of-range value (" << hnum << "), will be reset to default of " << HNUM << "\n";
    hnum=HNUM;
  }
  hspace={};
  if(hnum==3) {
    hspace.push_back(h3_1);
    hspace.push_back(h3_2);
    hspace.push_back(h3_3);
  } else if(hnum==4) {
    hspace.push_back(h4_1);
    hspace.push_back(h4_2);
    hspace.push_back(h4_3);
    hspace.push_back(h4_4);
  } else if(hnum==5) {
    hspace.push_back(h5_1);
    hspace.push_back(h5_2);
    hspace.push_back(h5_3);
    hspace.push_back(h5_4);
    hspace.push_back(h5_5);
  } else if(hnum==6) {
    hspace.push_back(h6_1);
    hspace.push_back(h6_2);
    hspace.push_back(h6_3);
    hspace.push_back(h6_4);
    hspace.push_back(h6_5);
    hspace.push_back(h6_6);
  } else if(hnum==7) {
    hspace.push_back(h7_1);
    hspace.push_back(h7_2);
    hspace.push_back(h7_3);
    hspace.push_back(h7_4);
    hspace.push_back(h7_5);
    hspace.push_back(h7_6);
    hspace.push_back(h7_7);
  } else if(hnum==8) {
    hspace.push_back(h8_1);
    hspace.push_back(h8_2);
    hspace.push_back(h8_3);
    hspace.push_back(h8_4);
    hspace.push_back(h8_5);
    hspace.push_back(h8_6);
    hspace.push_back(h8_7);
    hspace.push_back(h8_8);
  }
  make_LDvec(hnum+1,htimes);
  for(i=1;i<=hnum;i++) htimes[i] = timestep*hspace[i-1]; // Units are days

  outstream1.open(outfile1);
  outstream1.precision(17);
  status = read_horizons_statevecLD(planetfile,planetmjd,statevecs);
  if(status!=0) {
    cerr << "ERROR: read_horizons_statevecLD() returned error status " << status << "\n";
    return(status);
  }
  mjdnum = planetmjd.size();
  if(mjdnum<=0) {
    cerr << "ERROR: no lines read from the planet file\n";
    return(1);
  }
  // Check to make sure the timestep is an integer
  itimestep = round(timestep);
  if(fabs(timestep - static_cast<long double>(itimestep))>TIMEMATCHTOL) {
    cerr << "ERROR: timestep " << timestep << " does not appear to be an integer number of days\n";
    return(1);
  }
  // Find the first mjd in the planet file that is an integer multiple of timestep
  mjdct=0;
  iMJD = round(planetmjd[mjdct]);
  while(iMJD % itimestep != 0  &&  mjdct < mjdnum-1) {
    mjdct++;
    iMJD = round(planetmjd[mjdct]);
  }
  cout << "Output file will start at point " << mjdct << " of the input ephemeris file,\n";
  cout << "which corresponds to MJD " << planetmjd[mjdct] << "\n";

  while(mjdct+itimestep<mjdnum) {
    mjd0 = planetmjd[mjdct];
    for(i=1;i<=hnum;i++) {
      mjdnow = mjd0+htimes[i];
      status = planetpos02LD(mjdnow, polyorder, planetmjd, statevecs, outstate);
      if(status!=0) {
	cerr << "ERROR: planetpos02LD() returned error status " << status << "\n";
	return(status);
      }
      outstream1 << fixed << setprecision(11) << mjdnow << " ";
      outstream1 << scientific << setprecision(16) << outstate[0] << " " << outstate[1] << " " << outstate[2] << " " << outstate[3] << " " << outstate[4] << " " << outstate[5] << "\n";
      cout << i << " " << mjdnow << " " << outstate[0] << " " << outstate[1] << " " << outstate[2] << " " << outstate[3] << " " << outstate[4] << " " << outstate[5] << "\n";
    }
    mjdct+=itimestep;
  }
  if(mjdct<mjdnum) {
    /// Write out one more line, an integer step
     outstream1 << fixed << setprecision(11) << planetmjd[mjdct] << " ";
     outstream1 << scientific << setprecision(16) << statevecs[mjdct][0] << " "<< statevecs[mjdct][1] << " " << statevecs[mjdct][2] << " " << statevecs[mjdct][3] << " " << statevecs[mjdct][4] << " " << statevecs[mjdct][5] << "\n";
    cout << i << " " << planetmjd[mjdct] << " " << statevecs[mjdct][0] << " " << statevecs[mjdct][1] << " " << statevecs[mjdct][2] << " " << statevecs[mjdct][3] << " " << statevecs[mjdct][4] << " " << statevecs[mjdct][5] << "\n";
  }
  outstream1.close();

  // Now for the reverse output file
  outstream1.open(outfile2);
  outstream1.precision(17);
  // Find the last mjd in the planet file that is an integer multiple of timestep
  mjdct=planetmjd.size()-1;
  iMJD = round(planetmjd[mjdct]);
  while(iMJD % itimestep != 0  &&  mjdct >0) {
    mjdct--;
    iMJD = round(planetmjd[mjdct]);
  }
  while(mjdct-itimestep>=0) {
    mjd0 = planetmjd[mjdct];
    for(i=1;i<=hnum;i++) {
      mjdnow = mjd0-htimes[i];
      planetpos02LD(mjdnow, polyorder, planetmjd, statevecs, outstate);
      // Sign-flip time and velocity
      outstream1 << fixed << setprecision(11) << -mjdnow << " ";
      outstream1 << scientific << setprecision(16) << outstate[0] << " " << outstate[1] << " " << outstate[2] << " " << -outstate[3] << " " << -outstate[4] << " " << -outstate[5] << "\n";
      cout << i << " " << -mjdnow << " " << outstate[0] << " " << outstate[1] << " " << outstate[2] << " " << -outstate[3] << " " << -outstate[4] << " " << -outstate[5] << "\n";
    }
    mjdct-=itimestep;
  }
  if(mjdct>=0) {
    // Write one more line, an integer step.
    outstream1 << fixed << setprecision(11) << -planetmjd[mjdct] << " ";
    outstream1 << scientific << setprecision(16) << statevecs[mjdct][0] << " " << statevecs[mjdct][1] << " " << statevecs[mjdct][2] << " " << -statevecs[mjdct][3] << " " << -statevecs[mjdct][4] << " " << -statevecs[mjdct][5] << "\n";
  }
  outstream1.close();
  return(0);
}
