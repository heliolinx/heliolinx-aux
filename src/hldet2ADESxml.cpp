// March 12, 2025: hldet2ADESxml: Convert a file of heliolinc hldet
// format into an ADES file.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: hldet2ADESxml -infile hldetfile -head headerfile -timeshift timeshift -permID permID -provID provID -artSat artSat -trkSub trkSub -out outfile\n";
}

int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  string hldetfile,outfile,headerfile;
  string lnfromfile;
  long i=0;
  long status=0;
  long obsnum=0;
  int verbose=0;
  int year,month;
  double day;
  string signstring;
  long double timeshift = 0.0l;
  ofstream outstream1;
  ifstream instream1;
  string permID,provID,artSat;
  string trkSub="a000000";

  if(argc<9) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-hlfile" || string(argv[i]) == "-hldet" || string(argv[i]) == "-hldetfile") {
      if(i+1 < argc) {
	//There is still something to read;
	hldetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input hldet file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-head" || string(argv[i]) == "-header" || string(argv[i]) == "-headfile" || string(argv[i]) == "-headerfile" || string(argv[i]) == "--head" || string(argv[i]) == "--header" || string(argv[i]) == "--headfile" || string(argv[i]) == "--headerfile") {
      if(i+1 < argc) {
	//There is still something to read;
	headerfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input header file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-permID" || string(argv[i]) == "-permid") {
      if(i+1 < argc) {
	//There is still something to read;
	permID=argv[++i];
	i++;
      }
      else {
	cerr << "Input permID keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-provID" || string(argv[i]) == "-provid") {
      if(i+1 < argc) {
	//There is still something to read;
	provID=argv[++i];
	i++;
      }
      else {
	cerr << "Input provID keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-artSat" || string(argv[i]) == "-artsat") {
      if(i+1 < argc) {
	//There is still something to read;
	artSat=argv[++i];
	i++;
      }
      else {
	cerr << "Input artSat keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trkSub" || string(argv[i]) == "-trksub") {
      if(i+1 < argc) {
	//There is still something to read;
	trkSub=argv[++i];
	i++;
      }
      else {
	cerr << "Input trkSub keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timeshift" || string(argv[i]) == "-time" || string(argv[i]) == "--timeshift" || string(argv[i]) == "-timediff") {
      if(i+1 < argc) {
	//There is still something to read;
	timeshift = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--output_ADES_file") {
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
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
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

  // Catch required parameters if missing
  if(hldetfile.size()<=0) {
    cerr << "\nERROR: input hldet file is required\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cerr << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  cout << "input hldet file " << hldetfile << "\n";
  cout << "time shift " << timeshift << " seconds\n";
  cout << "output file " << outfile << "\n";

  // Read input observation file.
  detvec={};
  status=read_hldet_file(hldetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << hldetfile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  obsnum = detvec.size();
  cout << "Read " << obsnum << " data lines from observation file " << hldetfile << "\n";
  // Time-sort detvec, just to be sure
  sort(detvec.begin(), detvec.end(), early_hldet());
 
  outstream1.open(outfile);

  // Write introductory ADES stuff
  outstream1 << "<?xml version=\'1.0\' encoding=\'UTF-8\'?>\n";
  outstream1 << "<ades version=\"2022\">\n";

  outstream1 << "  <obsBlock>\n";
  outstream1 << "    <obsContext>\n";
  instream1.open(headerfile);
  if(!instream1) {
    cerr << "can't open input file " << headerfile << "\n";
    return(1);
  }
  while (!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(lnfromfile.size()>0) outstream1 << lnfromfile << "\n";
  }
  instream1.close();
  outstream1 << "      <comment>\n";
  outstream1  << fixed << setprecision(3) << "        <line> These observations formatted by Ari Heinze\'s program hldet2ADESxml </line>\n";
  outstream1 << "      </comment>\n";
  outstream1 << "    </obsContext>\n";
  outstream1 << "    <obsData>\n";
  // Write this data to the output file
  for(i=0; i<long(detvec.size()); i++) {
    outstream1 << "      <optical>\n";
    if(permID.size()>0) outstream1 << "        <permID>" << permID << "</permID>\n";
    else if(provID.size()>0) outstream1 << "        <provID>" << provID << "</provID>\n";
    else if(artSat.size()>0) outstream1 << "        <artSat>" << artSat << "</artSat>\n";
    outstream1 << "        <trkSub>" << trkSub << "</trkSub>\n";
    outstream1 << "        <obsSubID>" << detvec[i].idstring << "</obsSubID>\n";
    outstream1 << "        <mode>CCD</mode>\n";
    outstream1 << "        <stn>" << detvec[i].obscode << "</stn>\n";
    // Adjust the input MJD by timeshift, which could be positive (e.g., offset to mid-exposure)
    // or negative (e.g., to convert from TAI to UTC).
    mjd2mpcdate(detvec[i].MJD +(timeshift/SOLARDAY), year, month, day);
    long iday = floor(day);
    double fracday = day - double(iday);
    double hour = floor(24.0*fracday);
    double minute = floor(1440.0*fracday - 60.0*hour);
    double second = fracday*SOLARDAY - 3600.0*hour - 60.0*minute;
    outstream1 << "        <obsTime>" << year << "-";
    if(month<10) outstream1 << "0";
    outstream1 << month << "-";
    if(iday<10) outstream1 << "0";
    outstream1 << iday << "T";
    if(hour<10.0) outstream1 << "0";
    outstream1 << fixed << setprecision(0) << hour << ":";
    if(minute<10.0) outstream1 << "0";
    outstream1 << fixed << setprecision(0) << minute << ":";
    if(second<9.9995) outstream1 << "0";
    outstream1 << fixed << setprecision(3) << second << "Z</obsTime>\n";
    outstream1  << fixed << setprecision(6) << "        <ra>" << detvec[i].RA << "</ra>\n";
    outstream1  << fixed << setprecision(6) << "        <dec>" << detvec[i].Dec << "</dec>\n";
    if(detvec[i].sig_across>0 && detvec[i].sig_along>0) {
      outstream1  << fixed << setprecision(3) << "        <rmsRA>" << detvec[i].sig_across << "</rmsRA>\n";
      outstream1  << fixed << setprecision(3) << "        <rmsDec>" << detvec[i].sig_along << "</rmsDec>\n";
    }
    outstream1  << "        <astCat>Gaia2</astCat>\n";
    outstream1  << fixed << setprecision(3) << "        <mag>" << detvec[i].mag << "</mag>\n";
    outstream1  << fixed << setprecision(3) << "        <rmsMag>" << detvec[i].sigmag << "</rmsMag>\n";
    outstream1  << "        <band>" << detvec[i].band << "</band>\n";
    outstream1 << "      </optical>\n";
  }
  outstream1 << "    </obsData>\n";

  outstream1 << "  </obsBlock>\n";

  outstream1 << "</ades>\n";
  outstream1.close();
  return(0);
}
