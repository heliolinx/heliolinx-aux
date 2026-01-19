// February 15, 2025: hldet2mpc.cpp
// Converts a file in hldet format to MPC80 format.
// Assumes there is no header.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: hldet2ITFjson.cpp hldet_file jsonfile argperi eccen epoch_JD incl long_asc_node pericenter_dist percenter_JD\n";
}

int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  string hldetfile,jsonfile;
  long i=0;
  long j=0;
  long status=0;
  int year,month;
  double day=0.0;
  long iday=0;
  string datestring;
  string idstring;
  string oldidstring;
  string obscodestring;
  string catstring;
  vector <string> desigvec;
  char c;
  int verbose=1;
  int write_orbit = 1;
  
  if(argc!=10 && argc!=3) {
    show_usage();
    return(1);
  } if(argc!=10) write_orbit=0;
  
  hldetfile=argv[1];
  jsonfile=argv[2];

  cout.precision(17);  
  cout << "input observation file " << hldetfile << "\n";
  cout << "output file " << jsonfile << "\n";
  
  // Read input observation file.
  detvec={};
  status=read_hldet_file(hldetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read observation file " << hldetfile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  long obsnum = detvec.size();
  cout << "Read " << obsnum << " data lines from observation file " << hldetfile << "\n";
  
  ofstream outstream1 {jsonfile};

  outstream1 << "{\n";
  outstream1 << "  \"header\": {\n";
  outstream1 << "    \"name\": \"A. N. Heinze\",\n";
  outstream1 << "    \"email\": \"ariheinze@hotmail.com\",\n";
  outstream1 << "    \"comment\": \"A linkage from the ITF.\"\n";
  outstream1 << "  },\n";
  outstream1 << "  \"links\": {\n";
  outstream1 << "    \"link_0\": {\n";
  outstream1 << "      \"trksubs\": [\n";
  
  desigvec={};
  oldidstring="";
  for(i=0; i<obsnum; i++) {
    idstring = detvec[i].idstring;
    obscodestring = detvec[i].obscode;
    mjd2mpcdate(detvec[i].MJD,year,month,day);
    iday = floor(day);
    datestring = to_string(year);
    if(month<10) datestring += "0";
    datestring += to_string(month);
    if(iday<10) datestring += "0";
    datestring += to_string(iday);
    // cout << "Date check : " << year << " " << to_string(year) << " " << month << " " << to_string(month) << " " << iday << " " << to_string(iday) << " " << datestring << "\n";
    catstring = datestring + "," + idstring + "," + obscodestring;
    if(idstring!=oldidstring) {
      desigvec.push_back(catstring);
      oldidstring = idstring;
    }
  }
  sort(desigvec.begin(),desigvec.end());
  for(i=0;i<long(desigvec.size());i++) {
    j=0;
    c='0';
    datestring="";
    while(j<long(desigvec[i].size()) && c!=',') {
      c=desigvec[i][j];
      if(c!=',') datestring += c;
      j++;
    }
    if(j<long(desigvec[i].size())) c=desigvec[i][j];
    idstring="";
    while(j<long(desigvec[i].size()) && c!=',') {
      c=desigvec[i][j];
      if(c!=',') idstring += c;
      j++;
    }
    if(j<long(desigvec[i].size())) c=desigvec[i][j];
    obscodestring="";
    while(j<long(desigvec[i].size()) && c!=',') {
      c=desigvec[i][j];
      if(c!=',') obscodestring += c;
      j++;
    }
    outstream1 << "        [\n";
    outstream1 << "          \"" << idstring << "\",\n";
    outstream1 << "          \"" << datestring << "\",\n";
    outstream1 << "          \"" << obscodestring << "\"\n";
    outstream1 << "        ]";
    if(i<long(desigvec.size()-1)) {
      // This is not the last tracklet. Add a comma.
      outstream1 << ",\n";
    } else outstream1 << "\n";
  }
  if(write_orbit>0) {
    outstream1 << "      ],\n";
    outstream1 << "      \"orbit\": {\n";
    outstream1 << "        \"arg_pericenter\": " << argv[3] << ",\n";
    outstream1 << "        \"eccentricity\": " << argv[4] << ",\n";
    outstream1 << "        \"epoch\": " << argv[5] << ",\n";
    outstream1 << "        \"inclination\": " << argv[6] << ",\n";
    outstream1 << "        \"lon_asc_node\": " << argv[7] << ",\n";
    outstream1 << "        \"pericenter_distance\": " << argv[8] << ",\n";
    outstream1 << "        \"pericenter_time\": " << argv[9] << "\n";
    outstream1 << "      }\n";
  } else outstream1 << "      ]\n";
  outstream1 << "    }\n";
  outstream1 << "  }\n";
  outstream1 << "}\n";
  
  outstream1.close();
  
  return(0);
}

