// October 22, 2025: label_hldet01.cpp
// Reads two files of hldet format, an 'unlabled' file and a
// 'labeling' file. Constructs a 4-dimensional k-d tree of the
// 'unlabeled' file, which is expected to be larger than the
// 'labeling' file. Performs a nearest-neighbor query of the
// 'unlabeled' k-d tree around each point in the 'labeling'
// file. If there is a match within a given tolerance, the
// nearest neighbor is than labeled with string ID from the
// 'labeling' file, which supercedes the string ID originally
// read from the 'unlabeled' file, if any. Finally, the newly
// labeled 'unlabeled' file is written out. This output file
// is guaranteed to have the same number of lines as the
// input 'unlabeled' file, and the same string ID's except
// where matches were found.

#define DAY_TO_DEG_CONV 24.0 // Treats one day as equivalent to 24 degrees:
                             // in other words, one second of time is
                             // equivalent to one arcsecond on the sky.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: label_hldet01 -unlabeled unlabeled_file -label labeling_file -matchrad match_in_arcsec -timescale day_to_deg_conversion -outfile outfile\n";
}

int main(int argc, char *argv[])
{
  ofstream outstream1;
  vector <hldet> unlabeled_dets = {};
  vector <hldet> labeling_dets = {};
  string unlabeled_file,  labeling_file, outfile;
  double timescale = DAY_TO_DEG_CONV;
  double matchrad = 1.0;
  long ulnum,labnum,ulct,labct,i,j,bandlen,status;
  ulnum = labnum = ulct = labct = i = j = bandlen = status = 0;
  int verbose=0;
  long kdroot=0;
  long splitpoint=0;
  long nearest=0;
  long index=0;
  point4d_index onepoint = point4d_index(0,0,0,0,0);
  point4d_index querypoint = point4d_index(0,0,0,0,0);
  vector <point4d_index> poolvec;
  KD_point4d_index kdpoint = KD_point4d_index(onepoint,-1,-1,1,0);
  vector <KD_point4d_index> kdvec;
  vector <long> indexvec;
  vector <double> mjdvec;
  double mjdref = 0.0;
  double dist = 0.0;

  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-unlabeled" || string(argv[i]) == "-unlab" || string(argv[i]) == "-ul" || string(argv[i]) == "-nolab" || string(argv[i]) == "--unlabeled" || string(argv[i]) == "--unlab") {
      if(i+1 < argc) {
	//There is still something to read;
	unlabeled_file=argv[++i];
	i++;
      }
      else {
	cerr << "Unlabeled file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-labeled" || string(argv[i]) == "-labeling" || string(argv[i]) == "-lab" || string(argv[i]) == "-labeling_file" || string(argv[i]) == "-label" || string(argv[i]) == "--labeled") {
      if(i+1 < argc) {
	//There is still something to read;
	labeling_file=argv[++i];
	i++;
      }
      else {
	cerr << "Labeling file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-matchrad") {
      if(i+1 < argc) {
	//There is still something to read;
	matchrad=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Matching radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timescale") {
      if(i+1 < argc) {
	//There is still something to read;
	timescale=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Timescale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Labeling file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }
  cout.precision(17);  
  
  // Read input unlabeled file
  unlabeled_dets={};
  status=read_pairdet_file(unlabeled_file, unlabeled_dets, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read input unlabled file " << unlabeled_file << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  ulnum = unlabeled_dets.size();
  cout << "Read " << ulnum << " data lines from observation file " << unlabeled_file << "\n";
  
  // Read input labeling file
  labeling_dets={};
  status=read_pairdet_file(labeling_file, labeling_dets, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read input labeling file " << labeling_file << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  labnum = labeling_dets.size();
  cout << "Read " << labnum << " data lines from observation file " << labeling_file << "\n";

  // Find the median MJD for unlabeled detections
  mjdvec = {};
  for(ulct=0;ulct<ulnum;ulct++) mjdvec.push_back(unlabeled_dets[ulct].MJD);
  mjdref = dmedian(mjdvec);
  cout << "Median MJD of unlabled detections is " << mjdref << "\n";
  
  // Load pool of unlabeled sources
  poolvec = {};
  for(ulct=0;ulct<ulnum;ulct++) {
    // Note that timescale is supposed to indicate the number of
    // degrees on the sky that is the equivalent of one calendar day.
    // Then, since the trigonometric projection of sky-angles converts
    // degrees into radians, it is also necessary to divide the
    // time-quantity by DEGPRAD.
    onepoint = point4d_index((unlabeled_dets[ulct].MJD-mjdref)*timescale/DEGPRAD,cos(unlabeled_dets[ulct].RA/DEGPRAD)*cos(unlabeled_dets[ulct].Dec/DEGPRAD),sin(unlabeled_dets[ulct].RA/DEGPRAD)*cos(unlabeled_dets[ulct].Dec/DEGPRAD), sin(unlabeled_dets[ulct].Dec/DEGPRAD), ulct);
    poolvec.push_back(onepoint);
  }
  cout << "Finished loading pool of unlabeled sources\n";
  // Form KDtree
  kdvec={};
  kdroot = splitpoint = 0;
  splitpoint=medind_4d_index(poolvec,1);
  kdpoint = KD_point4d_index(poolvec[splitpoint],-1,-1,1,0);
  kdvec.push_back(kdpoint);
  kdtree_4d_index(poolvec,1,splitpoint,kdroot,kdvec);
  cout << "Finished constructing k-d tree sources\n";

  // Loop over pool of labeled sources, performing a nearest-neighbor
  // query of the unlabled k-d tree for each of them.
  for(labct=0;labct<labnum;labct++) {
    querypoint = point4d_index((labeling_dets[labct].MJD-mjdref)*timescale/DEGPRAD,cos(labeling_dets[labct].RA/DEGPRAD)*cos(labeling_dets[labct].Dec/DEGPRAD),sin(labeling_dets[labct].RA/DEGPRAD)*cos(labeling_dets[labct].Dec/DEGPRAD), sin(labeling_dets[labct].Dec/DEGPRAD), labct);
    nearest = kdnearest_4d_index(kdvec, querypoint);
    dist = ASECPRAD*sqrt(point4d_index_dist2(querypoint, kdvec[nearest].point));
    if(labct%1000==0) {
      cout << "Running query for labeling point " << labct << "\n";
      cout << querypoint.t << " " << querypoint.x << " " << querypoint.y << " " << querypoint.z << "\n";
      cout << "best match is " << kdvec[nearest].point.t << " " << kdvec[nearest].point.x << " " << kdvec[nearest].point.y << " " << kdvec[nearest].point.z << "\n";
      cout << "Distance is " << dist << "\n";
    }
    if(dist<=matchrad) {
      // Label the match with an index that maps all the way back to unlabeled_dets
      index = kdvec[nearest].point.index;
      // Change the label for unlabeled_dets
      for(i=0;i<SHORTSTRINGLEN;i++) unlabeled_dets[index].idstring[i] = labeling_dets[labct].idstring[i];
    }
  }
  
  outstream1.open(outfile);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(ulct=0;ulct<ulnum;ulct++) {
    outstream1 << fixed << setprecision(7) << unlabeled_dets[ulct].MJD << "," << unlabeled_dets[ulct].RA << "," << unlabeled_dets[ulct].Dec << ",";
    outstream1 << fixed << setprecision(4) << unlabeled_dets[ulct].mag << ",";
    outstream1 << fixed << setprecision(2) << unlabeled_dets[ulct].trail_len << "," << unlabeled_dets[ulct].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << unlabeled_dets[ulct].sigmag << ",";
    outstream1 << fixed << setprecision(3) << unlabeled_dets[ulct].sig_across << "," << unlabeled_dets[ulct].sig_along << ",";
    outstream1 << unlabeled_dets[ulct].image << "," << unlabeled_dets[ulct].idstring << "," << unlabeled_dets[ulct].band << ",";
    outstream1 << unlabeled_dets[ulct].obscode << "," << unlabeled_dets[ulct].known_obj << ","; 
    outstream1 << unlabeled_dets[ulct].det_qual << "," << unlabeled_dets[ulct].index << "\n"; 
  }
  outstream1.close();

  return(0);
}
