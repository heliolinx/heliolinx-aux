// June 20, 2023
// Read 2 files with 5 columns: a string ID, followed by four floating point numbers.
// Load the first file into a KD tree, and do a range-query on this KD tree around
// every point in the second file. Create an output file giving all the pair-matches
// that could be identified.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: 4D_KDtree_rangequery01 -pool file of potential matches -lonely file of things to be matched -range range for KD query -outfile output file\n";
}
  
int main(int argc, char *argv[])
{
  string pool_file, lonely_file, outfile;
  point4d_index onepoint = point4d_index(0,0,0,0,0);
  vector <point4d_index> poolvec;
  vector <point4d_index> lonelyvec;
  string desig;
  vector <string> pool_desig;
  vector <string> lonely_desig;
  double range = 0.0l;
  ifstream instream1;
  ofstream outstream1;
  double E, Lx, Ly, Lz;
  long i=0;
  long j=0;
  long kdroot=0;
  long splitpoint=0;
  KD_point4d_index kdpoint = KD_point4d_index(onepoint,-1,-1,1,0);
  vector <KD_point4d_index> kdvec;
  vector <long> indexvec;
  
  E = Lx = Ly = Lz = 0.0l;
  
  if(argc<9) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pool") {
      if(i+1 < argc) {
	//There is still something to read;
	pool_file=argv[++i];
	i++;
      }
      else {
	cerr << "Input match pool file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-lonely") {
      if(i+1 < argc) {
	//There is still something to read;
	lonely_file=argv[++i];
	i++;
      }
      else {
	cerr << "Input lonely file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-range") {
      if(i+1 < argc) {
	//There is still something to read;
	range = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Query range keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile") {
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

  // Read file of potential matches ('pool' file)
  instream1.open(pool_file);
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    instream1 >> desig >> E >> Lx >> Ly >> Lz;
    i = poolvec.size();
    if(i!=long(pool_desig.size())) {
      cerr << "Error: vector size mismatch " << poolvec.size() << " " << pool_desig.size() << "\n";
      return(1);
    }
    onepoint = point4d_index(E,Lx,Ly,Lz,i);
    poolvec.push_back(onepoint);
    pool_desig.push_back(desig);
  }
  instream1.close();

  // Read file of things to be matched ('lonely' file)
  instream1.open(lonely_file);
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    instream1 >> desig >> E >> Lx >> Ly >> Lz;
    i = lonelyvec.size();
    if(i!=long(lonely_desig.size())) {
      cerr << "Error: vector size mismatch " << lonelyvec.size() << " " << lonely_desig.size() << "\n";
      return(1);
    }
    onepoint = point4d_index(E,Lx,Ly,Lz,i);
    lonelyvec.push_back(onepoint);
    lonely_desig.push_back(desig);
  }
  instream1.close();

  // Form KDtree
  kdvec={};
  kdroot = splitpoint = 0;
  splitpoint=medind_4d_index(poolvec,1);
  kdpoint = KD_point4d_index(poolvec[splitpoint],-1,-1,1,0);
  kdvec.push_back(kdpoint);
  kdtree_4d_index(poolvec,1,splitpoint,kdroot,kdvec);
  cout << "Created a KD tree with " << kdvec.size() << " branches\n";

  // Perform range queries about each point in the 'lonely' vectors
  outstream1.open(outfile);
  for(i=0;i<long(lonelyvec.size());i++) {
    cout << "Querying around point " << i << ": " << lonely_desig[i] << "\n";
    indexvec={};
    kdrange_4d_index(kdvec, lonelyvec[i], range, indexvec);
    if(indexvec.size()>0) {
      // We found a match!
      cout << "MATCHES FOUND: ";
      outstream1 << lonely_desig[i] << " ";
      for(j=0;j<long(indexvec.size());j++) {
	cout << pool_desig[kdvec[indexvec[j]].point.index] << " ";
	outstream1 << pool_desig[kdvec[indexvec[j]].point.index] << " ";
      }
      cout << "\n";
      outstream1 << "\n";
    } else {
      cout << "no match found\n";
    }
  }

  outstream1.close();
  return(0);
}
