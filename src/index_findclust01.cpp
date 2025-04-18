// April 02, 2025: index_findclust01.cpp
//
// 
#include "solarsyst_dyn_geo01.h"
#include "cmath"


static void show_usage()
{
  cerr << "Usage: index_findclust01 -indexlist index_list -clust2det input cluster-to-detection file -out output file\n";
}

int main(int argc, char *argv[])
{
  string clust2detfile,indexlist,outfile;
  vector <longpair> inclust2det;
  long one_index=0;
  vector <long> indices;
  vector <long> clustind;
  long thiscluster=0;
  ifstream instream1;
  ofstream outstream1;
  int status=0;
  long clustct,i,j,matchnum,bestmatchnum,bestmatch,bestsize,fullmatch;
  clustct=i=j=matchnum=bestmatchnum=bestmatch=fullmatch=bestsize=0;
  int verbose=0;
  
  if(argc<7) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-indexlist" || string(argv[i]) == "-index" || string(argv[i]) == "-indices" || string(argv[i]) == "-list" || string(argv[i]) == "--indexlist" || string(argv[i]) == "--index_list") {
      if(i+1 < argc) {
	//There is still something to read;
	indexlist=argv[++i];
	i++;
      }
      else {
	cerr << "Index list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clust2det" || string(argv[i]) == "-c2d" || string(argv[i]) == "-inc2d" || string(argv[i]) == "-input_c2d" || string(argv[i]) == "--input_clust2det" ) {
      if(i+1 < argc) {
	//There is still something to read;
	clust2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outclust" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outclust" || string(argv[i]) == "--output_cluster_file") {
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

  // Catch required parameters if missing
  if(indexlist.size()<=0) {
    cout << "\nERROR: input index list file is required\n";
    show_usage();
    return(1);
  } else if(clust2detfile.size()<=0) {
    cout << "\nERROR: input cluster-to-detection file is required\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cout << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  cout << "input index list " << indexlist << "\n";
  cout << "input cluster-to-detection file " << clust2detfile << "\n";
  cout << "output file " << outfile << "\n";

  // Read input index list file
  instream1.open(indexlist);
  indices={};
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    one_index=-1;
    instream1 >> one_index;
    if(one_index>=0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) indices.push_back(one_index);
    else if(instream1.bad()) {
      cerr << "WARNING: problem reading file " << indexlist << " after line " << indices.size() << "\n";
    }
  }
  instream1.close();
  if(indices.size()<=0) {
    cerr << "ERROR: could not read any lines from index list " << indexlist << "\n";
    return(1);
  }
  cout << "Read " << indices.size() << " lines from index list file " << clust2detfile << "\n";
  sort(indices.begin(),indices.end());
  
  inclust2det={};
  // Read cluster-to-detection file
  status=read_longpair_file(clust2detfile, inclust2det, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read cluster-to-detection file " << clust2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
    return(1);
  }
  cout << "Read " << inclust2det.size() << " data lines from cluster-to-detection file " << clust2detfile << "\n";

  outstream1.open(outfile);
  clustind = {};
  thiscluster = inclust2det[0].i1;
  clustind.push_back(inclust2det[0].i2);
  for(clustct=0;clustct<long(inclust2det.size());clustct++) {
    if(inclust2det[clustct].i1==thiscluster) {
      // Add a new entry to the cluster index vector
      clustind.push_back(inclust2det[clustct].i2);
    } else {
      // Analyze this cluster
      sort(clustind.begin(),clustind.end());
      i=j=matchnum=0;
      while(i<long(indices.size()) && j<long(clustind.size())) {
	if(indices[i]==clustind[j]) {
	  // Found a match!
	  matchnum++;
	  i++;
	  j++;
	} else if(indices[i]>clustind[j]) j++;
	else if(indices[i]<clustind[j]) i++;
	else {
	  cerr << "Logically excluded comparison case: i,j,indices[i],clustind[j] = " << i << "," << j << "," << indices[i] << "," << clustind[j] << "\n";
	  return(2);
	}
      }
      if(matchnum>=long(indices.size())) {
	cout << "Match to all " << indices.size() << " lines found for cluster " << thiscluster << " which has a total of " << clustind.size() << " lines\n";
	outstream1 << "Match to all " << indices.size() << " lines found for cluster " << thiscluster << " which has a total of " << clustind.size() << " lines\n";
	for(j=0;j<long(clustind.size());j++) {
	  outstream1 << thiscluster << "," << clustind[j] << "\n";
	}
	fullmatch=1;
      }
      if(matchnum>bestmatchnum) {
	bestmatchnum=matchnum;
	bestmatch=thiscluster;
	bestsize=clustind.size();
      }
      // Set up for the next one
      clustind = {};
      thiscluster = inclust2det[clustct].i1;
      clustind.push_back(inclust2det[clustct].i2);
    }
  }
  if(fullmatch==1) {
    cout << "A fully inclusive match was found: cluster " << bestmatch << " with " << bestsize << " lines\n";
    outstream1 << "A fully inclusive match was found: cluster " << bestmatch << " with " << bestsize << " lines\n";
  } else {
    cout << "No complete match was found. Nearest approach was cluster " << bestmatch << ", which matched " << bestmatchnum << " out of " << indices.size() << " detections\n";
    cout << "This best-matching cluster had a total of " << bestsize << " lines\n";
    outstream1 << "No complete match was found. Nearest approach was cluster " << bestmatch << ", which matched " << bestmatchnum << " out of " << indices.size() << " detections\n";
    outstream1 << "This best-matching cluster had a total of " << bestsize << " lines\n";
  }
  outstream1.close();
  
  return(0);
}
