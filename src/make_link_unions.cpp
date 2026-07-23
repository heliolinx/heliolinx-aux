// July 17, 2026: make_link_unions.cpp
// Reads two files in hldet format, a 'submitted' file
// with linkages flagged by the stringID (column 11), which
// is supposed to correspond to a trkSub actually submitted
// to the MPC, and a 'new' file with linkages flagged by the
// index (column 16).
// For each linkage in the 'new' file, identify any overlapping
// linkages in the 'submitted' file. Create a new hldet vector
// containing the union of the new linkage with any overlapping
// linkages from the old file. De-duplicate this new linkage,
// retaining RA and Dec uncertainties where available. Output
// the de-duplicated union into a new file. 

#define DAY_TO_DEG_CONV 24.0 // Treats one day as equivalent to 24 degrees:
                             // in other words, one second of time is
                             // equivalent to one arcsecond on the sky.
#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define COLS_TO_READ 4

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: make_link_unions -submitted submitted_file -new newfile -matchrad match_in_arcsec -timescale day_to_deg_conversion -minobsadd minimum_number_of_added_observations -mintimeadd minimum_time_added_to_linkage(days) -outfile outfile\n";
}

int main(int argc, char *argv[])
{
  ofstream outstream1;
  ofstream outstream2;
  ifstream instream1;
  ifstream instream2;
  string submitted_file,newfile,outfile,stest,idstring,newstring;
  string lnfromfile;
  vector <hldet> detvec_sub;
  vector <hldet> detvec_new;
  vector <hldet> newlinkage;
  vector <hldet> templinkage;
  vector <double> MJDvec;
  vector <double> RAvec;
  vector <double> Decvec;
  vector <double> linkMJD;
  long detnum_sub,detnum_new;
  detnum_sub = detnum_new = 0;
  vector <long> sub_indexvec;
  vector <long> ID_indexvec;
  vector <string> IDvec;
  vector <vector <string>> file_IDvecs;
  double timescale = DAY_TO_DEG_CONV;
  double matchrad = 1.0;
  double mintimeadd = 0.0;
  long minobsadd = 0;
  double sub_timespan,max_timespan,new_timespan,merged_timespan;
  sub_timespan = max_timespan = new_timespan = merged_timespan = 0.0;
  long i,j,bandlen,status,sub_obsnum,max_obsnum;
  i = j = bandlen = status = sub_obsnum = max_obsnum = 0;
  int enough_added = 1;
  long kdroot=0;
  long splitpoint=0;
  point4d_index onepoint = point4d_index(0,0,0,0,0);
  point4d_index querypoint = point4d_index(0,0,0,0,0);
  vector <point4d_index> poolvec;
  KD_point4d_index kdpoint = KD_point4d_index(onepoint,-1,-1,1,0);
  vector <KD_point4d_index> kdvec;
  vector <long> indexvec;
  double median_MJD;
  double mjdref = 0.0;
  vector <string> linklines;
  vector <long> linkdets;
  string thisID;
  string dupID;
  vector <string> svec;
  vector <string> svec2;
  vector <vector <string>> dupIDvec;
  long detct,detct_new,thisindex;
  detct = detct_new = thisindex = 0;
  int verbose=0;
  long orig_linkagesize,merged_linkagesize,final_linkagesize;
  orig_linkagesize = merged_linkagesize = final_linkagesize = 0;
  long goodpt,istimedup,outdetnum,timedupnum,linkagenum,outlinkagenum;
  goodpt=istimedup=outdetnum=timedupnum=linkagenum=outlinkagenum=0;
  long overlapnum,fully_redundant_ct,non_augment_ct;
  overlapnum = fully_redundant_ct = non_augment_ct = 0;
    
  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-submitted" || string(argv[i]) == "-sub" || string(argv[i]) == "-subfile" || string(argv[i]) == "-sub_file" || string(argv[i]) == "--submitted_file" || string(argv[i]) == "--submitted") {
      if(i+1 < argc) {
	//There is still something to read;
	submitted_file=argv[++i];
	i++;
      }
      else {
	cerr << "submitted file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-new" || string(argv[i]) == "-newfile" || string(argv[i]) == "-new_file" || string(argv[i]) == "--new_file" || string(argv[i]) == "--newfile" || string(argv[i]) == "--new") {
      if(i+1 < argc) {
	//There is still something to read;
	newfile=argv[++i];
	i++;
      }
      else {
	cerr << "new file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-matchrad" || string(argv[i]) == "-match_rad" || string(argv[i]) == "-mrad" || string(argv[i]) == "-matching_radius" || string(argv[i]) == "--matching_radius" || string(argv[i]) == "-searchrad") {
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
    } else if(string(argv[i]) == "-timescale" || string(argv[i]) == "-time_scale" || string(argv[i]) == "-tscale" || string(argv[i]) == "-day2deg" || string(argv[i]) == "-day2degscale" || string(argv[i]) == "-day2deg_scale"|| string(argv[i]) == "--timescale" || string(argv[i]) == "--time_scale" || string(argv[i]) == "--day2deg" || string(argv[i]) == "--day2degscale" || string(argv[i]) == "--day2deg_scale") {
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
    } else if(string(argv[i]) == "-minobsadd" || string(argv[i]) == "-minobs" || string(argv[i]) == "-minptadd") {
      if(i+1 < argc) {
	//There is still something to read;
	minobsadd=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum number of added observations keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mintimeadd" || string(argv[i]) == "-minarcadd" || string(argv[i]) == "-min_timespan_increase") {
      if(i+1 < argc) {
	//There is still something to read;
	mintimeadd=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum timespan addition keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--out_file") {
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
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }
  cout.precision(17);  

  cout << "Input submitted file is called " << submitted_file << "\n";
  cout << "Input new linkage file is called " << newfile << "\n";
  cout << "Matching radius will be " << matchrad << " arcseconds\n";
  cout << "Timescale for converting 1 day of time to equivalent degrees will be " << timescale << "\n";
  cout << "Minimum number of added observations for the merged linkage is " << minobsadd << "\n";
  cout << "Minimum increase in timespan for the merged linakes is " << mintimeadd << "\n";
  cout << "Output file will be called " << outfile << "\n";

  //Read the input submitted linkage file
  detvec_sub={};
  status=read_pairdet_file(submitted_file, detvec_sub, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read input submitted linkage file " << submitted_file << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  detnum_sub = detvec_sub.size();
  cout << "Read " << detnum_sub << " data lines from submitted linkage file " << submitted_file << "\n";
  // Set all of the index values with negative numbers to distinguish them
  // from the new linkage file.
  for(detct=0;detct<detnum_sub;detct++) detvec_sub[detct].index = -1 - abs(detvec_sub[detct].index);
  
  //Read the input new linkage file
  detvec_new={};
  status=read_pairdet_file(newfile, detvec_new, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read input new linkage file " << newfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  detnum_new = detvec_new.size();
  cout << "Read " << detnum_new << " data lines from new linkage file " << newfile << "\n";

  // Load MJDvec, RAvec, Decvec for the submitted linkage file
  MJDvec=RAvec=Decvec={};
  for(i=0;i<detnum_sub;i++) {
    MJDvec.push_back(detvec_sub[i].MJD);
    RAvec.push_back(detvec_sub[i].RA);
    Decvec.push_back(detvec_sub[i].Dec);
  }
  // Find the median MJD
  median_MJD = dmedian(MJDvec);
  mjdref = median_MJD;
  
  // Create a k-d tree
  poolvec = {};
  for(detct=0;detct<detnum_sub;detct++) {
    onepoint = point4d_index((MJDvec[detct]-mjdref)*timescale/DEGPRAD,cos(RAvec[detct]/DEGPRAD)*cos(Decvec[detct]/DEGPRAD),sin(RAvec[detct]/DEGPRAD)*cos(Decvec[detct]/DEGPRAD), sin(Decvec[detct]/DEGPRAD), detct);
    poolvec.push_back(onepoint);
  }
  cout << "Finished loading pool of " << detnum_sub << " detections for submitted linkage file\n";
  // Form KDtree
  kdvec = {};
  kdroot = splitpoint = 0;
  splitpoint=medind_4d_index(poolvec,1);
  kdpoint = KD_point4d_index(poolvec[splitpoint],-1,-1,1,0);
  kdvec.push_back(kdpoint);
  kdtree_4d_index(poolvec,1,splitpoint,kdroot,kdvec);
  cout << "Finished constructing k-d tree for submitted linkage file, with " << kdvec.size() << " items\n";
  
  // Loop over the new linkage file, identifying any duplicates.
  outstream1.open(outfile);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  detct = detct_new = thisindex = overlapnum = 0;
  thisindex = detvec_new[detct_new].index;
  newlinkage={};
  sub_indexvec = ID_indexvec = {};
  IDvec = {};
  while(detct_new <= detnum_new) {
    if(detct_new < detnum_new && detvec_new[detct_new].index == thisindex) {
      // This is another point in the same linkage. Add it to the vector.
      newlinkage.push_back(detvec_new[detct_new]);
      // See if there is a duplicate of this point in the submissions file.
      querypoint = point4d_index((detvec_new[detct_new].MJD-mjdref)*timescale/DEGPRAD,cos(detvec_new[detct_new].RA/DEGPRAD)*cos(detvec_new[detct_new].Dec/DEGPRAD),sin(detvec_new[detct_new].RA/DEGPRAD)*cos(detvec_new[detct_new].Dec/DEGPRAD), sin(detvec_new[detct_new].Dec/DEGPRAD), detct_new);
      indexvec = {};
      status = kdrange_4d_index(kdvec, querypoint, matchrad/ASECPRAD, indexvec);
      if(indexvec.size()>0) {
	overlapnum++;
	for(i=0;i<long(indexvec.size());i++) sub_indexvec.push_back(kdvec[indexvec[i]].point.index);
      }
    } else {
      // We've finished one linkage and moved on to another one.
      linkagenum++;
      orig_linkagesize = newlinkage.size();
      max_obsnum = 0;
      enough_added = 1;
      // Calculate the timespan of the new linkage
      linkMJD={};
      for(i=0;i<orig_linkagesize;i++) linkMJD.push_back(newlinkage[i].MJD);
      sort(linkMJD.begin(),linkMJD.end());
      new_timespan = linkMJD[orig_linkagesize-1] - linkMJD[0];
      max_timespan = 0.0;
      cout << "Finished loading new linkage " << thisindex << " with timespan " << new_timespan << " days and " << orig_linkagesize << " detections, of which " << overlapnum << " overlapped with a submitted linkage\n";
      if(sub_indexvec.size()>0 && overlapnum<orig_linkagesize) {
	// Create a de-duplicated list of all the idstrings corresponding to duplicates
	IDvec={};
	ID_indexvec={};
	for(i=0;i<long(sub_indexvec.size());i++) {
	  int isnewID=1;
	  newstring = detvec_sub[sub_indexvec[i]].idstring;
	  //cout << "subindexvec, string: " << i << " " << sub_indexvec[i] << " " << detvec_sub[sub_indexvec[i]].idstring << " " << newstring << "\n";
	  for(j=0;j<long(IDvec.size());j++) {
	    if(IDvec[j]==newstring) isnewID=0;
	  }
	  if(isnewID==1) {
	    // No match was found to the new string. Add it to IDvec.
	    IDvec.push_back(newstring);
	    // And add the corresponding index to ID_indexvec.
	    ID_indexvec.push_back(sub_indexvec[i]);
	    //cout << "New ID found: " << IDvec.size()-1 << " " << IDvec[IDvec.size()-1] << " " << ID_indexvec[IDvec.size()-1] << "\n"; 
	  }
	}
	if(IDvec.size() != ID_indexvec.size()) {
	  cerr << "ERROR: IDvec and ID_indexvec do not have the same size: " << IDvec.size() << " vs. " << ID_indexvec.size() << "\n";
	  return(2);
	}
	cout << "The new linkage with index " << thisindex << " overlaps with " << IDvec.size() << " distinct previously submitted linkages\n";
	// Add all of these previously submitted linkages to the vector newlinkage,
	// and calculate the timespan and number of points for each of these overlapping, previously submitted linkages.
	for(i=0;i<long(IDvec.size());i++) {
	  //cout << "Checking overlapping linkage " << i << ", named " << IDvec[i] << ", at index " << ID_indexvec[i] << "\n";
	  linkMJD={};
	  sub_timespan = 0.0;
	  sub_obsnum = 0;
	  // Find the earliest entry for linkage IDvec[i] in detvec_sub.
	  detct = ID_indexvec[i];
	  newstring = detvec_sub[detct].idstring;
	  while(newstring == IDvec[i] && detct>0) {
	    detct--;
	    newstring = detvec_sub[detct].idstring;
	  }
	  if(detct>0 && newstring != IDvec[i]) {
	    // We've stepped off the start of the linkage. Increment detct to bring us back on
	    detct++;
	  }
	  newstring = detvec_sub[detct].idstring;
	  if(newstring == IDvec[i]) {
	    // We've successfully arrived at the beginning of this overlapping submitted linkage.
	    // Add this entry to the linkage
	    while(newstring == IDvec[i] && detct<detnum_sub) {
	      newlinkage.push_back(detvec_sub[detct]);
	      linkMJD.push_back(detvec_sub[detct].MJD);
	      detct++;
	      if(detct<detnum_sub) newstring = detvec_sub[detct].idstring;
	    }
	  } else {
	    cout << "Logical matching error: detct = " << detct << " newstring, detvec_sub, i, IDvec[i]: " << newstring << " " << detvec_sub[detct].idstring << " " << i << " " << IDvec[i] << "\n";
	    return(3);
	  }
	  // Calculate the timespan of this particular overlapping linkage
	  sort(linkMJD.begin(),linkMJD.end());
	  sub_obsnum = linkMJD.size();
	  sub_timespan = linkMJD[sub_obsnum-1] - linkMJD[0];
	  // See if this overlapping linkage has the largest number of points
	  // and/or the longest timespan of the previously submitted, overlapping linkages.
	  if(sub_timespan>max_timespan) max_timespan = sub_timespan;
	  if(sub_obsnum>max_obsnum) max_obsnum = sub_obsnum;
	  cout << "Overlapping linkage " << i << " (" << IDvec[i] << ") has timespan " << sub_timespan << " and " << sub_obsnum << " points, as compared to " << new_timespan << " and " << orig_linkagesize << " for the new linkage\n";
	  // Close loop over the previously submitted linkages that were found to overlap
	}
	merged_linkagesize = newlinkage.size();
	// Sort and then de-duplicate the merged linkage
	sort(newlinkage.begin(), newlinkage.end(), early_hldet());
	templinkage = newlinkage;
	newlinkage={};
	i=j=0;
	while(i<long(templinkage.size())) {
	  goodpt=j=i;
	  // See how many detections are redundant with detection i
	  while(j<long(templinkage.size()) && (fabs(templinkage[i].MJD-templinkage[j].MJD)*timescale + distradec01(templinkage[i].RA,templinkage[i].Dec,templinkage[j].RA,templinkage[j].Dec))<matchrad/3600.0) {
	    if(j!=goodpt && templinkage[goodpt].index<0 && templinkage[j].index==thisindex) {
	      // Currently, goodpt does not correspond to the new linkage,
	      // but the overlapping point j does. In order to prioritize new
	      // points, set goodpt to j
	      goodpt=j;
	    }
	    j++;
	  }
	  newlinkage.push_back(templinkage[goodpt]);
	  i=j;
	}
	final_linkagesize = newlinkage.size();
	cout << "Point counts for new linkage " << thisindex << " : " << orig_linkagesize << " " << merged_linkagesize << " " << final_linkagesize << "\n";
	// Calculate the timespan of the final linkage
	linkMJD={};
	for(i=0;i<final_linkagesize;i++) linkMJD.push_back(newlinkage[i].MJD);
	sort(linkMJD.begin(),linkMJD.end());
	merged_timespan = linkMJD[final_linkagesize-1] - linkMJD[0];
	cout << "Final merged linkage has timespan " << merged_timespan << " days, with " << final_linkagesize << " points, as compared max values of " << max_timespan << " and " << max_obsnum << " for linkages in the overlapping set\n";
	// Determine if enough points are added
	enough_added = 1;
	if(merged_timespan-max_timespan < mintimeadd || final_linkagesize-max_obsnum < minobsadd) {
	  cout << "The timespan of this merged linkage exceeds that of the longest constituent linkage by only " << merged_timespan-max_timespan << ",\n";
	  cout << "and has only " << final_linkagesize-max_obsnum << " additional points.\n";
	  cout << "As this does not meet the minimum values of " << mintimeadd << " and " << minobsadd << ", the linkage is rejected.\n";
	  enough_added = 0;
	}
	// Close case where some overlap is found
      } else if(overlapnum>=orig_linkagesize) {
	fully_redundant_ct++;
	cout << "New linkage " << thisindex << " was rejected because it was fully redundant with previously submitted linkages\n";
      }
      // Re-sort the linkage just in case
      sort(newlinkage.begin(), newlinkage.end(), early_hldet());
      // Make sure there are no time-duplicates
      istimedup=0;
      for(i=1;i<long(newlinkage.size());i++) {
	if(fabs(newlinkage[i].MJD-newlinkage[i-1].MJD)<IMAGETIMETOL/SOLARDAY) istimedup=1;
      }
      if(istimedup==0 && overlapnum<orig_linkagesize && enough_added==1) {
	outlinkagenum++;
	// This linkage is good. Set all the index values consistently
	for(i=0;i<long(newlinkage.size());i++) {
	  newlinkage[i].index = thisindex;
	}
	// Write it to the output file
	for(i=0;i<long(newlinkage.size());i++) {
	  outstream1 << fixed << setprecision(7) << newlinkage[i].MJD << "," << newlinkage[i].RA << "," << newlinkage[i].Dec << ",";
	  outstream1 << fixed << setprecision(4) << newlinkage[i].mag << ",";
	  outstream1 << fixed << setprecision(2) << newlinkage[i].trail_len << "," << newlinkage[i].trail_PA << ",";
	  outstream1 << fixed << setprecision(4) << newlinkage[i].sigmag << ",";
	  outstream1 << fixed << setprecision(4) << newlinkage[i].sig_across << "," << newlinkage[i].sig_along << ",";
	  outstream1 << newlinkage[i].image << "," << newlinkage[i].idstring << "," << newlinkage[i].band << ",";
	  outstream1 << newlinkage[i].obscode << "," << newlinkage[i].known_obj << ",";
	  outstream1 << newlinkage[i].det_qual << "," << newlinkage[i].index << "\n";
	  outdetnum++;
	}
	//Close case of no time duplicates
      } else if(istimedup>0 && enough_added==1) {
	cout << "New linkage " << thisindex << " was rejected because it contained points with duplicate times but different RA,Dec\n";
	for(i=0;i<long(newlinkage.size());i++) {
	  cout << fixed << setprecision(7) << newlinkage[i].MJD << "," << newlinkage[i].RA << "," << newlinkage[i].Dec << ",";
	  cout << fixed << setprecision(4) << newlinkage[i].mag << ",";
	  cout << fixed << setprecision(2) << newlinkage[i].trail_len << "," << newlinkage[i].trail_PA << ",";
	  cout << fixed << setprecision(4) << newlinkage[i].sigmag << ",";
	  cout << fixed << setprecision(4) << newlinkage[i].sig_across << "," << newlinkage[i].sig_along << ",";
	  cout << newlinkage[i].image << "," << newlinkage[i].idstring << "," << newlinkage[i].band << ",";
	  cout << newlinkage[i].obscode << "," << newlinkage[i].known_obj << ",";
	  cout << newlinkage[i].det_qual << "," << newlinkage[i].index << "\n";
	}
	timedupnum++;
      } else if(enough_added<1) {
	cout << "New linkage " << thisindex << " was rejected because it did not add enough to the best previously available linkage\n";
	non_augment_ct++;
      }
      if(detct_new<detnum_new) {
	// We've finished dealing with this linkage, but it isn't the last: start the next one.
	newlinkage = {}; // Wipe the old, finished linkage
	sub_indexvec = ID_indexvec = {};
	overlapnum = 0;
	newlinkage.push_back(detvec_new[detct_new]); // Start loading the new one.
	thisindex = detvec_new[detct_new].index; // Record the index of the new linkage
	// See if there is a duplicate of this point in the submissions file.
	querypoint = point4d_index((detvec_new[detct_new].MJD-mjdref)*timescale/DEGPRAD,cos(detvec_new[detct_new].RA/DEGPRAD)*cos(detvec_new[detct_new].Dec/DEGPRAD),sin(detvec_new[detct_new].RA/DEGPRAD)*cos(detvec_new[detct_new].Dec/DEGPRAD), sin(detvec_new[detct_new].Dec/DEGPRAD), detct_new);
	indexvec = {};
	status = kdrange_4d_index(kdvec, querypoint, matchrad/ASECPRAD, indexvec);
	if(indexvec.size()>0) {
	  overlapnum++;
	  for(i=0;i<long(indexvec.size());i++) sub_indexvec.push_back(kdvec[indexvec[i]].point.index);
	}
      }
      // Close case where we finished with a linkage.
    }
    detct_new++;
  }
  outstream1.close();
  cout << "Input file contained " << linkagenum << " new linkages totalling " << detnum_new << " individual detections\n";
  cout << "Of these, " << outlinkagenum << " linkages totaling " << outdetnum << " individual detections were written out,\n";
  cout << "while " << timedupnum << " linkages were rejected due to inconsistent duplication,\n";
  cout << fully_redundant_ct << " linkages were rejected because they were fully redundant with earlier submissions,\n";
  cout << "and " << non_augment_ct << " linkages were rejected because they did not sufficiently augment the best previously submitted linkage.\n";
    
  return(0);
}
