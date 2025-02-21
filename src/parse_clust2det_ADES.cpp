// December 12, 2024: parse_clust2det_ADES: Parse the output files from the
// heliolinc suite to produce a representation of the final linkages
// in MPC ADES format. to examine manually. Specifically, a hybrid file
// is produced with each cluster introduced by a summary line with
// 38 columns giving a comprehensive evaluation thereof, but then
// lines for each of the individual detections in that cluster follow,
// in the same format as lines from the input paired detection file.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: parse_clust2det_ADES -pairdet pairdet_file -insum input cluster summary file -clust2det input cluster-to-detection file -trackdiv tracklet_division_time -head header_file -out output file\n";
}

int main(int argc, char *argv[])
{
  string sumfile,clust2detfile,outfile,lnfromfile;
  vector <hlclust> inclustvec;
  vector  <longpair> inclust2det;
  vector <hldet> detvec;
  vector <hldet> cluster_detvec;
  vector <hldet> clustvec;
  vector <hldet> trackvec;
  double angvel,crosstrack,alongtrack,PA,poleRA,poleDec;
  angvel=crosstrack=alongtrack=PA=poleRA=poleDec=0.0;
  double arc,timespan;
  arc = timespan = 0;
  vector <double> angvelvec;
  vector <double> GCRvec;
  vector <double> PAvec;
  vector <double> arcvec;
  vector <double> timespanvec;
  vector <double> nightstepvec;
  vector <double> magvec;
  double min_nightstep,max_nightstep;
  long tracknum=0;
  double nightstep = 3.0l/24.0l;
  string pairdetfile,stest;
  ifstream instream1;
  ofstream outstream1;
  int status=0;
  long detnum=0;
  long detct=0;
  long clustct=0;
  long i=0;
  int verbose=0;
  long max_known_obj=0;
  double avg_det_qual=0.0l;
  double magrange,magmean,magrms,minGCR,maxGCR;
  magrange = magmean = magrms = minGCR = maxGCR = 0.0l;
  string stest1,stest2,header_file;
  int year,month;
  double day;
  string signstring;

  if(argc<9) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pairdets" || string(argv[i]) == "-pairdet" || string(argv[i]) == "-pd" || string(argv[i]) == "-pdet" || string(argv[i]) == "--pairdet" || string(argv[i]) == "--paireddetections" || string(argv[i]) == "--pairdetfile" || string(argv[i]) == "--pairdetections") {
      if(i+1 < argc) {
	//There is still something to read;
	pairdetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input pairdet file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-insum" || string(argv[i]) == "-inclust" || string(argv[i]) == "-clust" || string(argv[i]) == "-sum" || string(argv[i]) == "--input_summary" || string(argv[i]) == "--input_cluster" || string(argv[i]) == "--input_cluster_file" || string(argv[i]) == "--input_summary_file") {
      if(i+1 < argc) {
	//There is still something to read;
	sumfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster summary file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-trackdiv" || string(argv[i]) == "-nightstep") {
      if(i+1 < argc) {
	//There is still something to read;
	nightstep = stod(argv[++i]);
	nightstep /= 24.0l;
	i++;
      }
      else {
	cerr << "Input cluster list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-head" || string(argv[i]) == "-header" || string(argv[i]) == "--head" || string(argv[i]) == "--header" || string(argv[i]) == "--header_file") {
      if(i+1 < argc) {
	//There is still something to read;
	header_file=argv[++i];
	i++;
      }
      else {
	cerr << "Header file keyword supplied with no corresponding argument\n";
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
  if(pairdetfile.size()<=0) {
    cout << "\nERROR: input paired detection file is required\n";
    show_usage();
    return(1);
  } else if(sumfile.size()<=0) {
    cout << "\nERROR: input cluster summary file is required\n";
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

  cout << "input paired detection file " << pairdetfile << "\n";
  cout << "input cluster summary file " << sumfile << "\n";
  cout << "input cluster-to-detection file " << clust2detfile << "\n";
  cout << "output file " << outfile << "\n";


  // Read paired detection file
  detvec={};
  status=read_pairdet_file(pairdetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << detvec.size() << " data lines from paired detection file " << pairdetfile << "\n";
  
  // Read cluster summary file
  inclustvec={};
  status=read_clustersum_file(sumfile, inclustvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read input cluster summary file " << sumfile << "\n";
    cerr << "read_clustersum_file returned status = " << status << ".\n";
    return(1);
  }
  cout << "Read " << inclustvec.size() << " data lines from cluster summary file " << sumfile << "\n";
  inclust2det={};
  // Read cluster-to-detection file
  status=read_longpair_file(clust2detfile, inclust2det, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read cluster-to-detection file " << clust2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
    return(1);
  }
  cout << "Read " << inclust2det.size() << " data lines from cluster-to-detection file " << clust2detfile << "\n";

  status = parse_clust2det(detvec, inclust2det, cluster_detvec);
  if(status!=0) {
    cerr << "ERROR: parse_clust2det failed with error status " << status << "\n";
    return(status);
  }
  detnum = cluster_detvec.size();
  cout << "Wrote cluster detection vector with " << detnum << "entries\n";
  
  outstream1.open(outfile);
  // Write introductory ADES stuff
  outstream1 << "<?xml version=\'1.0\' encoding=\'UTF-8\'?>\n";
  outstream1 << "<ades version=\"2022\">\n";

  cout << "Writing " << inclustvec.size() << " clusters to output file " << outfile << "\n";
  detct=0;
  for(clustct=0 ; clustct<long(inclustvec.size()); clustct++) {
    // Echo header material to output file
    outstream1 << "  <obsBlock>\n";
    outstream1 << "    <obsContext>\n";
    instream1.open(header_file);
    if(!instream1) {
      cerr << "can't open input file " << header_file << "\n";
      return(1);
    }
    while (!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      getline(instream1,lnfromfile);
      if(lnfromfile.size()>0) outstream1 << lnfromfile << "\n";
    }
    instream1.close();

    // Load cluster-specific detection vector for this cluster
    clustvec={};
    if(detct<detnum) {
      if(cluster_detvec[detct].index!=clustct) {
	cerr << "ERROR: cluster counting mismatch, clustct = " << clustct << ", detct = " << detct << ", index = " << cluster_detvec[detct].index << "\n";
	return(2);
      }
      while(detct<detnum && cluster_detvec[detct].index==clustct) {
	clustvec.push_back(cluster_detvec[detct]);
	detct++;
      }
      cout << clustvec.size() << " points found for cluster " << clustct << "\n";
    }
    if(clustvec.size()>0) {
      // Time-sort clustvec, just to be sure
      sort(clustvec.begin(), clustvec.end(), early_hldet());
      // Run some analytics on this detection cluster
      avg_det_qual = 0.0l;
      max_known_obj=0;
      magvec={};
      for(i=0; i<long(clustvec.size()); i++) {
	avg_det_qual += double(clustvec[i].det_qual);
	if(clustvec[i].known_obj > max_known_obj) max_known_obj = clustvec[i].known_obj;
	if(clustvec[i].mag>0.0l) magvec.push_back(clustvec[i].mag);
      }
      avg_det_qual/=double(clustvec.size());
      // Loop over clustvec to extract individual tracklets
      trackvec={};
      nightstepvec={};
      trackvec.push_back(clustvec[0]);
      angvelvec = GCRvec = PAvec = timespanvec = arcvec = {};
      for(i=1; i<long(clustvec.size()); i++) {
	if((clustvec[i].MJD - clustvec[i-1].MJD) >= nightstep) {
	  // We're looking at a gap between two successive tracklets.
	  // We're interested in the distribution of such gaps,
	  // so we add it to the nightstepvec.
	  nightstepvec.push_back(clustvec[i].MJD - clustvec[i-1].MJD);
	}
	if((clustvec[i].MJD - clustvec[i-1].MJD) < nightstep) {
	  // Add a new point to this tracklet
	  trackvec.push_back(clustvec[i]);
	} else {
	  // A tracklet is finished. Analyze it.
	  tracknum = trackvec.size();
	  if(tracknum>1) {
	    greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
	    angvelvec.push_back(angvel);
	    if(tracknum>2) GCRvec.push_back(sqrt(DSQUARE(crosstrack)+DSQUARE(alongtrack)));
	    PAvec.push_back(PA);
	    timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
	    arc = timespan*angvel;
	    arcvec.push_back(arc*3600.0l);
	    timespanvec.push_back(timespan*24.0l);
	    // Wipe trackvec, and load the next point of the next tracklet.
	    trackvec = {};
	    trackvec.push_back(clustvec[i]);
	  } else if(tracknum==1) {
	    // The 'tracklet' is a singleton. Set all tracket vectors to error codes or zero.
	    angvelvec.push_back(-1.0l);
	    PAvec.push_back(-999.0l);
	    arcvec.push_back(0.0l);
	    timespanvec.push_back(0.0l);
	    // Wipe trackvec, and load the next point of the next tracklet.
	    trackvec = {};
	    trackvec.push_back(clustvec[i]);
	  }
	}
      }
      tracknum = trackvec.size();
      if(tracknum>1) {
	// Handle a final tracklet.
	greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
	angvelvec.push_back(angvel);
	if(tracknum>2) GCRvec.push_back(sqrt(DSQUARE(crosstrack)+DSQUARE(alongtrack)));
	PAvec.push_back(PA);
	timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
	arc = timespan*angvel;
	arcvec.push_back(arc*3600.0l);
	timespanvec.push_back(timespan*24.0l);
	// Wipe trackvec
	trackvec = {};
      } else if(tracknum==1) {
	// The final 'tracklet' is a singleton. Set all tracket vectors to error codes or zero.
	angvelvec.push_back(-1.0l);
	PAvec.push_back(-999.0l);
	arcvec.push_back(0.0l);
	timespanvec.push_back(0.0l);
	// Wipe trackvec.
	trackvec = {};
      }
      // Sort all of the tracklet statistics vectors
      tracknum = angvelvec.size();
      sort(angvelvec.begin(), angvelvec.end());
      sort(PAvec.begin(), PAvec.end());
      sort(timespanvec.begin(), timespanvec.end());
      sort(arcvec.begin(), arcvec.end());
      
      if(long(GCRvec.size())<1) minGCR = maxGCR = 0.0l;
      else if(long(GCRvec.size())==1) minGCR = maxGCR = GCRvec[0];
      else {
	sort(GCRvec.begin(), GCRvec.end());
	minGCR = GCRvec[0];
	maxGCR = GCRvec[GCRvec.size()-1];
      }

      // Sort nightstepvec
      sort(nightstepvec.begin(), nightstepvec.end());
      min_nightstep = max_nightstep = 0.0l;
      min_nightstep = nightstepvec[0];
      max_nightstep = nightstepvec[nightstepvec.size()-1];
      nightstepvec = {};

      // Analyze magnitude range.
      if(magvec.size()<=0) {
	magmean = 0.0;
	magrms = magrange = 99.9;
      } else if(magvec.size()==1) {
	magmean = magvec[0];
	magrms = magrange = 99.9;
      } else if(magvec.size()<=5) {
	// Sort magvec
	sort(magvec.begin(), magvec.end());
	dmeanrms01(magvec, &magmean, &magrms);
	// Magrange will be the full max-min
	magrange = magvec[magvec.size()-1] - magvec[0];
      } else {
	// Sort magvec
	sort(magvec.begin(), magvec.end());
	dmeanrms01(magvec, &magmean, &magrms);
	magrange = magvec[magvec.size()-2] - magvec[1];
	// Magrange will be the second-largest value minus
	// the second-smallest, offering some robustness
	// against outliers
      }
      magvec={};
      outstream1 << "      <comment>\n";
      outstream1  << fixed << setprecision(3) << "        <line> clusternum: " << inclustvec[clustct].clusternum << "  posRMS: " << inclustvec[clustct].posRMS << "  velRMS: " << inclustvec[clustct].velRMS << "  totRMS: " << inclustvec[clustct].totRMS << " </line>\n";
      outstream1 << fixed << setprecision(4)  << "        <line> astromRMS: " << inclustvec[clustct].astromRMS << "  timespan: " << inclustvec[clustct].timespan << "  uniquepoints: " << inclustvec[clustct].uniquepoints << "  obsnights: " << inclustvec[clustct].obsnights << " </line>\n";
      outstream1 << "        <line> metric: " << inclustvec[clustct].metric;
      outstream1 << fixed << setprecision(6) << "  orbit_a: " << inclustvec[clustct].orbit_a << "  orbit_e: " << inclustvec[clustct].orbit_e << "  orbit_MJD: " << inclustvec[clustct].orbit_MJD << " </line>\n";
      outstream1 << fixed << setprecision(2) << "        <line> orbitX: " << inclustvec[clustct].orbitX << "  orbitY: " << inclustvec[clustct].orbitY  << "  orbitZ: " << inclustvec[clustct].orbitZ << " </line>\n";
      outstream1 << fixed << setprecision(6) << "        <line> orbitVX: " << inclustvec[clustct].orbitVX << "  orbitVY: " << inclustvec[clustct].orbitVY << "  orbitVZ: " << inclustvec[clustct].orbitVZ << "  orbit_eval_count: " << inclustvec[clustct].orbit_eval_count << " </line>\n";
      outstream1 << fixed << setprecision(1)  << "        <line> avg_det_qual: " << avg_det_qual << "  max_known_obj: " << max_known_obj;
      outstream1 << fixed << setprecision(3) << "  minvel: " << angvelvec[0] << "  maxvel: " << angvelvec[tracknum-1] << " </line>\n";
      outstream1 << "        <line> minGCR: " << minGCR << "  maxGCR: " << maxGCR;
      outstream1 << fixed << setprecision(1) << "  minpa: " << PAvec[0] << "  maxpa: " << PAvec[tracknum-1] << " </line>\n";
      outstream1 << fixed << setprecision(3) << "        <line> mintimespan: " << timespanvec[0] << "  maxtimespan: " << timespanvec[tracknum-1] << "  minarc: " << arcvec[0] << "  maxarc: " << arcvec[tracknum-1] << " </line>\n";
      outstream1 << "        <line> stringID: " << clustvec[0].idstring << "  min_nightstep: " << min_nightstep << "  max_nightstep: " << max_nightstep << "  rating: " << inclustvec[clustct].rating << " </line>\n";
      outstream1 << "        <line> magmean: " << magmean << "  magrms: " << magrms << "  magrange: " << magrange << " </line>\n";
      outstream1 << "      </comment>\n";
      outstream1 << "    </obsContext>\n";
      outstream1 << "    <obsData>\n";
      // Write this data to the output file
      for(i=0; i<long(clustvec.size()); i++) {
	outstream1 << "      <optical>\n";
	// The temporary ID will just be 'a' followed by a 6-digit number
	if(clustvec[i].index<=999999) stest1 = to_string(clustvec[i].index);
	else {
	  stest1 = to_string(clustvec[i].index%1000000);
	}
	while(stest1.size()<6) {
	  stest1 = "0" + stest1;
	}
	outstream1 << "        <trkSub>a" << stest1 << "</trkSub>\n";
	outstream1 << "        <obsSubID>" << clustvec[i].idstring << "</obsSubID>\n";
	outstream1 << "        <mode>CCD</mode>\n";
	outstream1 << "        <stn>" << clustvec[i].obscode << "</stn>\n";
	// Convert the input MJD, assumed to be TAI, into a calendar date and time
	// in UTC, as required by the MPC.
	mjd2mpcdate(clustvec[i].MJD - (TAI_MINUS_UTC/SOLARDAY), year, month, day);
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
	outstream1  << fixed << setprecision(6) << "        <ra>" << clustvec[i].RA << "</ra>\n";
	outstream1  << fixed << setprecision(6) << "        <dec>" << clustvec[i].Dec << "</dec>\n";
	outstream1  << "        <astCat>Gaia2</astCat>\n";
	if(clustvec[i].sig_across>0 && clustvec[i].sig_along>0) {
	  // outstream1  << fixed << setprecision(3) << "        <rmsRA>" << clustvec[i].sig_across << "</rmsRA>\n";
	  // outstream1  << fixed << setprecision(3) << "        <rmsDec>" << clustvec[i].sig_along << "</rmsDec>\n";
	}
	outstream1  << fixed << setprecision(3) << "        <mag>" << clustvec[i].mag << "</mag>\n";
	// outstream1  << fixed << setprecision(3) << "        <rmsMag>" << clustvec[i].sigmag << "</rmsMag>\n";
	outstream1  << "        <band>L" << clustvec[i].band << "</band>\n";
	outstream1 << "      </optical>\n";
      }
      outstream1 << "    </obsData>\n";
    }
    outstream1 << "  </obsBlock>\n";
  }
  outstream1 << "</ades>\n";
  outstream1.close();
  return(0);
}
