// June 05, 2023: NEO_disc_odds.cpp

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define PROBSTEP 0.01l
 
static void show_usage()
{
  cerr << "Usage: NEO_disc_odds -numobs number of obs per night -numnights total number of nights -numobsold number of obs needed for an old-style tracklet -numobsnew number of obs needed for a new tracklet -nightneed number of nights needed for a new-style discovery -outfile output file\n";
}

// binom_dist: June 06, 2023:
// Given probability prob of success on a single trial,
// use the binomial formula to calculate the probability of
// k successes in n trials.
double binom_dist(double prob, int n, int k)
{ 
  return(factorialD(n)/factorialD(k)/factorialD(n-k) * intpowD(prob,k) * intpowD(1.0l-prob,n-k));
}

int main(int argc, char *argv[])
{
  int numobs = 4;
  int numnights = 3;
  int numobsold = 4;
  int numobsnew = 2;
  int nightneed = 3;
  int i=0;
  double detprob=0.0l;
  double nightprobold=0.0l;
  double nightprobnew=0.0l;
  double nightprobnone=0.0l;
  double helioprob=0.0l;
  double convprob=0.0l;
  double noprob=0.0l;
  double tempprob=0.0l;
  double detsig=5.0l;
  ofstream outstream1;
  string outfile = "NEOprobtemp01.txt";
  double deltamag=0.0l;
  double deltasig=0.0l;
  double sigma=0.0l;
  double cumhelio=0.0l;
  double cumconv=0.0l;
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-numobs") {
      if(i+1 < argc) {
	//There is still something to read;
	numobs = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Number of obs per night keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-numnights") {
      if(i+1 < argc) {
	//There is still something to read;
	numnights = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Total number of nights keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-numobsold") {
      if(i+1 < argc) {
	//There is still something to read;
	numobsold = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Keyword for the number of observations needed for an old-style tracklet supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-numobsnew") {
      if(i+1 < argc) {
	//There is still something to read;
	numobsnew = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Keyword for the number of observations needed for a new tracklet supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-nightneed") {
      if(i+1 < argc) {
	//There is still something to read;
	nightneed = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Keyword for the number of nights needed for a new-style discovery supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-detsig") {
      if(i+1 < argc) {
	//There is still something to read;
	detsig = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Detection significance threshold keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile = argv[++i];
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


  cout << "Observations are obtained in sets of " << numobs << " per night, over a total of " << numnights << " observing nights\n";
  cout << "An object will be discovered conventionally if there are " << numobsold << " or more detections on any given night\n";
  cout << "An object will be discovered with heliolinc if there are at least " << numobsnew << " detections per night on at least " << nightneed << " distinct nights\n";

  outstream1.open(outfile);
  outstream1 << "#deltamag sigma deltasig detprob prob_convdet prob_heliodet prob_nodet\n";
  deltamag = -0.5l;
  while(deltamag<=0.5+PROBSTEP) {
    // Calculate the sigma offset corresponding to this delta-magnitude
    sigma = detsig*pow(10.0,0.4*deltamag);
    deltasig = sigma-detsig;
    if(deltasig<0.0l) {
      detprob = 0.5l - erf(-deltasig/sqrt(2.0l))/2.0l;
    } else if(deltasig==0.0l) {
      detprob=0.5l;
    } else if(deltasig>0.0l) {
      detprob = 0.5l + erf(deltasig/sqrt(2.0l))/2.0l;
    }
    // Calculate the probability of no detection on a single night
    nightprobnone=0.0l;
    for(i=0;i<numobsnew;i++) {
      nightprobnone += binom_dist(detprob, numobs, i);
    }
    // Calculate the probability of getting enough detections for heliolinc,
    // but not enough for a conventional tracklet, on a single night.
    nightprobnew=0.0l;
    for(i=numobsnew;i<numobsold;i++) {
      nightprobnew += binom_dist(detprob, numobs, i);
    }
    // Calculate the probability of getting enough detections for a
    // conventional tracklet on a single night
    nightprobold=0.0l;
    for(i=numobsold;i<=numobs;i++) {
      nightprobold += binom_dist(detprob, numobs, i);
    }
    cout << setprecision(6) << "Single night vals for detprob=" << detprob << ": " << nightprobnone << " + " << nightprobnew << " + " << nightprobold << " = " << nightprobnone + nightprobnew + nightprobold << "\n";

    // Calculate the probability of at least one conventional detection
    // over numnights
    convprob=0.0l;
    for(i=1;i<=numnights;i++) {
      convprob += binom_dist(nightprobold, numnights, i);
    }
    // Calculate the probability of at least nightneed new-style
    // detections without any conventional detections
    helioprob=0.0l;
    for(i=nightneed; i<=numnights; i++)
      {
	// Probability of NO old-style detections
	tempprob = binom_dist(nightprobold, numnights, 0);
	// Conditional probabability of i new-style detections  
	tempprob *= binom_dist(nightprobnew/(nightprobnew+nightprobnone), numnights, i);
	helioprob += tempprob;
      }
    // Calculate the probability of fewer than nightneed new-style
    // detections without any conventional detections
    noprob=0.0l;
    for(i=0; i<nightneed; i++)
      {
	// Probability of NO old-style detections
	tempprob = binom_dist(nightprobold, numnights, 0);
	// Conditional probabability of i new-style detections  
	tempprob *= binom_dist(nightprobnew/(nightprobnew+nightprobnone), numnights, i);
	noprob += tempprob;
      }
    outstream1 << setprecision(6) << deltamag << " " << sigma << " " << deltasig << " " << detprob << " " << convprob << " " << helioprob << " " << noprob << " " << convprob+helioprob+noprob << "\n";
    cumhelio+=helioprob;
    cumconv+=convprob;
    deltamag+=PROBSTEP;
  }
  outstream1.close();
  cout << "Heliolinc boost: " << cumhelio/cumconv << "\n";
  
  return(0);
}
