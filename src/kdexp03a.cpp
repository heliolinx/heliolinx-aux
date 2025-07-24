// July 22, 2025
// Test new kdnearest function.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

int main()
{
  string infile;
  int i=0;
  double x=0.0;
  double y=0.0;
  long pct=0;
  long medpt=0;
  xy_index xyi = xy_index(0.0,0.0,0);
  vector <xy_index> xyvec1 = {};
  vector <xy_index> left = {};
  vector <xy_index> right = {};
  int dim=1;
  kdpoint root = kdpoint(xyi,-1,-1,dim);
  kdpoint kdtest = kdpoint(xyi,-1,-1,dim);
  vector <kdpoint> kdvec ={};

  cout << "Enter the name of your input file\n";
  cin >> infile;

  // Read input file
  ifstream instream1 {infile};
  pct=0;
  if(!instream1) {
    cerr << "ERROR: can't open input file " << infile << "\n";
    return(1);
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    //Read a line with x, y
    instream1 >> x;
    instream1 >> y;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      xyi = xy_index(x,y,pct);
      xyvec1.push_back(xyi);
      pct++;
    }
  }

  cout << "Printing vector with original indices\n";
  for(i=0;i<long(xyvec1.size());i++) {
    cout << xyvec1[i].x << " " << xyvec1[i].y << " " << xyvec1[i].index << " " << i << "\n";
  }

  // sort(xyvec1.begin(), xyvec1.end(), xyind_lower_y());

  medpt = medindex(xyvec1,dim);
  cout << "\nMedian point: " << xyvec1[medpt].x << " " << xyvec1[medpt].y << " " << xyvec1[medpt].index << " " << medpt << "\n";

  root = kdpoint(xyvec1[medpt],-1,-1,1);
  kdvec.push_back(root);

  kdtest = kdvec[0];
  cout << "Root " << kdtest.point.x << " " << kdtest.point.y << " " << kdtest.point.index << "\n";
  cout << "About to create k-d tree\n";
  kdtree01(xyvec1,dim,medpt,0,kdvec);
  cout << "Created k-d tree of size " << kdvec.size() << "\n";
  
  for(i=0;i<long(kdvec.size());i++)
    {
      kdtest=kdvec[i];
      // cout << "Current point " << kdtest.point.x << " " << kdtest.point.y << " " << kdtest.point.index <<  " " << kdtest.left << " " << kdtest.right << " " << kdtest.dim << "\n";
    }
      
  dim=1;
  int notdone=1;
  while(notdone!=0)
    {
      double x,y;
      vector <long> indexvec;
      int bestpoint=0;
      x=y=0.0;
      cout << "Enter an x,y point\n";
      cin >> x;
      cin >> y;
      bestpoint = kdnearest01(kdvec,x,y);
      cout << "Best point is " << bestpoint << ": " << kdvec[bestpoint].point.x << " " << kdvec[bestpoint].point.y << " " << kdvec[bestpoint].point.index <<  " " << kdvec[bestpoint].left << " " << kdvec[bestpoint].right << " " << kdvec[bestpoint].dim << "\n";
      cout << "Enter 1 to go again, or 0 to quit\n";
      cin >> notdone;
    }
  
  return(0);  
}
