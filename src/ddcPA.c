#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PVNUM 22
#define NML 256
#define MAXLINE 300
#define PI (double)3.141592653589793238
#define ASECPRAD (double)206264.806247 /*Arcseconds per radian*/
#define DEGPRAD ((double)180.0/PI) /*Degrees per radian*/
#define MAXMEM 40000000000.0 /*Refuse to allocate a matrix or vector that will
                             use more than 40GB of memory.*/
#define DSQUARE(x) ((double)(x)*(double)(x))
#define FSQKILL(fsq,tq) {printf("ERROR: fscanf read %d quantities when %d were required\n",fsq,tq); return(2);}
#define SQKILL(sq,tq) {printf("ERROR: scanf read %d quantities when %d were required\n",sq,tq); return(3);}


/*February 13, 2024: ddcPA.c: Closely based on ddc_pacalc03.c,
but does not use any routines from the bezalel01 library, for
greater portability and ease of compilation.

Suggested compile command:
gcc -std=c99 ddcPA.c -L. -lm -o ddcPA

Description of ancestor program ddc_pacalc03.c:
Like ddc_pacalc01.c, but reads the header as an independent
ASCII text file, rather than as an actual FITS header.

Description of ancestor program ddc_pacalc01.c:
Given an input ATLAS image, and the corresponding ddc file,
read the WCS information from the image header, and then, for
each line in the ddc file, calculate the celestial position
angle (degrees east from north) of each trail, using the
phi angle (measured CCW from the x-axis) recorded in the ddc
file. Also calculate the length of the trail in arcseconds.
Write an output file that includes the full ddc line, plus
the two additional at the end: celestial position angle
in degrees east from north, and trail length in arcseconds.*/

static void show_usage()
{
  printf("ddcPA -inheader image_name -inddc ddc_filename -outfile output_filename -verbose verbosity\n");
}

int readPVtextheader01(char *imname,double *CRvec,double **pvmat);
int poleswitch01(double *initpos,double *polepos, double *newpos,double oldpolera);
double *dvec(int np);
double **dmat(int nx,int ny);
int free_dmat(double **mt,int nx,int ny);
int linecount(char filename[]);
void readline03(FILE *fp1,int maxline,char *linehold);
int distradec02(double *ra,double *dec,double *dist,double *pa);

int main(int argc,char *argv[])
{
  double **pvmat,*CRvec;
  char headname[NML],infile[NML],outfile[NML],linehold[MAXLINE+1];
  strcpy(headname,"");
  strcpy(infile,"");
  strcpy(outfile,"test.ddcpa");
  int sq,status,ddclinenum,datalinenum,i,j;
  double xa,ya,xb,yb;
  FILE *fp1,*fp2;
  float d1,d2,d3,d4,d5;
  double major,phi,traillen,celpa;
  double dRA_dx,dRA_dy,dDec_dx,dDec_dy;
  double dRA_dxb,dRA_dyb,dDec_dxb,dDec_dyb;
  double dxb_dxa,dxb_dya,dyb_dxa,dyb_dya;
  double ddist_dRA,ddist_dDec,dpa_dRA,dpa_dDec;
  double ddist_dx,ddist_dy,dpa_dx,dpa_dy;
  double dRA,dDec,dx,dy;
  double projRA,projDec,pa,dist;
  double *initpos,*polepos,*newpos,oldpolera;
  double trueRA,trueDec;
  double dpa,ddist,rpa,cpa,cdist;
  double *racomp,*deccomp;
  int ddcheadnum=0;
  int verbose=0;
  
  initpos = dvec(2);
  polepos = dvec(2);
  newpos = dvec(2);
  racomp = dvec(2);
  deccomp = dvec(2);

  sq = status = ddclinenum = datalinenum = i = j = 0;
  
  pvmat = dmat(2,PVNUM);
  CRvec = dvec(8);  /*order is CRVAL1, CRVAL2, CRPIX1, CRPIX2,
                           CD1_1, CD1_2, CD2_1, and CD2_2.*/

  /* Parse the arguments */
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-inheader") == 0 || strcmp(argv[i], "-header") == 0 || strcmp(argv[i], "-inhead") == 0 || strcmp(argv[i], "-inimage") == 0 || strcmp(argv[i], "-inimg") == 0 || strcmp(argv[i], "-img") == 0) {
      strcpy(headname, argv[++i]);
    } else if (strcmp(argv[i], "-inddc") == 0 || strcmp(argv[i], "-infile") == 0 || strcmp(argv[i], "-ddc") == 0) {
      strcpy(infile, argv[++i]);
    } else if (strcmp(argv[i], "-outfile") == 0 || strcmp(argv[i], "-out") == 0 || strcmp(argv[i], "-outpa") == 0) {
      strcpy(outfile, argv[++i]);
    } else if (strcmp(argv[i], "-verbose") == 0 || strcmp(argv[i], "-verb") == 0 || strcmp(argv[i], "-VERB") == 0) {
      sq = sscanf(argv[++i],"%d",&verbose);
    } else {
      printf("ERROR: unrecognized input argument %s\n", argv[i]);
      show_usage();
      return(1);
    }
  }

  if(strlen(headname)<=0) {
    printf("ERROR: No valid name for input image/header file\n");
    show_usage();
    return(1);
  }
  if(strlen(infile)<=0) {
    printf("ERROR: No valid name for input ddc file\n");
    show_usage();
    return(1);
  }
  if(strcmp(outfile,"test.ddcpa")==0) {
    printf("Warning: no name for the output file has been specified.\n");
    printf("Default name test.ddcpa will be used\n");
  }

  status = readPVtextheader01(headname,CRvec,pvmat);
  if(status!=0) {
    printf("ERROR: readPV01 returns error status %d\n",status);
    return(status);
  }

  if((fp1 = fopen(infile,"r")) != NULL) fclose(fp1);
  else {
    printf("ERROR: unable to read input ddc file %s\n", infile);
    return(2);
  }
  ddclinenum = linecount(infile);
  /*Count the number of header lines*/
  fp1 = fopen(infile,"r");
  ddcheadnum=i=0;
  readline03(fp1,MAXLINE,linehold);
  i++;
  while(i<ddclinenum && linehold[0]=='#') {
    ddcheadnum++;
    readline03(fp1,MAXLINE,linehold);
    i++;
  }
  fclose(fp1);
  datalinenum = ddclinenum-ddcheadnum;
  printf("Your input file appears to have %d lines\n",ddclinenum);
  printf("including %d header lines\n",ddcheadnum);
  fp1 = fopen(infile,"r");
  fp2 = fopen(outfile,"w");
  /*Copy ddc header, up to the column header line*/
  for(i=1;i<ddcheadnum;i++) {
    strcpy(linehold,"");
    readline03(fp1,MAXLINE,linehold);
    for(j=0;j<=MAXLINE;j++) {
      if(linehold[j]=='\n') linehold[j]=0;
    }
    fprintf(fp2,"%s\n",linehold);
  }
  /*Copy column header line with new entries*/
  strcpy(linehold,"");
  readline03(fp1,MAXLINE,linehold);
  for(j=0;j<=MAXLINE;j++) {
    if(linehold[j]=='\n') linehold[j]=0;
  }
  fprintf(fp2,"%s   celpa     trailarc\n",linehold);
  for(i=1;i<=datalinenum;i++) {
       strcpy(linehold,"");
    readline03(fp1,MAXLINE,linehold);
    for(j=0;j<=MAXLINE;j++) {
      if(linehold[j]=='\n') linehold[j]=0;
    }
    /*Read x, y, major axis, and phi for this ddc line*/
    sq = sscanf(linehold,"%f %f %f %f %lf %lf %lf %f %lf", &d1, &d2, &d3, &d4, &xa, &ya, &major, &d5, &phi);
    if(sq!=9) SQKILL(sq,9);

    /*first apply CD matrix*/
    xb = (xa-CRvec[3])*CRvec[5] + (ya-CRvec[4])*CRvec[6];
    yb = (xa-CRvec[3])*CRvec[7] + (ya-CRvec[4])*CRvec[8];

    /*Calculate the projected RA*/
    projRA = pvmat[1][0] + pvmat[1][1]*xb + pvmat[1][2]*yb + pvmat[1][4]*xb*xb + pvmat[1][5]*xb*yb + pvmat[1][6]*yb*yb + pvmat[1][7]*xb*xb*xb + pvmat[1][8]*xb*xb*yb + pvmat[1][9]*xb*yb*yb + pvmat[1][10]*yb*yb*yb + pvmat[1][12]*xb*xb*xb*xb + pvmat[1][13]*xb*xb*xb*yb + pvmat[1][14]*xb*xb*yb*yb + pvmat[1][15]*xb*yb*yb*yb + pvmat[1][16]*yb*yb*yb*yb + pvmat[1][17]*xb*xb*xb*xb*xb + pvmat[1][18]*xb*xb*xb*xb*yb + pvmat[1][19]*xb*xb*xb*yb*yb + pvmat[1][20]*xb*xb*yb*yb*yb + pvmat[1][21]*xb*yb*yb*yb*yb + pvmat[1][22]*yb*yb*yb*yb*yb;

  /*Calculate the projected DEC*/
    projDec = pvmat[2][0] + pvmat[2][1]*yb + pvmat[2][2]*xb + pvmat[2][4]*yb*yb + pvmat[2][5]*yb*xb + pvmat[2][6]*xb*xb + pvmat[2][7]*yb*yb*yb + pvmat[2][8]*yb*yb*xb + pvmat[2][9]*yb*xb*xb + pvmat[2][10]*xb*xb*xb + pvmat[2][12]*yb*yb*yb*yb + pvmat[2][13]*yb*yb*yb*xb + pvmat[2][14]*yb*yb*xb*xb + pvmat[2][15]*yb*xb*xb*xb + pvmat[2][16]*xb*xb*xb*xb + pvmat[2][17]*yb*yb*yb*yb*yb + pvmat[2][18]*yb*yb*yb*yb*xb + pvmat[2][19]*yb*yb*yb*xb*xb + pvmat[2][20]*yb*yb*xb*xb*xb + pvmat[2][21]*yb*xb*xb*xb*xb + pvmat[2][22]*xb*xb*xb*xb*xb;

    /*Transform back to distance and PA relative to CRval coords*/
    /*We are trying to invert the following transformation:
      projRA[i] = (double)180.0/PI * tan(dist) * (-sin(pa));
      projDec[i] = (double)180.0/PI * tan(dist) * (cos(pa));
    */

    if(verbose>0) printf("projected RA and DEC: %f\t%f\n",projRA,projDec);
    /*First extract a distance that will really be tan(dist)*/
    dist = pow(pow(projRA,2.0)+pow(projDec,2.0),0.5);
    if(verbose>0) printf("tan(dist) = %f degrees\n",dist);

    /*Now take the arctangent.*/
    dist = atan(dist*PI/(double)180.0);

    /*Now calculate the pa*/
    if(projRA>=0.0) pa = PI/(double)2.0 - atan(projDec/projRA);
    else pa = (double)3.0*PI/(double)2.0 - atan(projDec/projRA);
    if(verbose>0) printf("PA = %f degrees\n",pa*180.0/PI);
    cpa = pa*DEGPRAD;
    /*Now reverse the chirality of PA to make it equivalent to RA
      in a coordinate system with the pole defined by CRval*/
    pa = (double)2.0*PI - pa;
    initpos[1] = pa;
    initpos[2] = PI/(double)2.0 - dist;
    polepos[1] = (double)0.0;
    polepos[2] = CRvec[2]*PI/(double)180.0;
    oldpolera = CRvec[1]*PI/(double)180.0;
    poleswitch01(initpos,polepos,newpos,oldpolera);

    /*Now all we need to do is output newpos in degrees*/
    trueRA = newpos[1]*(double)180.0/PI;
    trueDec = newpos[2]*(double)180.0/PI;
    if(verbose>0) printf("True coords: %f %f\n", trueRA,trueDec);

    /*Calculate derivatives of the WCS!*/
    /*derivatives of the CD matrix*/
    dxb_dxa = CRvec[5];
    dxb_dya = CRvec[6];
    dyb_dxa = CRvec[7];
    dyb_dya = CRvec[8];

    /*Now derivatives of projected RA*/
    
    dRA_dxb = pvmat[1][1] + 2.0*pvmat[1][4]*xb + pvmat[1][5]*yb + 3.0*pvmat[1][7]*xb*xb + 2.0*pvmat[1][8]*xb*yb + pvmat[1][9]*yb*yb + 4.0*pvmat[1][12]*xb*xb*xb + 3.0*pvmat[1][13]*xb*xb*yb + 2.0*pvmat[1][14]*xb*yb*yb + pvmat[1][15]*yb*yb*yb + 5.0*pvmat[1][17]*xb*xb*xb*xb + 4.0*pvmat[1][18]*xb*xb*xb*yb + 3.0*pvmat[1][19]*xb*xb*yb*yb + 2.0*pvmat[1][20]*xb*yb*yb*yb + pvmat[1][21]*yb*yb*yb*yb;

    dRA_dyb = pvmat[1][2] + pvmat[1][5]*xb + 2.0*pvmat[1][6]*yb + pvmat[1][8]*xb*xb + 2.0*pvmat[1][9]*xb*yb + 3.0*pvmat[1][10]*yb*yb + pvmat[1][13]*xb*xb*xb + 2.0*pvmat[1][14]*xb*xb*yb + 3.0*pvmat[1][15]*xb*yb*yb + 4.0*pvmat[1][16]*yb*yb*yb + pvmat[1][18]*xb*xb*xb*xb + 2.0*pvmat[1][19]*xb*xb*xb*yb + 3.0*pvmat[1][20]*xb*xb*yb*yb + 4.0*pvmat[1][21]*xb*yb*yb*yb + 5.0*pvmat[1][22]*yb*yb*yb*yb;

    dDec_dxb = pvmat[2][2] + pvmat[2][5]*yb + 2.0*pvmat[2][6]*xb + pvmat[2][8]*yb*yb + 2.0*pvmat[2][9]*yb*xb + 3.0*pvmat[2][10]*xb*xb + pvmat[2][13]*yb*yb*yb + 2.0*pvmat[2][14]*yb*yb*xb + 3.0*pvmat[2][15]*yb*xb*xb + 4.0*pvmat[2][16]*xb*xb*xb + pvmat[2][18]*yb*yb*yb*yb + 2.0*pvmat[2][19]*yb*yb*yb*xb + 3.0*pvmat[2][20]*yb*yb*xb*xb + 4.0*pvmat[2][21]*yb*xb*xb*xb + 5.0*pvmat[2][22]*xb*xb*xb*xb;

    dDec_dyb = pvmat[2][1] + 2.0*pvmat[2][4]*yb + pvmat[2][5]*xb + 3.0*pvmat[2][7]*yb*yb + 2.0*pvmat[2][8]*yb*xb + pvmat[2][9]*xb*xb + 4.0*pvmat[2][12]*yb*yb*yb + 3.0*pvmat[2][13]*yb*yb*xb + 2.0*pvmat[2][14]*yb*xb*xb + pvmat[2][15]*xb*xb*xb + 5.0*pvmat[2][17]*yb*yb*yb*yb + 4.0*pvmat[2][18]*yb*yb*yb*xb + 3.0*pvmat[2][19]*yb*yb*xb*xb + 2.0*pvmat[2][20]*yb*xb*xb*xb + pvmat[2][21]*xb*xb*xb*xb;

    /*Now multiply through by the derivatives of the CD matrix*/
    dRA_dx = dRA_dxb * dxb_dxa  +  dRA_dyb * dyb_dxa;
    dRA_dy = dRA_dxb * dxb_dya  +  dRA_dyb * dyb_dya;

    dDec_dx = dDec_dxb * dxb_dxa  +  dDec_dyb * dyb_dxa;
    dDec_dy = dDec_dxb * dxb_dya  +  dDec_dyb * dyb_dya;

    /*Convert the previous derivatives from degrees to radians*/
    dRA_dx /= DEGPRAD;
    dRA_dy /= DEGPRAD;
    dDec_dx /= DEGPRAD;
    dDec_dy /= DEGPRAD;

    /*Use these derivatives to get a local pixel scale in the x and y directions*/
    if(verbose>0) printf("Pixel scale in x dimension is %f\n", sqrt(DSQUARE(dRA_dx)+DSQUARE(dDec_dx))*ASECPRAD);
    if(verbose>0) printf("Pixel scale in y dimension is %f\n", sqrt(DSQUARE(dRA_dy)+DSQUARE(dDec_dy))*ASECPRAD);
    
    /*But these are really the derivatives of the Gnomonically projected
      RA and Dec.*/

    /*Calculate the derivatives of distance with respect to projected RA and Dec*/
    projRA/=DEGPRAD;
    projDec/=DEGPRAD;
    ddist_dRA = projRA/(sqrt(projRA*projRA + projDec*projDec)*(1.0 + projRA*projRA + projDec*projDec));
    ddist_dDec = projDec/(sqrt(projRA*projRA + projDec*projDec)*(1.0 + projRA*projRA + projDec*projDec));

    dpa_dRA = projDec/(projRA*projRA + projDec*projDec);
    dpa_dDec = -projRA/(projRA*projRA + projDec*projDec);


    /*Calculate ddist_dx and ddist_dy*/
    ddist_dx = ddist_dRA * dRA_dx  +  ddist_dDec * dDec_dx;
    ddist_dy = ddist_dRA * dRA_dy  +  ddist_dDec * dDec_dy;
    /*Calculate dpa_dx and dpa_dy*/
    dpa_dx = dpa_dRA * dRA_dx  +  dpa_dDec * dDec_dx;
    dpa_dy = dpa_dRA * dRA_dy  +  dpa_dDec * dDec_dy;
    //printf("ddist derivatives: %f %f dpa derivatives: %f %f\n", ddist_dx*ASECPRAD, ddist_dy*ASECPRAD, dpa_dx*sin(dist)*ASECPRAD, dpa_dy*sin(dist)*ASECPRAD);

    /*Now use the phi angle to get dx and dy*/
    dx = cos(phi/DEGPRAD);
    dy = sin(phi/DEGPRAD);

    /*Now use these to get the derivatives of dist and pa along the trail*/
    ddist = ddist_dx * dx  +  ddist_dy * dy;
    dpa = sin(dist)*dpa_dx * dx  +  sin(dist)*dpa_dy * dy;

    /*Now ddist and dpa should both be true arc angles on the sky, in radians,
      corresponding to a movement of one pixel along the trail*/
    traillen = sqrt(ddist*ddist + dpa*dpa)*major*ASECPRAD;

    /*Calculate celestial position angle from source position toward image center*/
    racomp[1] = trueRA;
    racomp[2] = CRvec[1];
    deccomp[1] = trueDec;
    deccomp[2] = CRvec[2];
    distradec02(racomp,deccomp,&cdist,&cpa);
    cpa-=180.0; /*Reverse the angle for continuity with dpa and rpa*/
    if(cpa<0.0) cpa+=360.0;

    if(dpa==0.0 && ddist>=0.0) rpa=0.0;
    else if(dpa==0.0 && ddist<0.0) rpa=180.0;
    else if(dpa>0.0) rpa = 90.0 - DEGPRAD*atan(ddist/dpa);
    else if(dpa<0.0) rpa = 270.0 - DEGPRAD*atan(ddist/dpa); 
    celpa = cpa + rpa;
    while(celpa>=360.0) celpa -= 360.0;
    while(celpa<0.0) celpa += 360.0;
    fprintf(fp2,"%s   %.2f %.2f\n",linehold,celpa,traillen);
  }

    
  fclose(fp1);
  fclose(fp2);

  free(CRvec);
  free_dmat(pvmat,2,PVNUM);
  free(racomp);
  free(deccomp);
  
  return(0);
}

/*dvec(): May 01, 2019: Duplicate functionality of Numerical Recipes dvector() function
without using copyrighted code.*/
double *dvec(int np)
{
  double *vt;
  float memuse;

  memuse=(float)np*(float)sizeof(double);
  if(np<1)
    {
      fprintf(stderr,"Error in dvec(): called with size %d\n",np);
      return(NULL);
    }
  else if(memuse>MAXMEM)
    {
      fprintf(stderr,"Error in dvec(): allocated vector would be too big (%d = %.0f bytes is greater than %.0f\n",np,memuse,MAXMEM);
      return(NULL);
    }
  vt = (double*) calloc(1+np,sizeof(double));
  return(vt);
}

/*dmat(): May 01, 2019: Duplicate functionality of Numerical Recipes dmatrix() function
without using copyrighted code.*/
double **dmat(int nx,int ny)
{
  double **mt;
  float memuse;
  int i;

  memuse=(float)nx*(float)ny*(float)sizeof(double);
  if(nx<1 || ny<1)
    {
      fprintf(stderr,"Error in mat(): called with dimensions %d %d\n",nx,ny);
      return(NULL);
    }
  else if(memuse>MAXMEM)
    {
      fprintf(stderr,"Error in mat(): allocated matrix would be too big (%dx%d = %.0f bytes\n",nx,ny,memuse);
      return(NULL);
    }

  mt = (double**)calloc(nx+1,sizeof(double*));
  for(i=0;i<=nx;i++) 
    {
      mt[i]=(double*)calloc(ny+1,sizeof(double));
      if(mt[i] == NULL) return(NULL);
    }
  return(mt);
}

/*free_dmat(): May 01, 2019: Free a matrix allocated with dmat().*/
int free_dmat(double **mt,int nx,int ny)
{ 
  int i;
  for(i=0;i<=nx;i++) free(mt[i]);
  free(mt);
  return(0);
}

/*readPVtextheader01: February 07, 2024: read the PV matrix
keywords out of the ASCII version of the header of an ATLAS image.*/
int readPVtextheader01(char *imname,double *CRvec,double **pvmat)
{
  int status;
  FILE *fp1;
  double pv1_0,pv1_1,pv1_2,pv1_4,pv1_5,pv1_6,pv1_7,pv1_8,pv1_9,pv1_10,pv1_12;
  double pv1_13,pv1_14,pv1_15,pv1_16,pv1_17,pv1_18,pv1_19,pv1_20,pv1_21,pv1_22;
  double pv2_0,pv2_1,pv2_2,pv2_4,pv2_5,pv2_6,pv2_7,pv2_8,pv2_9,pv2_10,pv2_12;
  double pv2_13,pv2_14,pv2_15,pv2_16,pv2_17,pv2_18,pv2_19,pv2_20,pv2_21,pv2_22;
  int linenum,linect,sq,readcount;
  char linehold[NML+1],sjunk1[NML],sjunk2[NML];
  
  if((fp1 = fopen(imname,"r")) != NULL) fclose(fp1);
  else {
    printf("ERROR: could not open header file %s\n", imname);
    return(2);
  }
  linenum = linecount(imname);
  fp1 = fopen(imname,"r");
  status=readcount=0;
  for(linect=1;linect<=linenum;linect++) {
    strcpy(linehold,"");
    readline03(fp1,NML,linehold);
    if(strncmp("CRVAL1  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+1);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("CRVAL2  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+2);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("CRPIX1  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+3);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("CRPIX2  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+4);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("CD1_1   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+5);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("CD1_2   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+6);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("CD2_1   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+7);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("CD2_2   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,CRvec+8);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_0   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_0);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_1   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_1);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_2   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_2);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_4   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_4);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_5   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_5);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_6   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_6);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_7   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_7);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_8   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_8);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_9   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_9);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_10  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_10);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_12  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_12);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_13  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_13);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_14  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_14);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_15  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_15);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_16  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_16);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_17  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_17);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_18  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_18);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_19  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_19);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_20  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_20);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_21  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_21);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV1_22  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv1_22);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_0   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_0);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_1   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_1);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_2   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_2);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_4   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_4);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_5   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_5);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_6   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_6);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_7   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_7);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_8   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_8);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_9   =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_9);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_10  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_10);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_12  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_12);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_13  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_13);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_14  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_14);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_15  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_15);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_16  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_16);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_17  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_17);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_18  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_18);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_19  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_19);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_20  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_20);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_21  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_21);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    } else if(strncmp("PV2_22  =",linehold,9)==0) {
      sq=sscanf(linehold,"%s %s %lf",sjunk1,sjunk2,&pv2_22);
      if(sq!=3) FSQKILL(sq,3);
      readcount++;
    }
  }
  if(readcount!=50) {
    printf("ERROR: read only %d of 50 required quantities from header file %s\n",readcount,imname);
    return(3);
  }
  pvmat[1][0] = pv1_0;
  pvmat[1][1] = pv1_1;
  pvmat[1][2] = pv1_2;
  pvmat[1][4] = pv1_4;
  pvmat[1][5] = pv1_5;
  pvmat[1][6] = pv1_6;
  pvmat[1][7] = pv1_7;
  pvmat[1][8] = pv1_8;
  pvmat[1][9] = pv1_9;
  pvmat[1][10] = pv1_10;
  pvmat[1][12] = pv1_12;
  pvmat[1][13] = pv1_13;
  pvmat[1][14] = pv1_14;
  pvmat[1][15] = pv1_15;
  pvmat[1][16] = pv1_16;
  pvmat[1][17] = pv1_17;
  pvmat[1][18] = pv1_18;
  pvmat[1][19] = pv1_19;
  pvmat[1][20] = pv1_20;
  pvmat[1][21] = pv1_21;
  pvmat[1][22] = pv1_22;

  pvmat[2][0] = pv2_0;
  pvmat[2][1] = pv2_1;
  pvmat[2][2] = pv2_2;
  pvmat[2][4] = pv2_4;
  pvmat[2][5] = pv2_5;
  pvmat[2][6] = pv2_6;
  pvmat[2][7] = pv2_7;
  pvmat[2][8] = pv2_8;
  pvmat[2][9] = pv2_9;
  pvmat[2][10] = pv2_10;
  pvmat[2][12] = pv2_12;
  pvmat[2][13] = pv2_13;
  pvmat[2][14] = pv2_14;
  pvmat[2][15] = pv2_15;
  pvmat[2][16] = pv2_16;
  pvmat[2][17] = pv2_17;
  pvmat[2][18] = pv2_18;
  pvmat[2][19] = pv2_19;
  pvmat[2][20] = pv2_20;
  pvmat[2][21] = pv2_21;
  pvmat[2][22] = pv2_22;

  return(0);
}


/*poleswitch01: given a double precision input vector
containing the position in radians of a source on a 
spherical coordinate system, another vector containing
the position of the pole of a new coordinate system
on the old coordinate system, calculates the position
of the point in the new coordinate system, and outputs
this in a third double precision vector.  The desired
RA for the old pole in the new coordinates is also required.
NOTE: the coords must be located in positions 1 and 2
of the double precision vectors, following the NRC
convention, NOT positions 0 and 1.*/
int poleswitch01(double *initpos,double *polepos, double *newpos,double oldpolera)
{
  double x,y,z,xp,yp,zp,thetapole,phipole,theta,phi,thetap,phip;
  int badphip;

  phi = initpos[1];
  theta = initpos[2];
  phipole = polepos[1];
  thetapole = polepos[2];


  z = sin(theta);
  x = cos(theta)*cos(phi-phipole);
  y = cos(theta)*sin(phi-phipole);


  zp = z*sin(thetapole) + x*cos(thetapole);
  xp = x*sin(thetapole) - z*cos(thetapole);
  yp = y;

  if(zp>1.0)
    {
      printf("WEIRD ERROR: POLESWITCH HAS z prime > 1.0!\n");
      printf("THIS VIOLATES BASIC TRIGONOMETRY\n");
      printf("zp-1.0 = %e\n",zp-(double)1.0);
      printf("theta = %f, phi = %f, thetapole=%f, phipole=%f\n",theta,phi,thetapole,phipole);
      printf("xyz, xp yp zp = %f %f %f %f %f %f\n",x,y,z,xp,yp,zp);
      printf("SETTING zp to exactly 1.0.\n");
      zp = (double)1.0;
    }
  thetap = asin(zp);

  phip=0.0;
  if(y==0.0)
    {
      if(x>=0.0)
	{
	  phip = 0.0;
	}
      else if(x<0.0)
	{
	  phip = PI;
	}
    }
  else if(y>0.0)
    {
      phip = PI/(double)2.0 - atan(xp/yp);
    }
  else if(y<0.0)
    {
      phip = (double)3.0*PI/(double)2.0 - atan(xp/yp);
    }

  fflush(stdout);

  phip+=(oldpolera-PI);

  badphip = 0;
  if(phip<0.0||phip>=2.0*PI)
    {
      badphip = 1;
    }
  while(badphip==1)
    {
      if(phip<0.0) phip+=(double)2.0*PI;
      else if(phip>=2.0*PI) phip-=(double)2.0*PI;
      badphip = 0;
      if(phip<0.0||phip>=2.0*PI)
	{
	  badphip = 1;
	}
    }
  newpos[1] = phip;
  newpos[2] = thetap;
  return(1);
}

/*linecount: September 08, 2015. Count the newline '\n' characters
in a file. Not very clever, but useful.*/
int linecount(char filename[])
{
  FILE *fp1;
  int c,lcount;

  if((fp1 = fopen(filename,"r"))!=NULL)
    {
      c = '0';
      lcount = 0;
      while(c!=EOF)
	{
	  c = getc(fp1);
	  if(c=='\n') lcount++;
	}
      fclose(fp1);
      return(lcount);
    }
  else
    {
      printf("ERROR in function linecount: FILE %s DOES NOT EXIST\n",filename);
      return(0);
    }
}


/*readline03: September 05, 2020:
Read the next line from a file that has
already been opened by the calling function. Will work
only of the lines are terminated with a newline '\n'*/
void readline03(FILE *fp1,int maxline,char *linehold)
{
  int i,c;
  i=0;
  c='0';
  while(i<maxline && c!='\n' && c!=EOF)
    {
      c=getc(fp1);
      linehold[i]=c;
      i++;
    }
  if(i>=maxline && c!='\n')
    {
      /*This line was too long. Write a newline at end of what we read.*/
      linehold[maxline-1]= '\n';
      /*Keep going to find the end of the line, but don't load
        anymore characters to the vector.*/
      while(c!='\n' && c!=EOF) c=getc(fp1);
    }
}


/*March 14, 2013: modified to use double precision
input and output variables.
*/
int distradec02(double *ra,double *dec,double *dist,double *pa)
{
  double ra1,ra2,dec1,dec2;
  double initpos[3],polepos[3],newpos[3],oldpolera;

  /*Get RA and DEC in radians*/
  ra1 = ra[1]*PI/180.0;
  ra2 = ra[2]*PI/180.0;
  dec1 = dec[1]*PI/180.0;
  dec2 = dec[2]*PI/180.0;

  /*poleswitch01: given a double precision input vector
    containing the position in radians of a source on a 
    spherical coordinate system, another vector containing
    the position of the pole of a new coordinate system
    on the old coordinate system, calculates the position
    of the point in the new coordinate system, and outputs
    this in a third double precision vector.  The desired
    RA for the old pole in the new coordinates is also required.
    NOTE: the coords must be located in positions 1 and 2
    of the double precision vectors, following the NRC
    convention, NOT positions 0 and 1.*/
  initpos[1] = ra2;
  initpos[2] = dec2;
  polepos[1] = ra1;
  polepos[2] = dec1;
  oldpolera = 0.0;
  poleswitch01(initpos,polepos,newpos,oldpolera);

  /*Calculate distance and position angle*/
  *dist = (PI/2.0 - newpos[2])*ASECPRAD;
  *pa = (2.0*PI-newpos[1])*180.0/PI;


  return(1);
}

