#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

double delta_bond = 0.0;
double delta_angle = 0.0;
double delta_dih = 0.0;
int numbin_bond = 0;
int numbin_angle = 0;
int numbin_dih = 0;
const double MY_PI = 3.14159265358979323846;
const double factor = 57.295779;  //180.0/MY_PI;
int* bondhist;
int* anglehist;
int* dihhist;

/* ********************************************************************************************/
void c_setup_hist_bond(double delta, int maxbin, double* bdist) {
    
    int ibin;
    
    delta_bond = delta;
    numbin_bond = maxbin;
    
    //Initialize arrays
    bondhist     = malloc(sizeof(int)*maxbin);
    
    for (ibin=0; ibin<maxbin ; ibin++) {
       bdist[ibin] = (double) ibin * delta_bond;    
       bondhist[ibin] = 0;
    }
    
}

/* ********************************************************************************************/
void c_setup_hist_angle(double delta, int maxbin, double* adist) {
    
    int ibin;
    
    delta_angle = delta;
    numbin_angle = maxbin;
    
    //Initialize arrays
    anglehist     = malloc(sizeof(double)*numbin_angle);
    
    for (ibin=0; ibin<maxbin ; ibin++) {
       adist[ibin] = (double) ibin * delta_angle ;
       anglehist[ibin] = 0;
    }
    
}

/* ********************************************************************************************/
void c_setup_hist_dih(double delta, int maxbin, double* ddist) {
    
    int ibin;
    
    delta_dih = delta;
    numbin_dih = maxbin;
    
    //Initialize arrays for histogram
    dihhist     = malloc(sizeof(double)*numbin_dih);
    
    for (ibin=0; ibin<maxbin ; ibin++) {
       ddist[ibin] = (double) ibin * delta_dih ;
       dihhist[ibin] = 0;
    }
    
}

/* ********************************************************************************************/

/* ********************************************************************************************/
double distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    
    double delx, dely, delz;
    double rsq, r;
    
    delx = x1 - x2;
    dely = y1 - y2;
    delz = z1 - z2;
    
    rsq = delx*delx+dely*dely+delz*delz;
    r = sqrt(rsq);
    
    return r;
    
}

/* ********************************************************************************************/
double angle(double x1, double y1, double z1,\
             double x2, double y2, double z2,\
             double x3, double y3, double z3) {
    
    double delx1, dely1, delz1;
    double delx2, dely2, delz2;
    double rsq1, rsq2, r1, r2;
    double c, a;

    //1st bond
    delx1 = x1 - x2;
    dely1 = y1 - y2;
    delz1 = z1 - z2;
    
    rsq1 = delx1*delx1+dely1*dely1+delz1*delz1;
    r1 = sqrt(rsq1);

    //2nd bond
    delx2 = x3 - x2;
    dely2 = y3 - y2;
    delz2 = z3 - z2;

    rsq2 = delx2*delx2+dely2*dely2+delz2*delz2;
    r2 = sqrt(rsq2);
    
    //angle
    c = delx1*delx2+dely1*dely2+delz1*delz2;
    c /= r1*r2;
    
    if (c>1.0) c = 1.0;
    if (c<-1.0) c = -1.0;
    
    a = acos(c)*180.0/MY_PI;
    
    return a;
}

/* ********************************************************************************************/
double dihedral(double x1, double y1, double z1,\
                double x2, double y2, double z2,\
                double x3, double y3, double z3,\
                double x4, double y4, double z4)
{
    //Dihedral angle between -180 to 180 with 0 being the cis conformation
    
    // Variables
    double delx1, delx2, delx3, delx2m;
    double dely1, dely2, dely3, dely2m;
    double delz1, delz2, delz3, delz2m;
    double ad = 0.0;
    double ax, ay, az;
    double bx, by, bz;
    double rasq, rbsq,rgsq, rg;
    double rginv, ra2inv, rb2inv, rabinv;
    double c,s;

    //printf("# =======================================\n");
    //printf("cJ1: %.3f %.3f %.3f \n",x1, y1, z1);
    //printf("cJ2: %.3f %.3f %.3f \n",x2, y2, z2);
    //printf("cJ3: %.3f %.3f %.3f \n",x3, y3, z3);
    //printf("cJ4: %.3f %.3f %.3f \n",x4, y4, z4);
    //printf("# =======================================\n");

    //1st bond
    delx1 = x1 - x2;
    dely1 = y1 - y2;
    delz1 = z1 - z2;
    
    //2nd bond
    delx2 = x3 - x2;
    dely2 = y3 - y2;
    delz2 = z3 - z2;
    
    delx2m = -delx2;
    dely2m = -dely2;
    delz2m = -delz2;

    //3rd bond
    delx3 = x4 - x3;
    dely3 = y4 - y3;
    delz3 = z4 - z3;
    
    //c,s calculation
    ax = dely1*delz2m - delz1*dely2m;
    ay = delz1*delx2m - delx1*delz2m;
    az = delx1*dely2m - dely1*delx2m;
    bx = dely3*delz2m - delz3*dely2m;
    by = delz3*delx2m - delx3*delz2m;
    bz = delx3*dely2m - dely3*delx2m;
    
    rasq = ax*ax + ay*ay + az*az;
    rbsq = bx*bx + by*by + bz*bz;    
    rgsq = delx2m*delx2m + dely2m*dely2m + delz2m*delz2m;
    rg = sqrt(rgsq);
 
    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0/rg;
    if (rasq > 0) ra2inv = 1.0/rasq;
    if (rbsq > 0) rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);

    c = (ax*bx + ay*by + az*bz)*rabinv;
    s = rg*rabinv*(ax*delx3+ay*dely3+az*delz3);
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    
    ad = factor * atan2(s,c);
    
    return ad;   
}

/* ********************************************************************************************/
void cross_product(double* x1, double* x2, double* x3)
{

    x1[0] =  x2[1]*x3[2] - x3[1]*x2[2];
    x1[1] = -x2[0]*x3[2] + x3[0]*x2[2];
    x1[2] =  x2[0]*x3[1] - x3[0]*x2[1];

}

/* ********************************************************************************************/
double dihedral2(double x1, double y1, double z1,\
                 double x2, double y2, double z2,\
                 double x3, double y3, double z3,\
                 double x4, double y4, double z4)
{
    //Dihedral angle between -180 to 180 with 0 being the cis conformation

    // Variables
    double del1[3], del2[3], del3[3];
    double n1[3], n2[3], n3[3], n1x2[3];
    double del2_l = 0.0;
    double ad = 0.0;
    double xx,yy;

    //1st bond
    del1[0] = x2 - x1;
    del1[1] = y2 - y1;
    del1[2] = z2 - z1;

    //2nd bond
    del2[0] = x3 - x2;
    del2[1] = y3 - y2;
    del2[2] = z3 - z2;

    //3rd bond
    del3[0] = x4 - x3;
    del3[1] = y4 - y3;
    del3[2] = z4 - z3;

    //c,s calculation
    cross_product(n1, del1, del2);
    cross_product(n2, del2, del3);
    del2_l = sqrt(del2[0]*del2[0]+del2[1]*del2[1]+del2[2]*del2[2]);
    n3[0] = del2[0]/del2_l;
    n3[1] = del2[1]/del2_l;
    n3[2] = del2[2]/del2_l;

    cross_product(n1x2, n1,n2);
    xx = n1x2[0]*n3[0]+n1x2[1]*n3[1]+n1x2[2]*n3[2];
    yy = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];

    ad  = factor * atan2(xx,yy);

    return ad;
}

/* ********************************************************************************************/
double dihedral3(double x1, double y1, double z1,\
                 double x2, double y2, double z2,\
                 double x3, double y3, double z3,\
                 double x4, double y4, double z4)
{
    //Dihedral angle between -180 to 180 with 0 being the cis conformation

    // Variables
    double del1[3], del2[3], del3[3];
    double n1[3], n2[3], n3[3], n1x2[3];
    double del2_l = 0.0;
    double ad = 0.0;
    double xx,yy;

    //1st bond
    del1[0] = x2 - x1;
    del1[1] = y2 - y1;
    del1[2] = z2 - z1;

    //2nd bond
    del2[0] = x3 - x2;
    del2[1] = y3 - y2;
    del2[2] = z3 - z2;

    //3rd bond
    del3[0] = x4 - x3;
    del3[1] = y4 - y3;
    del3[2] = z4 - z3;

    //c,s calculation
    cross_product(n1, del1, del2);
    cross_product(n2, del2, del3);
    del2_l = sqrt(del2[0]*del2[0]+del2[1]*del2[1]+del2[2]*del2[2]);
    n3[0] = del2[0]/del2_l;
    n3[1] = del2[1]/del2_l;
    n3[2] = del2[2]/del2_l;

    cross_product(n1x2, n1,n2);
    xx = n1x2[0]*n3[0]+n1x2[1]*n3[1]+n1x2[2]*n3[2];
    yy = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];

    ad  = factor * atan2(xx,yy);

    return ad;
}

/* ********************************************************************************************/


/* ********************************************************************************************/
void substract_vectors(double* a, double* b, double* c)
{
    
    a[0] = b[0] - c[0];
    a[1] = b[1] - c[1];
    a[2] = b[2] - c[2];
}

/* ********************************************************************************************/
void add_vectors(double* a, double* b, double* c)
{
    
    a[0] = b[0] + c[0];
    a[1] = b[1] + c[1];
    a[2] = b[2] + c[2];
}



/* ********************************************************************************************/
int c_bondDist(int natoms, int nbonds,int dim2, int* bl, double* x, double* y, double* z, int* bondhist)
{
   
  int ibond;
  int atom1 = -1;
  int atom2 = -1;
  int index = 0;
  int ibin = 0;
  double d = 0.0;
   
  for (ibond=0; ibond < nbonds; ibond++) {
    atom1 = bl[index];
    atom2 = bl[index+1];
    index += 2;
    d = distance(x[atom1-1],y[atom1-1],z[atom1-1], x[atom2-1], y[atom2-1], z[atom2-1]);
    ibin = d/delta_bond;
    
    //printf("%6.3f %6.3f %d\n",d,delta_bond,ibin);
    
    if (ibin>=numbin_bond) {
        printf ("\n\n\n ==============================PROBLEM==========================\n ");
        printf ("      Problem with BOND DISTRIBUTION. Increase maxbinBond\n");
        printf ("      ibin: %d numbin_bond: %d distance: %6.3f\n",ibin,numbin_bond,d);
        printf ("      atom1: %d atom2: %d \n",atom1,atom2);
        printf (" ==============================PROBLEM==========================\n\n\n ");
        //return 0;
        continue;
    }
    
    bondhist[ibin] += 1;
    
  }
 
    return 1;
   
}

/* ********************************************************************************************/
int c_angleDist(int natoms, int nangles,int dim2, int* al, double* x, double* y, double* z, int* anglehist)
{
    
  int iangle;
  int ibin;
  int atom1;
  int atom2;
  int atom3;
  int index;
  double a;
 
  index = 0;
 
  for (iangle=0; iangle < nangles; iangle++) {
    
    atom1 = al[index+0];
    atom2 = al[index+1];
    atom3 = al[index+2];    
    index += 3;
    
    a = angle(x[atom1-1], y[atom1-1], z[atom1-1],\
              x[atom2-1], y[atom2-1], z[atom2-1],\
              x[atom3-1], y[atom3-1], z[atom3-1]);
    ibin = a/delta_angle;
    
    //printf("%6.3f %6.3f %d\n",a,delta_angle,ibin);
    //printf("%d %d %d\n",atom1,atom2,atom3);
    
    if (ibin>=numbin_angle) {
        printf ("\n\n\n ==============================PROBLEM==========================\n ");
        printf ("      Problem with ANGLE DISTRIBUTION. Increase maxbinBond\n");
        printf ("      ibin: %d numbin_angle: %d angle: %6.3f\n",ibin,numbin_angle,a);
        printf ("      atom1: %d atom2: %d atom3: %d \n",atom1,atom2, atom3);
        printf (" ==============================PROBLEM==========================\n\n\n ");
        //return 0;
        continue;
    }
    
    anglehist[ibin] += 1;
    
  }
   
    return 1;

}

/* ********************************************************************************************/
int c_angleDist_omp(int natoms, int nangles,int dim2, int* al, double* x, double* y, double* z, int* anglehist)
{
    
  //export OMP_NUM_THREADS=2 --> The best results from 5.7 to 4.2 s
    
  int iangle;
  int ibin;
  int atom1;
  int atom2;
  int atom3;
  int index;
  int tid = 0;
  double a;
  int i,j;
  int local_anglehist[4][numbin_angle];
  
  int angleArray[nangles][3];
 
  for (i = 0; i<4; i++) {
        for (j=0; j<numbin_angle; j++) {
          local_anglehist[i][j] = 0;    
        }
  }
 
  //Atoms from angles
  index = 0;
  for (iangle=0; iangle < nangles; iangle++) {
    
    atom1 = al[index+0];
    atom2 = al[index+1];
    atom3 = al[index+2];    
    index += 3;
    
    angleArray[iangle][0] = atom1;
    angleArray[iangle][1] = atom2;
    angleArray[iangle][2] = atom3;
  }
  
  #pragma omp parallel
  {
    
    tid = omp_get_thread_num();
  
    #pragma omp for nowait private(atom1,atom2,atom3,a,ibin, iangle,tid)  
    for (iangle=0; iangle < nangles; iangle++) {
    
      atom1 = angleArray[iangle][0];
      atom2 = angleArray[iangle][1];
      atom3 = angleArray[iangle][2];
    
      a = angle(x[atom1-1], y[atom1-1], z[atom1-1],\
                x[atom2-1], y[atom2-1], z[atom2-1],\
                x[atom3-1], y[atom3-1], z[atom3-1]);
      int ibin = a/delta_angle;
    
      local_anglehist[tid][ibin]++;

    }
    
    
  }
  
   for (i=0; i<4 ;i++) {
     for (ibin=0; ibin < numbin_angle; ibin++) {
        anglehist[ibin] += local_anglehist[i][ibin];
     }
    }
   
    return 1;

}

/* ********************************************************************************************/
int c_dihDist(int natoms, int ndih,int dim3, int* dl, double* x, double* y, double* z, int* dihhist)
{
    
   int atom1, atom2, atom3, atom4;
   int index;
   int idih;
   int ibin;
   double ad = 0.0;

   index = 0;
 
   for (idih=0; idih < ndih; idih++) {
    
     atom1 = dl[index+0];
     atom2 = dl[index+1];
     atom3 = dl[index+2];
     atom4 = dl[index+3];
     index += 4;

     ad = dihedral(x[atom1-1], y[atom1-1], z[atom1-1],\
                   x[atom2-1], y[atom2-1], z[atom2-1],\
                   x[atom3-1], y[atom3-1], z[atom3-1],\
                   x[atom4-1], y[atom4-1], z[atom4-1]);
     //if ( atom1 == 50 && atom4 == 501) {
     //   printf("%d %d %d %d %f\n", atom1, atom2, atom3, atom4, ad);
     //   exit(0);
     //}
     //Dihedral from 0 to 360, with 180 representing the trans conformation
    
     //if ( atom1 == 13 && atom2 == 12 && atom3 == 11 && atom4 == 6) {
     //   printf("%d %d %d %d\n",atom1,atom2,atom3,atom4);
     //   printf("%6.3f\n",ad);
     //}

    //printf("%6.3f %6.3f %d\n",ad,delta_dih,ibin);
    //printf("%d %d %d %d\n",atom1,atom2,atom3,atom4);
    //printf("******\n");
     
     if (ad < 0) { ad += 360.0; }     
     ibin = ad/delta_dih;

     // Torsional time autocorrelation
     //ad_d = ad*MY_PI/180.0;
     //cosa = cos(ad_d);

    //printf("%6.3f %6.3f %d\n",ad,delta_dih,ibin);
    //printf("%6.3f %6.3f %6.3f\n",ad,ad_d);
    //printf("%d %d %d %d\n",atom1,atom2,atom3,atom4);
    //printf("******\n");

     if ((ibin>=numbin_dih) || (ibin < 0)) {
        printf ("\n\n\n ==============================PROBLEM==========================\n ");
        printf ("      Problem with TORSION DISTRIBUTION. Increase maxbinDih\n");
        printf ("      ibin: %d numbin_dih: %d dih: %6.3f\n",ibin,numbin_dih,ad);
        printf ("      atom1: %d atom2: %d atom3: %d atom4: %d\n",atom1,atom2, atom3, atom4);
        printf (" ==============================PROBLEM==========================\n\n\n ");
        //return 0;
        continue;
     }    

     dihhist[ibin] += 1;    

   }
   
   return 1;


}



/* ********************************************************************************************/
int c_dihDistFlory(int natoms, int ndih,int dim3, int* dl, double* x, double* y, double* z, int* dihhist)
{

   int i, ip1, ip2, im1, im2, j;
   int index;
   int idih;
   double phi1  = 0.0;
   double phi1p = 0.0;
   double phi1p_b = 0.0;
   double phi1_n  = 0.0;
   double phi1p_n = 0.0;
   double psi1  = 0.0;
   double psi1p = 0.0;
   char* filename ="alldihedrals.dat";
   char* filename2 ="alldiads.dat";
   FILE *fp;
   FILE *fp2;

   fp=fopen(filename, "a");
   fp2=fopen(filename2, "a");

   index = 0;

   for (idih=0; idih < ndih; idih++) {

     i   = dl[index+0];
     ip1 = dl[index+1];
     ip2 = dl[index+2];
     im1 = dl[index+3];
     im2 = dl[index+4];
     j   = dl[index+5];
     index += 6;


     phi1  = dihedral2(x[ip1-1], y[ip1-1], z[ip1-1],\
                      x[i-1]  , y[i-1]  , z[i-1]  ,\
                      x[im1-1], y[im1-1], z[im1-1],\
                      x[im2-1], y[im2-1], z[im2-1]);

     phi1p = dihedral2(x[im1-1], y[im1-1], z[im1-1],\
                      x[i-1]  , y[i-1]  , z[i-1]  ,\
                      x[ip1-1], y[ip1-1], z[ip1-1],\
                      x[ip2-1], y[ip2-1], z[ip2-1]);

     phi1p_b = dihedral2(x[ip2-1], y[ip2-1], z[ip2-1],\
                      x[ip1-1]  , y[ip1-1]  , z[ip1-1]  ,\
                      x[i-1], y[i-1], z[i-1],\
                      x[im1-1], y[im1-1], z[im1-1]);

     psi1 =  dihedral2(x[i-1]    , y[i-1]    , z[i-1]  ,\
                       x[ip1-1]  , y[ip1-1]  , z[ip1-1],\
                       x[im1-1]  , y[im1-1]  , z[im1-1],\
                       x[j-1]    , y[j-1]    , z[j-1]);

     psi1p=  dihedral2(x[i-1]    , y[i-1]    , z[i-1]  ,\
                       x[im1-1]  , y[im1-1]  , z[im1-1],\
                       x[ip1-1]  , y[ip1-1]  , z[ip1-1],\
                       x[j-1]    , y[j-1]    , z[j-1]);

     //printf("=========================================\n");
     //printf("phi1: %d %d %d %d %.3f\n",ip1-1, i-1, im1-1, im2-1, phi1);
     //printf("phi1p: %d %d %d %d %.3f\n",im1-1, i-1, ip1-1, ip2-1, phi1p);
     //printf("phi1p: %d %d %d %d %.3f\n",ip2-1, ip1-1, i-1, im1-1, phi1p_b);
     //printf("psi1: %d %d %d %d %.3f\n", i-1, ip1-1, im1-1, j-1, psi1);
     //printf("psi1p: %d %d %d %d %.3f\n",i-1, im1-1, ip1-1, j-1, psi1p);
     //printf("=========================================\n");
     //exit(0);


      if (psi1 > 0) {
         if (phi1>=0) {
           phi1_n = ( 180.0 - phi1);
         }
         else {
           phi1_n = -(phi1 + 180.0);
         }
      }
      else
      {
         if (phi1>=0) {
           phi1_n =  (phi1 - 180.0);
         }
         else {
           phi1_n =  (phi1 + 180.0);
         }
      }


      if (psi1p > 0) {
         if (phi1p>=0) {
           phi1p_n = (180.0 - phi1p);
         }
         else
         {
           phi1p_n = -(phi1p + 180.0);
         }
      }
      else
      {
         if (phi1p>=0) {
           phi1p_n = (phi1p - 180.0);
         }
         else
         {
           phi1p_n = (phi1p + 180.0);
         }
      }

     dihhist[1] += 1;
     dihhist[1] += 1;

     fprintf(fp, "%6.3f\n", phi1_n);
     fprintf(fp, "%6.3f\n", phi1p_n);

     // Write Diadas phi1, phi1p
     if ((phi1_n > -180.0) && (phi1_n < -60.0)) {
         fprintf(fp2, "%s", "u");   //gauche-
     }
     else if ((phi1_n > 60) && (phi1_n < 180.0)) {
         fprintf(fp2, "%s", "g");    //gauche+
     }
     else {
         fprintf(fp2, "%s", "t");    //trans
     }

     if ((phi1p_n > -180.0) && (phi1p_n < -60.0)) {
         fprintf(fp2, "%s\n", "u");   //gauche-
     }
     else if ((phi1p_n > 60) && (phi1p_n < 180.0)) {
         fprintf(fp2, "%s\n", "g");    //gauche+
     }
     else {
         fprintf(fp2, "%s\n", "t");    //trans
     }


   }

   fclose(fp);
   fclose(fp2);

   return 1;


}


/* ********************************************************************************************/
int c_tacticity(int natoms, int ndih,int dim3, int* dl, double* x, double* y, double* z, int* dihhist) {

    //                           (at4)
    //                             |
    //   (at7)---(at5)---(at2)---(at1) --- (at3)---(at6)


    int index;
    int iat1, iat2, iat3, iat4, iat5, iat6, iat7;
    //int idih, j;
    int idih;
    double imp = 0.0;
    double d1 = 0.0;
    double d2 = 0.0;
    double d3 = 0.0;
    char* filename ="tacticity_D.dat";
    char* filename2 ="tacticity_R.dat";
    FILE *fpD;
    FILE *fpR;

    fpD=fopen(filename, "a");
    fpR=fopen(filename2, "a");

    index = 0;
    //j=1;

    for (idih=0; idih < ndih; idih++) {

        iat1 = dl[index+0];
        iat2 = dl[index+1];
        iat3 = dl[index+2];
        iat4 = dl[index+3];
        iat5 = dl[index+4];
        iat6 = dl[index+5];
        iat7 = dl[index+6];
        index += 7;

        imp = dihedral(x[iat1-1],y[iat1-1], z[iat1-1],\
                       x[iat2-1],y[iat2-1], z[iat2-1],\
                       x[iat3-1],y[iat3-1], z[iat3-1],\
                       x[iat4-1],y[iat4-1], z[iat4-1]);

        d1 = dihedral(x[iat3-1],y[iat3-1], z[iat3-1],\
                      x[iat1-1],y[iat1-1], z[iat1-1],\
                      x[iat2-1],y[iat2-1], z[iat2-1],\
                      x[iat5-1],y[iat5-1], z[iat5-1]);

        d2 = dihedral(x[iat1-1],y[iat1-1], z[iat1-1],\
                      x[iat2-1],y[iat2-1], z[iat2-1],\
                      x[iat5-1],y[iat5-1], z[iat5-1],\
                      x[iat7-1],y[iat7-1], z[iat7-1]);

        d3 = dihedral(x[iat6-1],y[iat6-1], z[iat6-1],\
                      x[iat3-1],y[iat3-1], z[iat3-1],\
                      x[iat1-1],y[iat1-1], z[iat1-1],\
                      x[iat2-1],y[iat2-1], z[iat2-1]);

        if (d1 < 0) { d1 += 360.0; }
        if (d2 < 0) { d2 += 360.0; }
        if (d3 < 0) { d3 += 360.0; }
        if (imp < 0) {

            fprintf(fpR, "%.3f %.3f %.3f\n",d1, d2, d3 );
        } else {
            fprintf(fpD, "%.3f %.3f %.3f\n",d1, d2, d3 );
        }

        dihhist[1] += 1;
        dihhist[1] += 1;

    }

    fclose(fpD);
    fclose(fpR);
    return 1;


}