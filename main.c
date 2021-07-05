//=================================================================C
//  2-Dimensional Cavity Flow by Stream function and Vorticity     C
//       with FDM and SOR                                          C
//-----------------------------------------------------------------C
//    Written by Yasunori Ushiro,   2007/06/04                     C
//        ( Tokyo Polytechnic University )                         C
//    後 保範（東京工芸大学）                                      C
//=================================================================C
#include <stdio.h>
#include <math.h>
#include <time.h>
// Global Define
#define NDX  301
#define NDY  301
double PHI[NDX][NDY], OMG[NDX][NDY], OMGW[NDX][NDY] ;
double U[NDX][NDY], V[NDX][NDY] ;
FILE   *FT1 ;
// Function Definition
void PHISOR(int NX, int NY, double H2, double ALP, double EPS) ;
void OUTUVP(int NX, int NY) ;
// Main Program
void main()
{ int     i, j, k, NX, NY, NT, MT, ID, T1, T2, NI1, NI2 ;
  double  Re, DT, T, H, DH, H2, D2, EPS, ALP, CPU ;
  char    FLOW1[13]={'F','L','O','W','1','-','?','?','.','t','x','t',} ;
  char    NTABL[10]={'0','1','2','3','4','5','6','7','8','9'} ;
//  Initial Data
  ID = 0 ;
  printf("Type In  NX,Re \n") ;
  scanf("%d %lf",&NX,&Re) ; 
  printf("DT,NT(Total NO.),MT(Output Interval) \n") ;
  scanf("%lf %d %d",&DT,&NT,&MT) ;
  if(NX >= NDX) { NX = NDX - 1 ; }
  NY = NX ;
  printf("# Cavity Flow Numerical Analysis \n") ;
  printf("NX=%d NY=%d Re=%f DT=%f  NT,MT=%d %d \n",NX,NY,Re,DT,NT,MT) ;
//  Initial Constant
  T  = 0.0 ;
  DH = NX ;
  H  = 1.0/DH ;
  H2 = H*H ;
  D2 = DH*DH ;
//  Set Initial Condition
  for (i=0; i<=NX; i++) {
    for (j=0; j<=NY; j++) {
      OMG[i][j] = 0.0 ;
      PHI[i][j] = 0.0 ;  }
    }
//  SOR Parameter (ALP) 
  ALP = 1.0 + log(NX*1.0)/log(NDX*1.2) ;
  EPS = 1.0e-4 ;
//  Main Loop
  T1 = clock() ;
  for (k=1; k<=NT; k++) {
    T = T + DT ;
//   Set Boundary Condition
    for (i=0; i<=NX; i++) {
      OMG[i][0]  = -2.0*PHI[i][1]*D2 ; 
      OMG[i][NY] = -2.0*(PHI[i][NY-1] + H)*D2 ;  }
    for (j=1; j<NY; j++) {
      OMG[0][j]  = -2.0*PHI[1][j]*D2 ;
      OMG[NX][j] = -2.0*PHI[NX-1][j]*D2 ;  }
//   Compute Omega
    for (i=1; i<NX; i++) {
      for (j=1; j<NY; j++) {
        OMGW[i][j] = OMG[i][j] + DT*D2*( (
                  - (PHI[i][j+1]-PHI[i][j-1])*(OMG[i+1][j]-OMG[i-1][j]) 
                  + (PHI[i+1][j]-PHI[i-1][j])*(OMG[i][j+1]-OMG[i][j-1])
                  ) / 4.0 + (OMG[i][j-1]+OMG[i-1][j]-4.0*OMG[i][j]
                  + OMG[i+1][j]+OMG[i][j+1]) / Re ) ;  }
       } 
//   Copy back Omega
    for (i=1; i<NX; i++) {
      for (j=1; j<NY; j++) {
        OMG[i][j] = OMGW[i][j] ; }
     }
//   Compute Phi
    PHISOR(NX, NY, H2, ALP, EPS) ;
//   Output U,V,Phi
    if( k%MT == 0) {
      NI1 = ID/10 ;
      NI2 = ID - NI1*10 ;
      FLOW1[6] = NTABL[NI1] ;
      FLOW1[7] = NTABL[NI2] ;
      FT1 = fopen(FLOW1,"w") ;
      fprintf(FT1,"# Cavity Flow Numerical Analysis \n") ;
      fprintf(FT1,"NX=%d NY=%d Re=%f DT=%f NT,MT=%d %d \n",NX,NY,Re,DT,NT,MT) ; 
      fprintf(FT1,"#Time=%f,  Step=%d \n",T,k) ; 
      printf("#Time=%f,  Step=%d \n",T,k) ; 
      OUTUVP(NX, NY) ;
      fclose(FT1) ;
      if(ID < 99) { ID = ID + 1 ; }
     }
   }
  T2 = clock() ;
  CPU  = (double)(T2 - T1)/CLOCKS_PER_SEC ;
  printf(" NX=%4d,  NT=%5d,  Time(s)=%9.2f \n",NX,NT,CPU) ; 
  } 
//=================================================================C
void PHISOR(int NX, int NY, double H2, double ALP, double EPS)
//=================================================================C
//  Solve Ax=b by SOR with 2 dimensional FDM                       C
//    Given Omega ( Acceleration factor )                          C
//-----------------------------------------------------------------C
//    PHI,OMG ;  Global Array                                      C
//    NX      I*4, In,  Grid Numbers on X-axis                     C
//    NY      I*4, In,  Grid Numbers on Y-axis                     C
//    H2      R*8, In,  H2=H**2                                    C
//    ALP     R*8, In,  SOR Acceleration factor                    C
//    EPS     R*8, In,  if ||r||/||b|| <= EPS --> return           C
//-----------------------------------------------------------------C
//    Written by Yasunori Ushiro,   2007/06/04                     C
//        ( Tokyo Polytechnic University )                         C
//=================================================================C
{ int     i, j, k ;
  double  BN, RN, ERR, W1, R ;
//  Get 2-Norm B=D2*OMG
  BN = 0.0 ;
  for (i=1; i<NX; i++) {
    for (j=1; j<NY; j++) {
      W1 = H2*OMG[i][j] ;
      BN = BN + W1*W1 ;  }
    } 
//   Main Loop
  for (k=1; k<= NX*NY; k++) {
    RN = 0.0 ;
    for (i=1; i<NX; i++) {
      for (j=1; j<NY; j++) {
        R = (H2*OMG[i][j] + PHI[i][j-1] + PHI[i-1][j] + PHI[i+1][j]
            + PHI[i][j+1] )/4.0 - PHI[i][j] ;
        PHI[i][j] = PHI[i][j] + ALP*R ;
        W1 = R*4.0 ;  
        RN = RN + W1*W1 ;  }
      }
//   if(ERR <= EPS) return
    ERR = sqrt(RN/BN) ;
    if(ERR <= EPS) break ;
   }
 }
//=================================================================C
void OUTUVP(int NX, int NY)
//=================================================================C
//  Compute U,V and Output U,V,Phi                                 C
//    U=d(PHI)/dy, V=-d(PHI)/dx                                    C
//-----------------------------------------------------------------C
//    PHI,U,V ; Global Array                                       C
//    NX      I*4, In,  Grid Numbers on X-axis                     C
//    NY      I*4, In,  Grid Numbers on Y-axis                     C
//-----------------------------------------------------------------C
//    Written by Yasunori Ushiro,   2007/06/04                     C
//        ( Tokyo Polytechnic University )                         C
//=================================================================C
{ int     i, j ;
  double  XV, YV, X, Y, DX, DY ;
//   Initial 
  XV = NX/2.0 ;
  YV = NY/2.0 ;
//   Compute U,V
  for (i=1; i<NX; i++) {
    for (j=1; j<NY; j++) {
      U[i][j] = (PHI[i][j+1] - PHI[i][j-1])*YV ;
      V[i][j] = (PHI[i-1][j] - PHI[i+1][j])*XV ;  }
    }
//    Boundary on Y=0,1
  for (i=0; i<=NX; i++) {
    U[i][NY] = 1.0 ;
    V[i][NY] = 0.0 ;
    U[i][0]  = 0.0 ;
    V[i][0]  = 0.0 ;  }
//    Boundary on X=0,1
  for (j=1; j<NY; j++) {
    U[0][j]  = 0.0 ;
    V[0][j]  = 0.0 ;
    U[NX][j] = 0.0 ;
    V[NX][j] = 0.0 ;  }
//   Output
  DX = 1.0/NX ;
  DY = 1.0/NY ;
  for (j=0; j<=NY; j++) {
    Y = j*DY ;
    for (i=0; i<=NX; i++) {
      X = i*DX ;
      fprintf(FT1,"%9.5f %9.5f %9.5f %9.5f %9.5f \n",
                  X,Y,PHI[i][j],U[i][j],V[i][j]) ;  }
    }
  }