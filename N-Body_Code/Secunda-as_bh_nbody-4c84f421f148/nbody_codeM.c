#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#define glimax          11
#define glnmax          1200
#define glncol          7
#define MASSJUP         0.0009548
#define PI    3.141592653589793
#define TWOPI           6.283185307179586
#define TORAD   PI/180.0
#define TODEG           180.0/PI
#define G               0.01720209895  // the Gaussian grav. constant
#define G_CGS   6.674e-8
#define SUN2GR    1.989e33
#define AU2CM   1.496e13
#define SDCONV    1.12521e-7  // convert Msun/AU^2 to gram/cm^2
#define SPEEDOFLIGHT    3e10  // cm per second
#define UNO   1.0;
#define TWO             2.0;
#define UNO2    0.5;
#define ERROR   1.0e-13;
#define CERO    0.0;
#define MAXNBOD         200              // maximum number of bodies except of the central one
#define BETA            50.0
#define ASPECT_RATIO  0.05             //not used, calculated based on temperature
#define NPLMAX          10
#define FLT_MAX         3.4028E+38F      // maximum float value

#define NTS     1000000   //number of Time Steps in time.dat from gas profiles (Brandon)
#define NRAD    5049           //radial resolution in gas simulation
#define NStart    3           //number of planets in elements.dat
#define MODES   200          //num of harmonics in turb potential
#define gamma   4.25e-4         //strength of turb potential

char * gasdirect  ="./SGtimefiles/";//location of gas profiles


// MAXNBOD = glnmax/6


typedef double vektor[glnmax];

 long nv,  NBOD, NEW_NBOD, IMIN,IBMIN, JMIN, JBMIN;
 int CLOSE,BINFORM,COLLISION,BODY_CREATE, PRINT_TIME, ESCAPE, GASDRAG, MIGRATION_TYPE1, ADIABATIC, NPL, NPHI, LamIndex[MODES], UPDATE_OPACITY; //NRAD
 
double xbe, tbmove,tbmerge, thard, htrybe, hkovki, max_time, xki,htemp2print, BODY0, BODY[MAXNBOD], BODY_AFTER[MAXNBOD], BODY_BAFTER[MAXNBOD], RHO, Lambda[6*MODES], LamTime[MODES],ProbLog[257],Dyn_Disp[MAXNBOD],Disp_Array[MAXNBOD],I_BVECTOR[6],J_BVECTOR[6],OPACITY[MAXNBOD];
 double RESCMIN, RESCMAX, SDEXP, SIGMA_0, DISP_BEGIN, DISP_DURATION, EPOCH, ADINDEX, EXPSD, DCRITICAL,R_inner,R_outer,D_r;
 
 
 vektor ybe, yki, yki_after, yki_bafter,adj_pos, yki_after_esc;
 
 
 double glx[glimax];
 double gld[glnmax][glncol];
 
FILE *ptr_myfile, *surf, *fin1,*fin2,*fin3, *foutb, *foutm1, *foutm2, *foute, *fout_sd, *fout2, *fout3, *fout4, *fout5, *fin_time, *abkap, *hconv;


/*structure for printing into output.dat*/
struct rec
        {
        float i,j,mass,time,ax,ex,inc,MAn,omP,OmN,x,y,z,vx,vy,vz,oldhtemp;        
        };


struct ev
  {
  float alpha,beta,d_torque;
  };

struct kap
  {
  float a,b,rho,kk;
  };

struct hcon
  {
  float h;
  }; 

struct sigR
  {
  float sig, r;
  };

int int_min(i,j)
int i,j;
{
  if (i < j)
     return i;
  else
     return j; 
}


int int_max(i,j)
int i,j;
{
  if (i>j)
      return i;
  else 
      return j;
}

double my_rand()
{
  return (double) rand()/RAND_MAX;
}

//calculates size of body of given mass, assumes uniform density of
//2g/cm^3
// JMB: this has been replaced for the black hole case.
double radius(i,body)
int i; double *body;
{

    /*   double sugar;
    
    sugar = body[i]*SUN2GR*1.5/(RHO*TWOPI);
    sugar = pow(sugar,0.3333333333);
    sugar = sugar/AU2CM;
   
    return sugar; */  

    double schwarz;
    schwarz = 2.0*body[i]*SUN2GR*G_CGS/(SPEEDOFLIGHT*SPEEDOFLIGHT);
    schwarz = schwarz/AU2CM; // consistent with all AU units. ugh.
    return schwarz;
    /*  JMB 4/14/17 */
}

//calculates the mutual hill radius
double rhill(m1,m2,M0)
double m1,m2,M0;
{
  double a1, a2, sfact, ffact3;
  a1=pow(I_BVECTOR[0]*I_BVECTOR[0]+I_BVECTOR[1]*I_BVECTOR[1]+I_BVECTOR[2]*I_BVECTOR[2],0.5);
  a2=pow(J_BVECTOR[0]*J_BVECTOR[0]+J_BVECTOR[1]*J_BVECTOR[1]+J_BVECTOR[2]*J_BVECTOR[2],0.5);

  ffact3=(m1+m2)/(3*M0);
  sfact=(a1+a2)/2;
  return pow(ffact3,1/3.)*sfact;

}

//want sigma in g/cm^2, masses in msun, a in parsec, abin in AU
double t_hard(temp,sigma,M0,mbin,a,abin)
double temp,sigma,M0,mbin,a,abin;
{
  double factor1,factor2,factor3,factor4,factor5,factor6;
  sigma=sigma/SDCONV; //convert from AU/MSUN^2 to g/cm^2
  a=a*4.84814e-6; //convert from AU to pc
  factor1=pow((temp/2000),3);
  factor2=1./(sigma/100);
  factor3=sqrt(3e6/M0);
  factor4=(2/mbin)*(2/mbin);
  factor5=sqrt(pow((a/0.1),3));
  factor6=abin/50;
  return (900*factor1*factor2*factor3*factor4*factor5*factor6);
}

double Solve_Kep_eq( double e, double M )
{
  double E0, E1, error;

    if( M != 0.0 )  {
    E0 = M + e*(sin(M))/(1-sin(M+e)+sin(M));
    do{
      E1 = E0 - ( E0-e*sin(E0)-M )/( 1-e*cos(E0) );
      error = fabs(E1-E0);
      E0 = E1;
    }
    while( error > 1e-12 );
  }
    else
    E1 = M;

  return E1;
}

void get_e_i(body0,body1,r,x,y,z,vx,vy,vz,e,i)
double body0,body1,r,x,y,z,vx,vy,vz;
double *e, *i;
{
  double mu,h,cx,cy,cz,c,v2,exc2,exc,inc,cinc;
  
  mu=G*G*(body0+body1);
  v2=vx*vx+vy*vy+vz*vz;
  cx=y*vz -z*vy;
  cy=z*vx - x*vz;
  cz=x*vy-y*vx;
  c=sqrt(cx*cx+cy*cy+cz*cz);
  h=0.5*v2-mu/r;
  exc2=1.0+2.0*h*c*c/(mu*mu);
  if(fabs(exc2) <= 1.0e-10) exc2=0.0;
  exc=sqrt(exc2);
  cinc=cz/c;
  inc=acos(cinc);
  *e=exc;
  *i=inc;
}

void element(tomeg0,tomeg1,x,y,z,vx,vy,vz,ax,exc,ink,omegaP,OmegaN,kozan)
double tomeg0,tomeg1,x,y,z,vx,vy,vz;      // coordinates and velocities of the orbiting body
double *ax;     //semi-major axis
double *exc;    // eccentricity
double *ink;  // inclination          
double *omegaP; // argument of perihelion
double *OmegaN; // longitude of node
double *kozan;  // mean anomaly
{
  double mu,r,u,h,cx,cy,cz,c,lapx,lapy,lapz,lap;
  double lax, lexc, link, cink, sOmN, cOmN, lOmegaN, lomegaP, 
     seb2, sqrexc;
  double TST2,TST3;
  double swf, cwf, wf, ecu, esu, uu, ecc1, tnu, ff, w, lkozan, lambda;

  TST2 = 2.0e-5;
  TST3 = 1.0e-10;  
  
  
  mu = G*G*(tomeg0+tomeg1);    // tomeg0 - mass of the central body
                               // tomeg1 - mass of the orbiting body
  r = sqrt(x*x + y*y + z*z);
  seb2 = vx*vx + vy*vy + vz*vz;
  u = x*vx + y*vy + z*vz;


  h = 0.5*seb2 - mu/r; //energia
  
  cx = y*vz - z*vy; //AGMX
  cy = z*vx - x*vz; //AGMY
  cz = x*vy - y*vx; //AGMZ
  c = sqrt(cx*cx + cy*cy + cz*cz);          //absolut value of the angular momentum
  
  lapx = vy*cz - mu*x/r - vz*cy;                      // Laplace-vector
  lapy = vz*cx - mu*y/r - vx*cz;
  lapz = vx*cy - mu*z/r - vy*cx;
  lap = sqrt(lapx*lapx + lapy*lapy + lapz*lapz);
  
  
  
  lax = -mu/(2.0*h);
  
  sqrexc = 1.0+2.0*h*c*c/(mu*mu); 
  if (fabs(sqrexc) <= TST3) sqrexc = 0.0; 
  lexc= sqrt(sqrexc);
  
  
  
  cink = cz/c;
  link = acos(cink);
  
  if (link < TST2) 
      link = 0.0;
   
  if (link != 0.0) {
      if (cz > 0.0) {
          cOmN = - cy;
    sOmN =   cx; 
      } else {
        cOmN =  cy;
        sOmN = -cx; 
  }
          
  
  
  lOmegaN = atan2(sOmN,cOmN);
  
  } else lOmegaN = 0.0;
  
  if (lOmegaN < 0.0) {
                      lOmegaN = 2.0*PI + lOmegaN; 
          } 
  
  swf = (-x*sin(lOmegaN)+y*cos(lOmegaN))/(r*cos(link));
  cwf = ( x*cos(lOmegaN)+y*sin(lOmegaN))/r;
  
  wf  =  atan2(swf,cwf); 
                        // node-body distance
  
       if (wf < 0.0) wf = 2.0*PI + wf;
       
       
  if (lexc != 0.0) {                           // uu eccentric anomaly
                    ecu = 1.0 + 2.0*h*r/mu;
                    esu = sqrt(-2.0*h)*u/mu;
         uu = atan2(esu,ecu);
  }
  
  if (uu < 0.0) uu = 2.0*PI + uu;
  
  ecc1 = sqrt((1.0+lexc)/(1.0-lexc));
  tnu  = ecc1*tan(uu/2.0);
  ff   = 2.0*atan(tnu);                       // ff true anomaly 
  
  if ((ff >= 0.0) && (uu >= PI)) ff = ff + PI;
  if ((ff <  0.0) && (uu <  PI)) ff = ff + PI;
  if  (ff <  0.0)                ff = 2.0*PI + ff;
  if  (ff >= 2.0*PI)             ff = ff - 2.0*PI;
  
  
  
  lkozan = uu - esu;                          // lkozan mean anomaly
  
  if (lkozan < 0.0) {
                     lkozan = 2.0*PI + lkozan;
                     if (lkozan >=2.0*PI) {
                                           lkozan = lkozan - 2.0*PI;
         } else {
           uu =     0.0;
                 ff =     0.0;
                 lkozan = 0.0;
           }
  }
  
  
  
 //    w PERIHELION LONGITUDE

       w = wf - ff + lOmegaN;
       
       if (w < 0.0)    w = 2.0*PI + w;
       if (w >= 2.0*PI) w = w - 2.0*PI;
       
//     lambda MEAN LONGITUDE

       lambda = w + lkozan;
       
       if (lambda >= 2.0*PI) lambda = lambda - 2.0*PI;
       
//     ARGUMENTO OF PERICENTER             

       lomegaP = w-lOmegaN;
       
       if (lomegaP < 0.0) lomegaP = lomegaP + 2.0*PI; 
  
  *ax = lax;
  *exc= lexc;
  *ink= link * TODEG;
  *omegaP= lomegaP * TODEG;
  *OmegaN= lOmegaN * TODEG; 
  *kozan = lkozan * TODEG;
  
}


void coordinate(M0,UM,SMA,ECC,INC,W,OMN,MEAN,x,y,z,vx,vy,vz)
double M0, UM, SMA, ECC, INC, W, OMN, MEAN;
double *x, *y, *z, *vx, *vy, *vz;
{
double MU,EXAN,SMA3,VELM2,VELM,CU,SU,CUE,PX,QX,PY,QY,PZ,QZ,SMA1,RADI,UDOT;

// M0 central mass, UM mass of the other body, W arg of pericenter
// OMN  longitude of the ascending node

//   INITIALIZATION OF VARIABLES

      SMA3  = SMA*SMA*SMA; // semi-major axis
      MU    = G*G*(M0+UM); // UM mass of the smaller body
      VELM2 = MU/SMA3;    
      VELM  = sqrt(VELM2);   // units of 1/time
      //printf("INC Three = %g \n",INC); *JIA*9/25*
      INC = INC*TORAD;   // TORAD = to radians
      //printf("INC Four = %g \n",INC); *JIA*9/25*
      //printf("MEAN One = %g \n",MEAN); *JIA*9/25*
      MEAN = MEAN*TORAD;
      //printf("MEAN Two = %g \n",MEAN); *JIA*9/25*
      //printf("W One = %g \n",W); *JIA*9/25*
      W = W*TORAD;
      //printf("OMN One = %g \n",OMN); *JIA*9/25*
      OMN = OMN*TORAD;
      EXAN = Solve_Kep_eq(ECC,MEAN);
	      
//    DIRECTOR COSINES    
    
      CU  = cos(EXAN);
      printf("CU One = %g \n",CU);
      SU  = sin(EXAN);
      printf("SU One = %g \n",SU);
      printf("EXAN One = %g \n",EXAN);
      CUE = CU - ECC;
      printf("CUE One = %g \n",CUE);
      printf("SMA One = %g \n",SMA);
      PX  = cos(W)*cos(OMN) - sin(W)*sin(OMN)*cos(INC);
      printf("PX One = %g \n",PX);
      //printf("INC Five = %g \n",INC); *JIA*9/25*
      QX  =-sin(W)*cos(OMN) - cos(W)*sin(OMN)*cos(INC);
      //printf("INC Six = %g \n",INC); *JIA*9/25*
      PY  = cos(W)*sin(OMN) + sin(W)*cos(OMN)*cos(INC);
      //printf("INC Seven = %g \n",INC); *JIA*9/25*
      QY  =-sin(W)*sin(OMN) + cos(W)*cos(OMN)*cos(INC);
      //printf("INC Eight = %g \n",INC); *JIA*9/25*
      PZ  = sin(W)*sin(INC);
      printf("PZ One = %g \n",PZ);
      //printf("INC Nine = %g \n",INC); *JIA*9/25*
      QZ  = cos(W)*sin(INC);
      printf("QZ One = %g \n",QZ);
      //printf("INC Ten = %g \n",INC); *JIA*9/25*
      SMA1 = SMA*sqrt(1.0-ECC*ECC);
      printf("SMA1 One = %g \n",SMA1);
      RADI = (1.0-ECC*CU)*SMA;
      UDOT = SMA*VELM/RADI;

//     HELIOCENTRIC COORDINATES AND VELOCITIES X AND XDOT [AU] & [AU/DAY]

      *x = SMA*PX*CUE+SMA1*QX*SU;           // X
      *y = SMA*PY*CUE+SMA1*QY*SU;           // Y
      *z = SMA*PZ*CUE+SMA1*QZ*SU;           // Z
      *vx = UDOT*(-SMA*PX*SU+SMA1*QX*CU);   // velocity components
      *vy = UDOT*(-SMA*PY*SU+SMA1*QY*CU);   
      *vz = UDOT*(-SMA*PZ*SU+SMA1*QZ*CU);   
      printf("*x One = %g \n",*x);
      printf("*z One = %g \n",*z);
}  
/*
//calculates speed of sound/keplerian velocity
// has been adjusted to cgs units
double get_aspect_ratio(r,T,Mcenter)
  double r, T,Mcenter;
{
  double vKep, k, adiabat, mH, cs, asp_ratio;
  k=1.3806503e-16;//Boltzman in cgs
  adiabat=1.67;// gamma for ideal gass
  mH=1.67262e-24; //mass of hydrogen atom in g
  r = r*AU2CM;  // convert r from AU to cm
  cs=sqrt(adiabat*k*T/mH); // cm/s
  //  printf("cs: %g\n",cs);
  //printf("Mcenter: %g\n",Mcenter);
  vKep=sqrt(G_CGS*Mcenter*SUN2GR/r);  // cm/s
  //printf("vKep: %g\n",vKep);
  //printf("r: %g\n",r);
  asp_ratio=cs/vKep;
  //printf("aspect ratio: %g\n",asp_ratio);
  return asp_ratio;
}
*/

double get_sigma(r,radii,surfd)
double r, *radii,*surfd;
{
  int index;
  double d_r1,r_remainder,s_dens;

  if (r>R_outer){
     r=R_outer -0.01;
   }
  
  //  d_r1=1.0/D_r;
  // index=((r-R_inner)*d_r1);
  double minvalue=FLT_MAX;
  double testvalue=FLT_MAX;
  int i;
  for (i=0;i<NRAD;i++){      
      testvalue = fabs(radii[i]-r);
      if (testvalue < minvalue) {
	  minvalue = testvalue;
	  index = i;
	  }
      }

  r_remainder=(r-radii[index])*(r-radii[index])/NRAD;
  
  s_dens = surfd[index]+r_remainder*(surfd[index+1]-surfd[index]);
  return s_dens;
}

//linear interpolation between two nearest datapoints in gas profile to location r
//alpha is local gradient of density,beta is local gradient of temp. tables calculated in Wlads code
void get_local_gas_properties(r,radii,aspR,expal,surfd,temp,beta,s_dens,T,aspR_local,alpha_local,beta_local)
double r, *radii,*expal,*surfd,*temp,*beta,*aspR;
double *s_dens,*T,*alpha_local,*beta_local,*aspR_local;
{
  int index;
  double d_r1,r_remainder;

  //  printf("Rinner %g Router %g\n",R_inner,R_outer);
  if (r>R_outer){
     r=R_outer -0.01;
   }
  // r is in AU

  d_r1=1.0/D_r;
  index=((r-R_inner)*d_r1); // linear interp gives a misleading
			    // (wrong) index
  r_remainder=(r-radii[index])*(r-radii[index])/NRAD;
  // printf("dr1 %g D_r %g\n",d_r1,D_r);
  //printf("r %g radii[index] %g \n",r,radii[index]); 

  /* radii[index]  does not  represent the  correct radius  at all  */

  /*  need to minimize radius value with values in radii array
      loop through the array and find the min value
      JMB 5/31/17
  */
  double minvalue=FLT_MAX;
  double testvalue=FLT_MAX;
  int i;
  for (i=0;i<NRAD;i++){      
      testvalue = fabs(radii[i]-r);
      if (testvalue < minvalue) {
	  minvalue = testvalue;
	  index = i;
	  }
      }

  r_remainder=(r-radii[index])*(r-radii[index])/NRAD;

  *alpha_local = expal[index]+r_remainder*(expal[index+1]-expal[index]);
  
  *s_dens = surfd[index]+r_remainder*(surfd[index+1]-surfd[index]);

  *T = temp[index]+r_remainder*(temp[index+1]-temp[index]);
  
  *beta_local = beta[index]+r_remainder*(beta[index+1]-beta[index]);

  *aspR_local = aspR[index]+r_remainder*(aspR[index+1]-aspR[index]);
}

//derives speed of sound from local temp. This is used to fnd scale height of disk
double get_speed_of_sound(r,radii,aspR,Mcenter)
double r,*radii,*aspR,Mcenter;
{
  
  int index;
  double d_r1,r_remainder,H,k,mH,cs,adiabat;

  if (r>R_outer){
     r=R_outer -0.01;
   }
  
  d_r1=1.0/D_r;
  // index=((r-R_inner)*d_r1);
  double minvalue=FLT_MAX;
  double testvalue=FLT_MAX;
  int i;
  for (i=0;i<NRAD-1;i++){      
      testvalue = fabs(radii[i]-r);
      if (testvalue < minvalue) {
	  minvalue = testvalue;
	  index = i;
	  }
      }
  r_remainder=(r-radii[index])*(r-radii[index])/NRAD;
  //printf("Radii=%g\n",radii[index]);
  H = aspR[index]+r_remainder*(aspR[index+1]-aspR[index]);
  //printf("speed of sound T: %g, r_r: %g index: %d radius: %g \n",T,r_remainder,index,r);
  /*k=1.3806503e-16;//Boltzman in cgs
  adiabat=ADINDEX; // monotomic ideal gass
  mH=1.67262e-24; //mass of hydrogen atom in g
  cs=sqrt(adiabat*k*T/mH);   // cm per second !
  return cs*5.77547828e-5; // cm/s to AU/day 
  */
  //old calculation doesnt take into account radiation pressure changing the scale height
  cs= H*sqrt(G*G*Mcenter/r);
  return cs; // cm/s to AU/day 
}
 
/* major units hairyness here.  JMB 4/13/17  */
void get_migration_damping(body0,bodyi,x,y,z,vx,vy,vz,r,sig,T,aspR,alpha,beta,Torque,Damp,Theta)
double body0,bodyi,x,y,z,vx,vy,vz,r,sig,T,aspR,alpha,beta,*Torque,*Damp,Theta;
{
    double ecc, inc, entropy_exp, h, rte,rti,d_torque, torque_rn2,rv,torque_norm;

  get_e_i(body0,bodyi,r,x,y,z,vx,vy,vz,&ecc,&inc);
  h=aspR;
  //printf("e i h =%lg  %lg   %lg\n",ecc, inc, h);


  Theta=((Theta/SUN2GR)*(pow(AU2CM,2)))*(60*60*24); //converting this from cgs to AU, Msun

  rte=h*(r*body0*sqrt(r*body0))*(h*h*h-0.14*ecc*ecc*h+0.06*ecc*ecc*ecc+0.18*ecc*inc*inc)/(0.78*bodyi*sig*G);
  //printf("inc Eleven = %g \n",inc); *JIA*9/25*
  rti=h*(body0*sqrt(body0))*(h*h*h-0.3*inc*inc*h+0.24*inc*inc*inc+0.14*ecc*ecc*inc)/(0.544*bodyi*sig*G*sqrt(r));
  //printf("inc Twelve = %g \n",inc); *JIA*9/25*
  entropy_exp = beta - alpha*(ADINDEX-1.0);
  d_torque = ((1.0/ADINDEX)*(-0.85 - alpha - 1.7*beta + 7.9*entropy_exp/ADINDEX)*Theta*Theta + (-0.85-alpha-0.9*beta))/((Theta+1)*(Theta+1));
  printf("inc Thirteen = %g \n",inc);// *JIA*9/25*
  if (inc > 3.14)       //*JIA*9/25*
    {
      d_torque=0;
      printf("d_torque = %g \n",d_torque);
    }

//everything file is taking too much memory space
/*  struct ev my_recorde;
  my_recorde.alpha=alpha;
  my_recorde.beta=beta;
  my_recorde.d_torque=d_torque;
  fwrite(&my_recorde, sizeof(struct ev), 1, fout5);*/

  printf("d_torque = %g   Theta %g \n",d_torque,Theta);
  /* prior torque_rn2 statement */
  torque_rn2=d_torque*G*G*bodyi*bodyi*sig/(h*h*body0*r);  // these quantities do
						    // not match the
						    // paper JMB
						    // 4/13/17
  /* JMB 4/18/17 fixed the above statement to include another factor
     of bodyi.  torque_rn2 = torque/r2 ; dimensions were missing a
     factor of mass, which can be recovered by looking at
     normalization equation in paper. Horn et al results need to be examined  */
  rv=x*vx+y*vy+z*vz;  // r dot v
  //printf("x,y: %g %g \n",x,y);
  /* Torque[i] below is actually a Force. x and y components have
     been derived by JMB/MML and confirmed to be correct 4/18/17 */
  Torque[0]= -1.0*torque_rn2*y; 
  Torque[1]= 1.0*torque_rn2*x;
  //printf("t =%lg  x %g  y %g\n",Torque[1],x,y);

  /* dimensions here are ridiculous */
  Damp[0]=-2.0*rv*x/rte;  // L
  Damp[1]=-2.0*rv*y/rte;  // L
  Damp[2]=-z*2.0*rv/rte -vz/rti; // 1st term L, 2nd term L2/T2  WTF
}

 
  /*calculates damping timescale for eccentricity and inclination, following creswell and Nelson 2007. No damping force on semimajor axis. Instead, a dimensionless torque is calculated based on temp, density and entropy gradients, following Paardekooper et al*/
void damping_timescales(r,s_dens, alpha,beta, exc, inc,h, index, taue, taui, dless_t)
double r, s_dens,alpha,beta,exc,inc,h;
int index;
double *taue, *taui, *dless_t; // damping timescales
{
  double entropy_exp, dless_torque,te,ti,twave,e_h,i_h;
  int i, i_int; 
  
  /* find torques following Cresswell and Nelson 2007 */
  
  e_h=exc/h;
  i_h=inc/h;
  //printf("inc Fourteen = %g \n",inc); *JIA*9/25*     
  twave=h*h*h*h/(s_dens*sqrt(r)*G);
  te= twave*(1.0-0.14*e_h*e_h + 0.06*e_h*e_h*e_h + 0.18*e_h*i_h*i_h)/0.78;
  ti=twave*(1.0-0.3*i_h*i_h + 0.24*i_h*i_h*i_h + 0.14*e_h*e_h*i_h)/0.544;

  entropy_exp = beta - alpha*(ADINDEX-1.0);

  switch (index) {
          case 0:
	    //temp_exp = - 1.0;
          dless_torque = - 0.85 - alpha - 0.9*beta;
          break;
  

          case 1:
	  //this is adiabatic unsaturated torque
          dless_torque = (1.0/ADINDEX)*(-0.85 - alpha - 1.7*beta + 7.9*entropy_exp/ADINDEX);
          break;
  }

 
  *dless_t = dless_torque;
  
  /* damping timescales for the eccentricity and inclination */   
  *taue = te;
  *taui = ti;
}


// procedure not used in this version
// confirmed JMB 4/14/17

double Vol_Density(a,h,x,y,z,sig_0)
double a, h, x, y, z, sig_0;
{
   double r;
  
   r=sqrt(x*x + y*y);
   return sig_0*pow(r,a-1.0)/(h*sqrt(2.0*PI)*1.49e13)*exp(z*z/(2.0*h*h*r*r));  
}


// procedure not used in this version
// confirmed JMB 4/14/17

double mean_free_path(rho_gas)
double rho_gas;
{
 double mu_H2, mcross, m_prot, lambda;
 
   mu_H2 = 2.3;       /* mean molecular weight of the H2 molecule*/
  mcross = 2e-15;     /* collisional cross section of the H2 molecule*/
  m_prot = 1.6726e-24;/* proton's mass in gramm */
  
  /* mean free path of the H2 molecule in cm */
  
  lambda = mu_H2*m_prot/(rho_gas*mcross); 
  
  return lambda;
}



// porcedure not used
// confirmed JMB 4/14/17. called in epstein drag (not used)

void disk_dispersion(rho,time,t0,deltat,newrho)
double rho,time,t0,deltat,*newrho;
{

  if (time <= t0) 
      *newrho = rho;
      
  if ((time > t0) && (time <= t0+deltat))
      *newrho = rho*(1.0-(time-t0)/deltat);
      
  if (time > t0 + deltat)
      *newrho = 0.0;
  
}


double Normal(m, s)
   double m , s;
/* ========================================================================
 * Returns a normal (Gaussian) distributed real number.
 * NOTE: use s > 0.0
 *
 * Uses a very accurate approximation of the normal idf due to Odeh & Evans, 
 * J. Applied Statistics, 1974, vol 23, pp 96-97.
 * ========================================================================
 */
{ 
  const double p0 = 0.322232431088;     const double q0 = 0.099348462606;
  const double p1 = 1.0;                const double q1 = 0.588581570495;
  const double p2 = 0.342242088547;     const double q2 = 0.531103462366;
  const double p3 = 0.204231210245e-1;  const double q3 = 0.103537752850;
  const double p4 = 0.453642210148e-4;  const double q4 = 0.385607006340e-2;
  double u, t, p, q, z;

  u   = my_rand();
  if (u < 0.5)
    t = sqrt(-2.0 * log(u));
  else
    t = sqrt(-2.0 * log(1.0 - u));
  p   = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
  q   = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
  if (u < 0.5)
    z = (p / q) - t;
  else
    z = t - (p / q);
  return (m + s * z);
}
//used for making random intergers with log distribution
 void set_rand_log()
 {
   double p,factor;
   int i;
   p=0.999;
   factor=1/(log(1.0-p));
   ProbLog[0]=0.0;
   for (i=1;i<257;i++){
     ProbLog[i]=ProbLog[i+1] - factor*pow(p,i)/i;
   }
   ProbLog[256]=1.0;
 }

 //generates a random interger with a logarithmic distribution
 double rand_log()
 {
   double rand;
   int index,i;
   index=0;
   rand=my_rand();
   for (i=0;i<256;i++){
     if(ProbLog[i]<=rand &&ProbLog[i+1]>rand){
       index=i+1;
       break;
     }
   }
   return index;
 }

//finds gradient of turbulent potential at planet's location. Lambda is a global variable. 
 void get_Fstoc(x,y,time,r,Mcenter,s_dens,FS)
  double x,y,r,Mcenter, time,s_dens,*FS;
{
  double theta, Lambr,Lambp,rc,phic,t0,delT,zeta,sigma,omega,rc2,vkep,cs,sc_height;
  int j,m;
  theta=atan2(y,x);
  if (theta < 0.0){
      theta = theta+ TWOPI;
   }
  Lambr=0.0;
  Lambp=0.0;
  for (j=0;j<MODES;j++){
      m=Lambda[j];
      if (m < 7){
	rc=Lambda[j+MODES];
	phic=Lambda[j+2*MODES];
	t0=Lambda[j+3*MODES];
	delT=Lambda[j+4*MODES];
	zeta=Lambda[j+5*MODES];
	sigma=PI*rc/(4*m);
	rc2=sqrt(rc);
	omega=sqrt(G*G*Mcenter)/(rc*rc2);
	vkep=omega*rc;
	cs=TWOPI*rc/(m*delT);//true based on definition of delT
	sc_height=rc*cs/vkep;
	if (sigma > sc_height){
	  sigma=sc_height; //Using Wlad's idea that largest perturbation should equal scale height.
	}
	if (fabs(phic-theta) < (TWOPI*rc/m)){
	   Lambr=Lambr+(1+2*r*(r-rc)/(sigma*sigma))*zeta*exp(-(r-rc)*(r-rc)/(sigma*sigma))*cos(m*theta-phic-omega*(time-t0))*sin(PI*(time-t0)/delT);
	   Lambp=Lambp+m*zeta*exp(-(r-rc)*(r-rc)/(sigma*sigma))*sin(m*theta-phic-omega*(time-t0))*sin(PI*(time-t0)/delT);
	}
      }
   }
  Lambr=Lambr*gamma*64.0*s_dens*G*G/(TWOPI*TWOPI);//radial force
  Lambp=Lambp*gamma*64.0*s_dens*G*G/(TWOPI*TWOPI);//azimuthal force
  FS[0]= Lambr*cos(theta) - Lambp*sin(theta);//translated to x-y plane
  FS[1]=Lambr*sin(theta) + Lambp*cos(theta);
  FS[2]= 0.0; //Lambr*0.05;  //force out of plane
}

/*checks if any harmonics in turb potential have expired. LamIndex is an ordered list of all existing modes and the time at which they are set to expire*/
/* keeping this function in AU/day units JMB 4/13/17 */
void Lambda_Update(time,radii,aspR,Mcenter)
double time,*radii,*aspR,Mcenter;
{
  double r,phi,delT,zeta,t0,cs;
  int i,j,index,m;
  while (time - LamTime[0] > 0.0){
        m=rand_log();
        r=R_inner + (1000-R_inner-0.1)*my_rand();
        phi=TWOPI*my_rand();
        zeta=Normal(0.0,1.0);
	cs=get_speed_of_sound(r,radii,aspR,Mcenter); // AU per day
        delT=TWOPI * r/(m*cs);
	index = LamIndex[0];
	t0=time;
	Lambda[index]=m;
	Lambda[index+MODES]=r;
	Lambda[index+2*MODES]=phi;
	Lambda[index+3*MODES]=t0;
	Lambda[index+4*MODES]=delT;
	Lambda[index+5*MODES]=zeta;
	j=0;
	for (i=1;i<MODES;i++){
	  if (LamTime[i] < (t0+delT)){
		LamIndex[i-1]=LamIndex[i];
		LamTime[i-1]=LamTime[i];
		j=i;
	  }
	}
	LamIndex[j]=index;
	LamTime[j]=t0+delT;
  }
}
//called at beginning of simulation. 
void Generate_Lambda(r,aspR,Mcenter)
double *r,*aspR, Mcenter;
{
int i;
  for (i=0; i<MODES;i++){
    LamIndex[i]=i;
  }
  Lambda_Update(0.0001,r,aspR,Mcenter);
}
//finds minimum distance between two point particles deemed to be
//coliding...doesn't work
// called in collision function JMB 4/14/17  so it is used? does it work?
double r_min(x,v,m1,m2)
double *x,*v,m1,m2;
{
  double v2,v1,mu,r,epsilon,sma,e1,e2,e3,e,rmin,vx,vy,vz,x1,x2,x3;
  vx=v[0];
  vy=v[1];
  vz=v[2];
  x1=x[0];
  x2=x[1];
  x3=x[2];
  v2=vx*vx+vy*vy+vz*vz;
  v1=sqrt(v2);
  mu=G*G*(m1*m2);
  r=sqrt(x1*x1+x2*x2+x3*x3);
  epsilon=0.5*v2 - mu/r;
  sma=fabs(0.5*mu/epsilon);
  e1=v2*x1/mu - v1*vx*x1/mu - x1/r;
  e2=v2*x2/mu - v1*vy*x2/mu - x2/r;
  e1=v2*x3/mu - v1*vz*x3/mu - x3/r;
  e=sqrt(e1*e1+e2*e2+e3*e3);
  rmin=sma*(e-1);
  return rmin;
}
//returns ratio of kinetic energy of two coliding bodies to the
//binding energy of the resulting body 
// JMB: modify for BHs. used in collision function.
double Impact_energy(m1,m2,v,r)
double m1,m2,*v,r;
{
  double KE, BE,mu;
  mu=(m1*m2)/(m1+m2);
  KE=0.5*mu*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  BE= 0.5*G*G*(m1*m2)/r;
  return (KE/BE);

}
//no longer used (verified 4/14/17 JMB)
double Stochastic(F_old,h,sigma,Mcenter,sma)//,F_out)
double F_old,h,sigma,Mcenter,sma;
{
	double rms, tau,F_out;
    rms=PI*G_CGS*sigma*0.5/SDCONV*.000499;//end constant converts cgs to au/day
    //printf("sma = %lg\n",sma);
    //printf("rms = %lg\n",rms);
    tau=1.0*sqrt(Mcenter*G*G / (4*PI*sma*sma*sma));
    //printf("tau = %lg\n",tau);
    F_out = 0.0*(Normal(0.0,rms)+F_old*exp(-h*tau));
    //printf("stochastic = %lg\n",F_out);
    return F_out;
}


//returns semi-major axis of body]

void get_sma(cbody,ntot,bodies,radii,y,sma)
     double cbody;
     int ntot;
     double  *bodies, *radii, *y, *sma;
{
     double xx,yy,zz,vx,vy,vz,v2,v1,mu,r,epsilon;
     int i;
      for (i=0;i<ntot;i++) {
	xx=y[6*i];
	yy=y[6*i+1];
	zz=y[6*i+2];
	vx=y[6*i+3];
	vy=y[6*i+4];
	vz=y[6*i+5];
        v2=vx*vx+vy*vy+vz*vz;
        v1=sqrt(v2);
        mu=G*G*(bodies[i]*cbody);
        r=sqrt(xx*xx+yy*yy+zz*zz);
        epsilon=0.5*v2 - mu/r;
        sma[i]=fabs(0.5*mu/epsilon);
//	element(cbody,bodies[i],xx,yy,zz,vx,vy,vz,&ax,&ex,&inc,&omP,&OmN,&MAn);
//	sma[i]=ax;
      }
}



// procedure not used
// Verified JMB 4/14/17
double ALPHA(double radius)
{
 double eta = 0.0019*sqrt(radius);
 return sqrt(1.0 - 2.0*eta);
}

//not used :-)
// Verified JMB 4/14/17
double epstein_drag(rho_gas,p_rad,c_s)
double rho_gas, p_rad, c_s;
{ 
 double v_therm, t_stop;
 /* parameters*/
 double k_Boltzman, T, mu_H2, m_prot, c_s2;
 
  k_Boltzman = 1.3807e-16; /* Boltzman constant [cgs] */
       mu_H2 = 2.3;        /* mean molecular weight of the H2 molecule*/
      m_prot = 1.6726e-24; /* proton's mass [g] */
 
    c_s2 = c_s*c_s;
       T = c_s2*mu_H2*m_prot/k_Boltzman; 
 v_therm = sqrt(8.0*k_Boltzman*T/(M_PI*mu_H2*m_prot));
 
  //rho_gas = disk_dispersion(rho_gas,xbe,DISP_BEGIN,DISP_DURATION);
  disk_dispersion(rho_gas,xbe,DISP_BEGIN,DISP_DURATION,&rho_gas);
  if (rho_gas==0.0) {
      t_stop = 1.0e15;
  } else {
          t_stop = p_rad*RHO/(rho_gas*v_therm); /* in secundum */
          t_stop = t_stop/86400.0;
  }
  
  return t_stop;
}

// not used 
// Verified JMB 4/14/17

double stokes_drag(v_rel,rho_gas,p_rad,body_mass,c_sound,lambda)
double v_rel, rho_gas, p_rad, body_mass, c_sound, lambda;
{
 double mach, Reyn, C_drag, p_rad2, v_rel2, F_drag, nu, t_stop;


    v_rel =  v_rel*1.49e13/86400.0; /* [AU/day] --> [cm/sec] */
    body_mass = body_mass*SUN2GR;   /* Solar mass --> gramm  */
    
   /* calculation of the Reynolds number*/
   
   mach = v_rel/c_sound;       /* Mach number of the body*/
     nu = lambda*c_sound/3.0;  /* molecular viscosity of the gas */
   Reyn = 2.0*v_rel*p_rad/nu;  /* Reynolds number */
     
   /* calculation of the (Stokes) drag coefficient C_drag valid 
      for high Mach numbers, also in supersonic case            */
      
   if (mach < 0.5) {
       if (Reyn<=1.0) {
           C_drag = 24.0/Reyn;
       } else if (Reyn > 1.0 && Reyn <= 800.0) { 
                  C_drag = 24.0/pow(Reyn,0.6);
	      } else if (Reyn > 800.0) {
                         C_drag = 0.44;
	             }
   } else if (mach >= 0.5) {
              C_drag=1.096-0.638/(1.0+exp((mach-0.961)/0.127));
          }
   //printf("Mach %lg and Reynolds numbers %lg\n",mach,Reyn);
   p_rad2 = p_rad*p_rad;
   v_rel2 = v_rel*v_rel;
   
   //rho_gas = disk_dispersion(rho_gas,xbe,DISP_BEGIN,DISP_DURATION);
   disk_dispersion(rho_gas,xbe,DISP_BEGIN,DISP_DURATION,&rho_gas);
   
   F_drag = 0.5*C_drag*M_PI*p_rad2*rho_gas*v_rel2; 
   
   if (F_drag == 0.0) {
       t_stop = 1.0e15;
   } else {
            t_stop = body_mass*v_rel/F_drag; /* stopping time in sec */
            t_stop = t_stop/86400.0;         /* sec  -->  day */
   }
   
   return t_stop;
}

void read_op_files(indexs,sevens,six5s,sixs,five5s,fives,four5s,fours,three5s,threes,two5s,twos,one5s,ones,zero5s,zeros,pzero5s,pones)
/*read in opacity tables once, check them at every time step*/

double *indexs,*sevens, *six5s, *sixs,*five5s,*fives, *four5s, *fours, *three5s, *threes, *two5s, *twos,*one5s,*ones,*zero5s,*zeros,*pzero5s,*pones;
{
double index,seven, six5, six,five5,five, four5, four, three5, three, two5,two,one5,one,zero5,zero,pzero5,pone;
FILE *fin1,*fin2,*fin3,*fin4,*fin5,*fin6,*fin7,*fin8,*fin9,*fin10,*fin11,*fin12,*fin13,*fin14,*fin15,*fin16,*fin17;
FILE *fin18,*fin19,*fin20,*fin21,*fin22,*fin23,*fin24,*fin25,*fin26,*fin27,*fin28,*fin29,*fin30,*fin31,*fin32,*fin33,*fin34,*fin35;

fin1=fopen("./kfs/kappaaf-7.txt","r");
fin2=fopen("./kfs/kappaaf-6.5.txt","r");
fin3=fopen("./kfs/kappaaf-6.txt","r");
fin4=fopen("./kfs/kappaaf-5.5.txt","r");
fin5=fopen("./kfs/kappaaf-5.txt","r");
fin6=fopen("./kfs/kappaaf-4.5.txt","r");
fin7=fopen("./kfs/kappaaf-4.txt","r");
fin8=fopen("./kfs/kappaaf-3.5.txt","r");
fin9=fopen("./kfs/kappaaf-3.txt","r");
fin10=fopen("./kfs/kappaaf-2.5.txt","r");
fin11=fopen("./kfs/kappaaf-2.txt","r");
fin12=fopen("./kfs/kappaaf-1.5.txt","r");
fin13=fopen("./kfs/kappaaf-1.txt","r");
fin14=fopen("./kfs/kappaaf-0.5.txt","r");
fin15=fopen("./kfs/kappaaf0.txt","r");
fin16=fopen("./kfs/kappaaf0.5.txt","r");
fin17=fopen("./kfs/kappaaf1.txt","r");

fin18=fopen("./kfs/kappari-7.txt","r");
fin19=fopen("./kfs/kappari-6.5.txt","r");
fin20=fopen("./kfs/kappari-6.txt","r");
fin21=fopen("./kfs/kappari-5.5.txt","r");
fin22=fopen("./kfs/kappari-5.txt","r");
fin23=fopen("./kfs/kappari-4.5.txt","r");
fin24=fopen("./kfs/kappari-4.txt","r");
fin25=fopen("./kfs/kappari-3.5.txt","r");
fin26=fopen("./kfs/kappari-3.txt","r");
fin27=fopen("./kfs/kappari-2.5.txt","r");
fin28=fopen("./kfs/kappari-2.txt","r");
fin29=fopen("./kfs/kappari-1.5.txt","r");
fin30=fopen("./kfs/kappari-1.txt","r");
fin31=fopen("./kfs/kappari-0.5.txt","r");
fin32=fopen("./kfs/kappari0.txt","r");
fin33=fopen("./kfs/kappari0.5.txt","r");
fin34=fopen("./kfs/kappari1.txt","r");
fin35=fopen("./kfs/kappari_ind.txt","r");



int I=0;  
for (I;I<26;I++){
        fscanf(fin1,"%lg %lg\n",&index,&seven);
        indexs[I]=index;
        sevens[I]=seven;
        fscanf(fin2,"%lg %lg\n",&index,&six5);
        six5s[I]=six5;
        fscanf(fin3,"%lg %lg\n",&index,&six);
        sixs[I]=six;
        fscanf(fin4,"%lg %lg\n",&index,&five5);
        five5s[I]=five5;
        fscanf(fin5,"%lg %lg\n",&index,&five);
        fives[I]=five;
        fscanf(fin6,"%lg %lg\n",&index,&four5);
        four5s[I]=four5;
        fscanf(fin7,"%lg %lg\n",&index,&four);
        fours[I]=four;
        fscanf(fin8,"%lg %lg\n",&index,&three5);
        three5s[I]=three5;
        fscanf(fin9,"%lg %lg\n",&index,&three);
        threes[I]=three;
        fscanf(fin10,"%lg %lg\n",&index,&two5);
        two5s[I]=two5;
        fscanf(fin11,"%lg %lg\n",&index,&two);
        twos[I]=two;
	fscanf(fin12,"%lg %lg\n",&index,&one5);
        one5s[I]=one5;
        fscanf(fin13,"%lg %lg\n",&index,&one);
        ones[I]=one;
        fscanf(fin14,"%lg %lg\n",&index,&zero5);
        zero5s[I]=zero5;
        fscanf(fin15,"%lg %lg\n",&index,&zero);
        zeros[I]=zero;
        fscanf(fin16,"%lg %lg\n",&index,&pzero5);
        pzero5s[I]=pzero5;
        fscanf(fin17,"%lg %lg\n",&index,&pone);
        pones[I]=pone;
  }
int J=26;  
for (J;J<96;J++){
        fscanf(fin35,"%lg\n",&index);
        indexs[J]=index;
        fscanf(fin18,"%lg\n",&seven);
        sevens[J]=seven;
        fscanf(fin19,"%lg\n",&six5);
        six5s[J]=six5;
        fscanf(fin20,"%lg\n",&six);
        sixs[J]=six;
        fscanf(fin21,"%lg\n",&five5);
        five5s[J]=five5;
        fscanf(fin22,"%lg\n",&five);
        fives[J]=five;
        fscanf(fin23,"%lg\n",&four5);
        four5s[J]=four5;
        fscanf(fin24,"%lg\n",&four);
        fours[J]=four;
        fscanf(fin25,"%lg\n",&three5);
        three5s[J]=three5;
        fscanf(fin26,"%lg\n",&three);
        threes[J]=three;
        fscanf(fin27,"%lg\n",&two5);
        two5s[J]=two5;
        fscanf(fin28,"%lg\n",&two);
        twos[J]=two;
        fscanf(fin29,"%lg\n",&one5);
        one5s[J]=one5;
        fscanf(fin30,"%lg\n",&one);
        ones[J]=one;
        fscanf(fin31,"%lg\n",&zero5);
        zero5s[J]=zero5;
        fscanf(fin32,"%lg\n",&zero);
        zeros[J]=zero;
        fscanf(fin33,"%lg\n",&pzero5);
        pzero5s[J]=pzero5;
        fscanf(fin34,"%lg\n",&pone);
        pones[J]=pone;
  }

/*for (int j=0;j<97;j++){
        printf("%lg\n",indexs[j]);
//      printf("%lg\n",sevens[j]);
}*/
}



//new function to get kappa based more closely on the tables feeling 90% confident that all units are now correct
double get_kappa_agn(r,aspR_local,T,sig,BODY0,indexs,sevens,six5s,sixs,five5s,fives,four5s,fours,three5s,threes,two5s,twos,one5s,ones,zero5s,zeros,pzero5s,pones)
  double r,T,sig,BODY0,aspR_local;
  double *indexs,*sevens, *six5s, *sixs,*five5s,*fives, *four5s, *fours, *three5s, *threes, *two5s, *twos,*one5s,*ones,*zero5s,*zeros,*pzero5s,*pones;
  {
    double kk,H,rho,a,b;
    H=aspR_local*r*AU2CM;
    //printf("H= %g\n",H);
    rho=sig/(SDCONV*(2*H)); // convert from MSUN/AU^2 to g/cm^2
    //printf("rho = %g\n",rho);
    b = log10(rho/pow((T/1e6),3));
    a=log10(T);
    //printf("R: %g T: %g\n",b,a);
    if(a<2.9){
        kk=0.76;
    }
    else{
      int I=0;
      for (I; I<96;I++){
//      printf("pones=%g",pones[I]);
        if (b<=-6.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=sevens[I];
                        break;
                }
        }

        else if (b<-6.25){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=six5s[I];
                        break;
                }
        }

        else if (b<-5.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=sixs[I];
                        break;
                }
        }

        else if (b<-5.25){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=five5s[I];
                        break;
                }
        }

        else if (b<-4.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=fives[I];
                        break;
                }
        }

        else if (b<-4.25){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=four5s[I];
                        break;
                }
        }

        else if (b<-3.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=fours[I];
                        break;
                }
        }
        else if (b<-3.25){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=three5s[I];
                        break;
                }
        }

        else if (b<-2.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=threes[I];
                        break;
                }
        }

        else if (b<-1.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=twos[I];
                        break;
                }
        }

        else if (b<-1.25){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=one5s[I];
                        break;
                }
        }

        else if (b<-0.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=ones[I];
                        break;
                }
        }

        else if (b<-0.25){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=zero5s[I];
                        break;
                }
        }
        else if (b<0.25){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=zeros[I];
                        break;
                }
        }
        else if (b<0.75){
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=pzero5s[I];
                        break;
                }
        }

        else if (b<1.25){
        //      printf("indexs[I]=%g indexs[I-1]=%g\n",indexs[I],indexs[I-1]);
                if (a>=indexs[I]&&a<indexs[I-1]){
                        kk=pones[I];
                        break;
                }
        }
      }
}

//commenting out writing kappa file - takes too much memory
/*struct kap my_recordk;
my_recordk.a=a;
my_recordk.b=b;
my_recordk.rho=rho;
my_recordk.kk=kk;

fwrite(&my_recordk, sizeof(struct kap), 1,abkap);
//printf("kk= %g\n",kk);
fflush(abkap);*/
double kk_fin=pow(10.0,kk);
return kk_fin;
}


/*
// creating AGN opacity function JMB 5/30/17
// based on Iglesias & Rogers 96
// cutting a few corners, can be improved
// grid is more coarse here than in paper.
double get_kappa_agn_old(r,T,sig,BODY0)
     double r,T,sig,BODY0;
{
    double T1,T2,T3,T4,T5,T6,T7,T8,T9;
    double kk,H,rho,logR;

    T9 = pow(10.0,6);
    T8 = pow(10.0,5.75);
    T7 = pow(10.0,5.5);
    T6 = pow(10.0,5.25);
    T5 = pow(10.0,5.0);
    T4 = pow(10.0,4.75);
    T3 = pow(10.0,4.5);
    T2 = pow(10.0,4.25);
    T1 = pow(10.0,4.0);
    H=0//get_aspect_ratio(r,T,BODY0)*r;
    printf("radius= %g",r);
    rho=sig/(2*H);
    printf("rho: %g \n",rho);
    logR = log10(rho/pow((T/1e6),3));
    //printf("LogR: %g \n",logR);

	
    if (T <= T1){  // T1 = 1e4 K
	if (logR <= -8.0){
	    kk = pow(10.0,-0.518);
	}
	else if (logR <= -7.0){
	    kk = pow(10.0,-0.513);
	    }
	else if (logR <= -6.0){
	    kk = pow(10.0,-0.411);
	    }
	else if (logR <= -5.0){
	    kk = pow(10.0,-0.034);
	    }
	else if (logR <= -4.0){
	    kk = pow(10.0,0.651);
	    }
	else if (logR <= -3.0){
	    kk = pow(10.0,1.262);
	    }
	else if (logR <= -2.0){
	    kk = pow(10.0,1.587);
	    }
	else if (logR <= -1.0){
	    kk = pow(10.0,1.806);
	    }
	else if (logR <= 0.0){
	    kk = pow(10.0,2.04);
	    }
	else {
	    kk = pow(10.0,2.348);
	    }
    }
    else if (T <= T2){  // T2 = 1e4.25 K
	if (logR <= -8.0){
            kk = pow(10.0,-0.474);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.439);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,-0.353);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,-0.096);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,0.479);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,1.337);
            }
        else if (logR <= -2.0){
            kk = pow(10.0,2.315);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,3.223);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,3.794);
            }
        else {
            kk = pow(10.0,4.141);
            }
	}
    else if (T <= T3){ // T3 =1e 4.5 K
	if (logR <= -8.0){
            kk = pow(10.0,-0.415);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.380);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,-0.258);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,0.013);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,0.529);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,1.404);
            }
        else if (logR <= -2.0){
            kk = pow(10.0,2.456);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,3.548);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,4.524);
            }
        else {
            kk = pow(10.0,5.202);
            }
	}
   else if (T <= T4){ // T4 =1e4.75 K
	if (logR <= -8.0){
            kk = pow(10.0,-0.410);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.370);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,-0.281);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,0.028);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,0.712);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,1.664);
            }
       else  if (logR <= -2.0){
            kk = pow(10.0,2.688);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,3.674);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,4.561);
            }
        else {
            kk = pow(10.0,5.203);
            }
	}
   else if (T <= T5){ // T5 =1e5 K
	if (logR <= -8.0){
            kk = pow(10.0,-0.339);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.315);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,-0.239);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,-0.004);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,0.565);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,1.416);
            }
        else if (logR <= -2.0){
            kk = pow(10.0,2.416);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,3.434);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,4.290);
            }
        else {
            kk = pow(10.0,4.788);
            }
	}
 else if (T <= T6){ // T6 =1e5.25 K
     //printf("we are in T6\n");
	if (logR <= -8.0){
            kk = pow(10.0,-0.215);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.089);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,0.091);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,0.385);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,0.841);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,1.531);
            }
        else if (logR <= -2.0){
            kk = pow(10.0,2.405);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,3.195);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,3.768);
            }
        else {
	    //printf("we are in the else\n");
	    kk = pow(10.0,4.166);
            }
	}

else if (T <= T7){ // T7 =1e5.5 K
	if (logR <= -8.0){
            kk = pow(10.0,-0.457);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.423);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,-0.306);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,-0.014);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,0.571);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,1.349);
            }
       else  if (logR <= -2.0){
            kk = pow(10.0,2.083);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,2.643);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,3.083);
            }
        else {
            kk = pow(10.0,3.508);
            }
	}

else if (T <= T8){ // T8 =1e5.75 K
	if (logR <= -8.0){
            kk = pow(10.0,-0.467);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.460);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,-0.440);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,-0.355);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,0.053);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,0.986);
            }
        else if (logR <= -2.0){
            kk = pow(10.0,1.780);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,2.269);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,2.671);
            }
        else {
            kk = pow(10.0,3.103);
            }
	}

else if (T <= T9){ // T9 =1e6 K
	if (logR <= -8.0){
            kk = pow(10.0,-0.464);
	    }
        else if (logR <= -7.0){
            kk = pow(10.0,-0.457);
            }	
	else if (logR <= -6.0){
            kk = pow(10.0,-0.441);
            }
	else if (logR <= -5.0){
            kk = pow(10.0,-0.383);
            }
        else if (logR <= -4.0){
            kk = pow(10.0,-0.137);
            }
        else if (logR <= -3.0){
            kk = pow(10.0,0.585);
            }
        else if (logR <= -2.0){
            kk = pow(10.0,1.464);
            }
        else if (logR <= -1.0){
            kk = pow(10.0,2.039);
            }
        else if (logR <= 0.0){
            kk = pow(10.0,2.441);
            }
        else {
            kk = pow(10.0,2.840);
            }
	}
 else if (T > T9){ // S&G doesn't go this high.  this is a placeholder.
	      // JMB 5/31/17
     kk = pow(10.0,2.840);
     }
    //printf("kappa = %lg \n",kk);
 return kk;
    }


// needs AGN disk upgrade badly JMB 4/14/17
double get_kappa(r,T,sig,BODY0)
     double r,T,sig,BODY0;  
{
  double T1,T2,T3,T4,T5,T6,T7,T8,T9;
  double kk,a,b,kap,H,rho,logkk,logk;
  kk=0.0;
  T1=132.0;
  T2=170.0;
  T3=375.0;
  T4=390.0;
  T5=580.0;
  T6=680.0;
  T7=960.0;
  T8=1570.0;
  T9=3730.0;

      if (T <= T1){
            kap=2e-4; 
            a=0.0;
            b= 2.1;
            kk=kap*pow(T,b);
        }
      else if ((T > T1) && (T < T2)){
            kap=3.0;
            a=0.0;
            b=-0.01;
            kk=kap*pow(T,b);
        }
      else if ((T > T2) && (T < T3)){
            kap=0.01;
            a=0.0;
            b= 1.1;
            kk=kap*pow(T,b);
	}
      else if ((T > T3) && (T < T4)) {
            kap=5e4;
            a=0;
            b=-1.5;
            kk=kap*pow(T,b);
	}
      else if ((T > T4) && (T < T5)) {
            kap=0.1;
            a=0.0;
            b= 0.7;
            kk=kap*pow(T,b);
	}
      else if ((T > T5) && (T < T6)) {
            kap=2e15;
            a=0.0;
            b=-5.2;
            kk=kap*pow(T,b);
	}
      else if ((T > T6) && (T < T7)) {
            kap=0.02;
            a=0.0;
            b= 0.8;
            kk=kap*pow(T,b);
	}
      else if ((T > T7) && (T < T8)) {
            logk=81.3010;
            a=1.0;
            b=-24.0;
            H=get_aspect_ratio(r,T,BODY0)*r;
            rho=sig/(2*H);
            logkk=logk+a*log10(rho)+b*log10(T);
            kk=pow(10,logkk);
            kap=1e33;
	}
      else if ((T > T8)){  // this is true for AGN disks.  likely need
			   // to update.
            kap=1e-8;
            a=2.0/3;
            b=3.0;	  
            H=get_aspect_ratio(r,T,BODY0)*r; 
					   // 1.0 was MCenter which is
					   // not 1.0
            rho=sig/(2*H);
            kk=kap*pow(rho,a)*pow(T,b);
	}
    return kk;
}
*/

double Opacity(r,aspR_local,sig,T,BODY0,aindexs,asevens,asix5s,asixs,afive5s,afives,afour5s,afours,athree5s,athrees,atwo5s,atwos,aone5s,aones,azero5s,azeros,apzero5s,apones)
     double r,sig,T,BODY0,aspR_local; // sig, T transposed.  fixed JMB 4/11/17
     double *aindexs,*asevens,*asix5s,*asixs,*afive5s,*afives,*afour5s,*afours,*athree5s,*athrees,*atwo5s,*atwos,*aone5s,*aones,*azero5s,*azeros,*apzero5s,*apones;
{
  double omega, kap,taup,taueff,theta, cv;

  //kap=get_kappa_agn_old(r,T,sig,BODY0);  // opacity tables for planet disks
  kap=get_kappa_agn(r,aspR_local,T,sig,BODY0,aindexs,asevens,asix5s,asixs,afive5s,afives,afour5s,afours,athree5s,athrees,atwo5s,atwos,aone5s,aones,azero5s,azeros,apzero5s,apones);  // use opacity tables for AGN

/*
  struct sigR my_recordsr;
  my_recordsr.sig=sig;
  my_recordsr.r=r;

  fwrite(&my_recordsr, sizeof(struct sigR), 1, surf);
  fflush(surf);
 */
  //printf("unconverted sig= %g\n", sig);
  sig = sig/SDCONV;  // now in g/cm2 - AS moved this up here b/c kap is in cm2/g so need taup to be dimensionless
  //printf("converted sig= %g r= %g\n", sig,r*AU2CM);
  //printf("kappa= %g\n", kap);
  taup=0.5*kap*sig;
  //printf("taup: %g\n",taup);
  taueff=0.375*taup+sqrt(3.0)*0.25+0.25/taup;
  //printf("taueff: %g\n",taueff);



  //printf("sig: %g\n",sig);
  //printf("temp: %g\n",T);
  cv = 12.5e7; // JMB 4/12/17 updated to monatomic gas 8.314e7*3/2; 
  r = r*AU2CM;
  // BODY0 = 1e8  solarmasses
  // r = in  AU from input file
  // T = 807294
  // sig = 5518.8  g/cm2  // reasonable
  //omega=G*sqrt(BODY0)/(24*3600*r*sqrt(r));  // missing mass factor, presumably
				// because Mass = 1 Msun?
                                // BODY0 Mass added 4/17/17
  //printf("omega first: %g\n",omega);
  //printf("sig: %g\n",sig);
  //printf("r %g\n",r);
  omega = sqrt(G_CGS*BODY0*SUN2GR/(pow(r,3.0)));// cgs units!

  theta=cv*sig*omega*taueff/(4.0*PI*5.67e-5 *T*T*T);
  //printf("theta: %g\n",theta);
return theta;
}

// displacement array?  what is this? maybe displacement of bodies due to
// other bodies?  this calculates the distances between bodies.
// 4/14/17 JMB
void Set_Disp_Array(nbody,y)
int nbody;
double *y;
{
  int i,j, arr_ind;
  double dx,dy,dz,d2;
  for (i=0;i<(nbody-1);i++) {
    for (j=i+1;j<(nbody);j++) {
      dx=y[6*j]-y[6*i];
      dy=y[6*j+1]-y[6*i+1];
      dz=y[6*j+2]-y[6*i+2];
      d2=dx*dx+dy*dy+dz*dz;
      arr_ind=2*(i*(2*nbody-(i+3))+2*j-2);
      Disp_Array[arr_ind]=dx;
      Disp_Array[arr_ind+1]=dy;
      Disp_Array[arr_ind+2]=dz;
      Disp_Array[arr_ind+3]=d2;
      //printf("dist %d \n",arr_ind);
    }
  }
}

// set up an array of distances between all bodies and all other
// bodies (and their components!)  ???  JMB 4/18/17
void Set_Dyn_Disp(nbody)
int nbody;
{
  int i,j,arr_ind;
  for (i=0;i<(nbody-1);i++) {
    for (j=i+1;j<(nbody);j++) {
      arr_ind=2*(i*(2*nbody-(i+3))+2*j-2);
      Dyn_Disp[arr_ind]=Disp_Array[arr_ind];
      Dyn_Disp[arr_ind+1]=Disp_Array[arr_ind+1];
      Dyn_Disp[arr_ind+2]=Disp_Array[arr_ind+2];
      Dyn_Disp[arr_ind+3]=Disp_Array[arr_ind+3];
    }
  }
}

void Set_temp_Disp(disp,ntot)
int ntot;
double *disp;
{
  int i;
  for (i=0;i<(ntot*(ntot-1)*2);i++) {
    disp[i]=Dyn_Disp[i];
  }
}

void Update_Dyn_Disp(disp,ntot)
int ntot;
double *disp;
{
  int i;
  for (i=0;i<(ntot*(ntot-1)*2);i++) {
    Dyn_Disp[i]=disp[i];
  }
}

void Set_equal_Disp(disp1,disp2,ntot)
int ntot;
double *disp1,*disp2;
{
  int i;
  for (i=0;i<(ntot*(ntot-1)*2);i++) {
    disp1[i]=disp2[i];
  }
}

void Update_temp_Disp(disp,nbody,y,dt)
int nbody;
double *y,*disp,dt;
{
  int i,j,arr_ind;
  double vix,viy, viz,vjx,vjy,vjz,dx,dy,dz;
  for (i=0;i<(nbody-1);i++) {
    vix=y[6*i]*dt;
    viy=y[6*i+1]*dt;
    viz=y[6*i+2]*dt;
    for (j=i+1;i<(nbody);i++) {
      vjx=y[6*j]*dt;
      vjy=y[6*j+1]*dt;
      vjz=y[6*j+2]*dt;
      arr_ind=2*(i*(2*nbody-(i+3))+2*j-2);
      disp[arr_ind]+=vix+vjx;
      disp[arr_ind+1]+=viy+vjy;
      disp[arr_ind+2]+=viz+vjz;
      dx=disp[arr_ind];
      dy=disp[arr_ind+1];
      dz=disp[arr_ind+2];
      disp[arr_ind +3]=dx*dx+dy*dy+dz*dz;
      }
  }
}

void Update_Disp(disp1,disp2,ntot,dy,dt)
int ntot;
double *disp1, *disp2, *dy, dt;
{
  int i, j ,arr,i1,i2,i3,j1,j2,j3;
  double x,y,z,mag;
  for (i=0;i<(ntot-1);i++) {
    i1=6*i; i2=i1+1; i3=i2+1;
    for (j=i+1;j<ntot;j++)  {
      j1=6*j; j2=j1+1; j3=j2+1;
      arr=2*(i*(2*ntot-(i+3))+2*j-2);
      x=disp2[arr]+(dy[j1]-dy[i1])*dt;
      y=disp2[arr+1]+(dy[j2]-dy[i2])*dt;
      z=disp2[arr+2]+(dy[j3]-dy[i3])*dt;
      mag=x*x+y*y+z*z;
      disp2[arr]=disp1[arr];
      disp2[arr+1]=disp1[arr+1];
      disp2[arr+2]=disp1[arr+2];
      disp2[arr+3]=disp1[arr+3];
      disp1[arr]=x;
      disp1[arr+1]=y;
      disp1[arr+2]=z;
      disp1[arr+3]=mag;
      Dyn_Disp[arr]=x;
      Dyn_Disp[arr+1]=y;
      Dyn_Disp[arr+2]=z;
      Dyn_Disp[arr+3]=mag;      
    }
  }
}

/* calculation of the mutual gravitational forces in a star/Sun centered 
   coordinate system */
// verified that units
void Gravitational_Force(nbody,cbody,body,y,force,mindist,i_dmin,j_dmin)
int nbody; 
double cbody, *body, *y, *force, *mindist;
int *i_dmin, *j_dmin;
{
 int i, i1, i2, i3, i4, i5, i6, imin;
 int j, j1, j2, j3, j4, j5, j6, jmin;
 int arr_ind;
 double cent_acc[3],xi,yi,zi,ri2,rin3,gsum,gcent,xij,yij,zij,rij2,rijn3,gi,gj,grav;

      double  dmin2 = 10000.0;
               imin = 0;
               jmin = 0;
grav=G*G;
cent_acc[0]=0.0;
cent_acc[1]=0.0;
cent_acc[2]=0.0;
 for (i=0;i<nbody;i++) {   // loop over all orbiting bodies

      i1 = 6*i; i2=i1+1; i3=i2+1; i4=i3+1; i5=i4+1; i6=i5+1;
      xi=y[i1];
      yi=y[i2];
      zi=y[i3];
      ri2=xi*xi+yi*yi+zi*zi;
      rin3 = 1.0/(ri2*sqrt(ri2));  // 1/r3
      gsum=grav*body[i]*rin3;
      cent_acc[0]+=gsum*xi;    // summing the force(!) on the central body
      cent_acc[1]+=gsum*yi;
      cent_acc[2]+=gsum*zi;

      force[i1] = y[i4];  // setting force = velocity?
      force[i2] = y[i5]; 
      force[i3] = y[i6];

      gcent=grav*rin3*cbody;  // G*M_BH/r3

      force[i4] = - gcent*xi;  // force on orbiting body
      force[i5] = - gcent*yi; 
      force[i6] = - gcent*zi;
  }
  for (i=0;i<nbody;i++) {
    force[i*6+3]+=-cent_acc[0];  // cent_acc is a force.
    force[i*6+4]+=-cent_acc[1];
    force[i*6+5]+=-cent_acc[2];
  }
  // printf("GForce \n");
  for (i=0;i<(nbody-1);i++) {
    i4=6*i+3; i5=i4+1; i6=i5+1;
    for (j=i+1;j<nbody;j++)  {
      j4=6*j+3; j5=j4+1; j6=j5+1;
      arr_ind=2*(i*(2*nbody-(i+3))+2*j-2);
      xij=Dyn_Disp[arr_ind];
      yij=Dyn_Disp[arr_ind+1];
      zij=Dyn_Disp[arr_ind+2];
      rij2=Dyn_Disp[arr_ind+3];
      rijn3=1.0/(rij2*sqrt(rij2));  // 1/ distance between bodies ^3
      gi=grav*body[i]*rijn3;   
      gj=grav*body[j]*rijn3;
      force[i4]+=gj*xij;  // forces from bodies on each other
      force[i5]+=gj*yij;
      force[i6]+=gj*zij;
      force[j4]-=gi*xij;
      force[j5]-=gi*yij;
      force[j6]-=gi*zij;
      if (rij2 < dmin2) {
         dmin2 = rij2;
	 imin = i;
	 jmin = j;
      }  
    }
  }
      *mindist = sqrt(dmin2);
      *i_dmin  = imin;
      *j_dmin  = jmin;
}

   

/* adding additional forces: originating from migration, and drag - this is not implemented in
   this verion                                                                               */
/*  pretty sure it is implemented JMB 4/13/17 */ 
void N_Body_Force(Body0,Ntot,Body,radii,time,sigma,alpha,temp,beta,Y,F_stoc,Force)
double  Body0,time;
int     Ntot;
double *Body, *radii, *sigma, *alpha, *temp, *beta,*Y, *F_stoc, *Force;
{
  int    i, i1, i2, i3, i4, i5, i6, imin, jmin;
  
  double gsum, ri, ri2, r_vr, dmin, dcrit,bcap,
         h, h2, h4, coefc, alfa, vcirc, anomf, 
         vxcirc, vycirc, tau_a, tau_e, tau_i, const_taua, const_taue, const_taui, 
         torque_const,torque,s_dens, dless_torque,T,alpha_local,beta_local, rrel, ratio;
  double Torque[2], Damp[3], xrel[3], vrel[3]; 
    
    imin = 0;
    jmin = 0;
    dmin =10000.0;  /* where does this come from?  JMB 6/1/17
		       Looks like it's overwritten in the next call
		       though */

    /* Gravitational force between all the bodies */
    Gravitational_Force(Ntot,Body0,Body,Y,Force,&dmin,&imin,&jmin);


    /* Migration of the planets: --> IFLAGMIGR = 1 
                    no migration --> IFLAGMIGR = 0 
       using the approach of Cresswell & Nelson (2008)                  */

    if (MIGRATION_TYPE1 == 1 ) {
	
	
	
        for (i=0;i<Ntot;i++) {
             //gsum = G*G*(Body0 + Body[i]);

             i1 = 6*i;
             i2 = i1+1;
             i3 = i2+1;
             i4 = i3+1;
             i5 = i4+1;
             i6 = i5+1;
             ri = sqrt(Y[i1]*Y[i1] + Y[i2]*Y[i2] + Y[i3]*Y[i3]);
             //ri2= ri*ri;
             //r_vr = Y[i1]*Y[i4] + Y[i2]*Y[i5] + Y[i3]*Y[i6];
	     
             //vcirc  =  sqrt(gsum/ri);
             //anomf  =  atan2(Y[i2],Y[i1]);
             //vxcirc = -vcirc*sin(anomf);
             //vycirc =  vcirc*cos(anomf);
             
	     //get_local_gas_properties(ri,radii,alpha,sigma,temp,beta, &s_dens,&T,&alpha_local,&beta_local);

             //get_migration_damping(Body0,Body[i],Y[i1],Y[i2],Y[i3],Y[i4],Y[i5],Y[i6],ri,s_dens,T,alpha_local,beta_local,Torque,Damp);

	     //get_Fstoc(Y[i1],Y[i2],time,ri,Body0,s_dens, FS);
	     //element(Body0,Body[i],Y[i1],Y[i2],Y[i3],Y[i4],Y[i5],Y[i6],&sma,&exc,&inc,&omegaP,&OmegaN,&kozan);
	     
	     //h=get_aspect_ratio(ri,T,Body0);
	     //h2=h*h;
	     //h4=h2*h2;
	     //torque_const = G*G/(Body0*h2);
	     //const_taue = h4*pow(Body0,1.5)/0.78;  
	     //const_taui = h4*pow(Body0,1.5)/0.544;

	     //damping_timescales(ri,s_dens,alpha_local,beta_local,exc,inc,h,ADIABATIC,&tau_e,&tau_i,&dless_torque);
             
	     //torque = dless_torque*torque_const*Body[i]*s_dens;
	    // tau_a = tau_a*const_taua/Body[i];
	     //tau_e = tau_e*Body0*sqrt(Body0)/Body[i];
             //tau_i = tau_i*Body0*sqrt(Body0)/Body[i];

		//printf("               SIGMA = %lg\n",s_dens);
		//printf("torque = %lg\n",Torque[1]);
		//printf("F = %lg\n",Force[i4]);
                                   //  ++alpha		                       	++Damp                    ++Stoch                           
             Force[i4] = Force[i4] +  F_stoc[i];
             Force[i5] = Force[i5]+ F_stoc[i+Ntot];
                                   //                              ++Damp                     ++Damp    ++Stoch
	     Force[i6] = Force[i6]+ F_stoc[i+2*Ntot];


             //coefc = 0.5/tau_a + 1.0/tau_e;
             //alfa = 1.0/tau_e/coefc;

             //Force[i4] = Force[i4] - coefc*(Y[i4]-alfa*vxcirc);
             //Force[i5] = Force[i5] - coefc*(Y[i5]-alfa*vycirc);
             //Force[i6] = Force[i6] - coefc* Y[i6];

             
        }
    }
    /* migration is added */


    /*first time we go through binary protocal*/
      //commenting this out for now

      I_BVECTOR[0]=Y[6*imin];
      I_BVECTOR[1]=Y[6*imin+1];
      I_BVECTOR[2]=Y[6*imin+2];
      I_BVECTOR[3]=Y[6*imin+3];
      I_BVECTOR[4]=Y[6*imin+4];
      I_BVECTOR[5]=Y[6*imin+5];
      J_BVECTOR[0]=Y[6*jmin];
      J_BVECTOR[1]=Y[6*jmin+1];
      J_BVECTOR[2]=Y[6*jmin+2]; 
      J_BVECTOR[3]=Y[6*jmin+3];
      J_BVECTOR[4]=Y[6*jmin+4];
      J_BVECTOR[5]=Y[6*jmin+5];
      bcap=rhill(Body[imin],Body[jmin],Body0);
      if (dmin<bcap){
	 xrel[0]=I_BVECTOR[0]-J_BVECTOR[0];
	 xrel[1]=I_BVECTOR[1]-J_BVECTOR[1];
         xrel[2]=I_BVECTOR[2]-J_BVECTOR[2];
         vrel[0]=I_BVECTOR[3]-J_BVECTOR[3];
         vrel[1]=I_BVECTOR[4]-J_BVECTOR[4];
         vrel[2]=I_BVECTOR[5]-J_BVECTOR[5];
         rrel=sqrt(xrel[0]*xrel[0]+xrel[1]*xrel[1]+xrel[2]*xrel[2]);
	 ratio=Impact_energy(Body[imin],Body[jmin],&vrel,rrel); 
	 //printf("r=%g dmin=%g KE/BE=%g \n",rrel,dmin,ratio);
	 if (ratio<1){
         	COLLISION=1;
        	IMIN = int_min(imin,jmin);
        	JMIN = int_max(imin,jmin);
	 }
      }
  /* 
      else{
	BINFORM=0;
      }
*/

    //collision is still commented out because I am now using the old binary method because it
    //uses the old binary conditions just changed BINFORM=1 to COLLISION=1
    /* check for collision                                                         */
    /* critical distance = few times (= BETA) of the sum of the radii of the closest bodies */
    //dcrit = BETA*(radius(imin,Body) + radius(jmin,Body));
    //DCRITICAL = dcrit;
    /* criterion for collision: the distance of the closest bodies < critical distance */
    /*if (dmin < dcrit) {

    // collision --> IFLAGCOLL = 1
        COLLISION = 1;  // setting this = 1 triggers collision function later
	I_VECTOR[0]=Y[6*imin];
	I_VECTOR[1]=Y[6*imin+1];
	I_VECTOR[2]=Y[6*imin+2];
	I_VECTOR[3]=Y[6*imin+3];
	I_VECTOR[4]=Y[6*imin+4];
	I_VECTOR[5]=Y[6*imin+5];
	J_VECTOR[0]=Y[6*jmin];
	J_VECTOR[1]=Y[6*jmin+1];
	J_VECTOR[2]=Y[6*jmin+2];
	J_VECTOR[3]=Y[6*jmin+3];
	J_VECTOR[4]=Y[6*jmin+4];
	J_VECTOR[5]=Y[6*jmin+5];
	IMIN = int_min(imin,jmin);
	JMIN = int_max(imin,jmin);
    } */
   
}


/* calculation of the new velocity of the body, elimination of the colliding body 
   collisions: always accretional completely non-elastic, the lowest index body will be 
   the new one the masses are added         
                                            */
void print_positions(nbody,y,index,F_out,xbe)  // prints x,y,z, vx,vy,vz,
					   // and fx,fy,vz I
					   // believe...  JMB
					   // 5/29/17
double *y;
double *F_out;
int nbody,*index;
double xbe; // time?
{
 int i;
    for (i=0; i<nbody; i++){
	fprintf(fout4," %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg \n",index[i],xbe,y[6*i],y[6*i+1],y[6*i+2],y[6*i+3],y[6*i+4],y[6*i+5],F_out[i],F_out[i+nbody],F_out[i+2*nbody]);
          fflush(fout4);
	  /*	  fprintf(fout5,"   \n",stuff);
		  fflush(fout5); */
	}
}



void bpos(y,ycom,truePOS,nbod)
double *y, *ycom, *truePOS;
int nbod;
{
int icol,jcol;

  icol=IMIN;
  jcol=JMIN;

  /* just set all values of ybnew and body_bnew as the same*/
  int i=0;
  for (i; i<nbod; i++) {
         truePOS[6*i+0] = y[6*i+0];
         truePOS[6*i+1] = y[6*i+1];
         truePOS[6*i+2] = y[6*i+2];
         truePOS[6*i+3] = y[6*i+3];
         truePOS[6*i+4] = y[6*i+4];
         truePOS[6*i+5] = y[6*i+5];
      }

  truePOS[6*icol+0]=y[6*icol+0]-(0.25)*ycom[6*icol+0]; 
  truePOS[6*icol+1]=y[6*icol+1]-(0.25)*ycom[6*icol+1]; 
  truePOS[6*icol+2]=y[6*icol+2]-(0.25)*ycom[6*icol+2]; 

  truePOS[6*jcol+0]=y[6*jcol+0]-(0.25)*ycom[6*jcol+0];
  truePOS[6*jcol+0]=y[6*jcol+1]-(0.25)*ycom[6*jcol+1];
  truePOS[6*jcol+0]=y[6*jcol+2]-(0.25)*ycom[6*jcol+2]; 

}

//this method creates an array that replaces the two bodies that are in a
//binary with their center of mass vectors and masses at their indexes
//to be used in the migration, dampening and stochastic force calculations  
void binary(time,nbod,body,y,index,y_bnew,body_bnew,y_before)
double time; int nbod; double *body, *y, *y_before; double *y_bnew; double *body_bnew; int *index;
{
 int icol, jcol, i, i1;
 double tmp1, tmp2, x[3],v[3],rmin,r,v1c[3],v2c[3],V1,V2,ratio;

      /* IMIN and JMIN external static integers are given values in N_Body_Force() */

      icol = IBMIN;
      jcol = JBMIN;

/*      x[0]=I_BVECTOR[0]-J_BVECTOR[0];
      x[1]=I_BVECTOR[1]-J_BVECTOR[1];
      x[2]=I_BVECTOR[2]-J_BVECTOR[2];
      v[0]=I_BVECTOR[3]-J_BVECTOR[3];
      v[1]=I_BVECTOR[4]-J_BVECTOR[4];
      v[2]=I_BVECTOR[5]-J_BVECTOR[5];

      r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      rmin=r_min(x,v,body[icol],body[jcol]);
*/
      tmp1 = body[icol]/(body[icol]+body[jcol]);
      tmp2 = body[jcol]/(body[icol]+body[jcol]);
/*
      v1c[0]=I_VECTOR[3] - tmp1*I_VECTOR[3]-tmp2*J_VECTOR[3];
      v1c[1]=I_VECTOR[4] - tmp1*I_VECTOR[4]-tmp2*J_VECTOR[4];
      v1c[2]=I_VECTOR[5] - tmp1*I_VECTOR[5]-tmp2*J_VECTOR[5];
      v2c[0]=J_VECTOR[3] - tmp1*I_VECTOR[3]-tmp2*J_VECTOR[3];
      v2c[1]=J_VECTOR[4] - tmp1*I_VECTOR[4]-tmp2*J_VECTOR[4];
      v2c[2]=J_VECTOR[5] - tmp1*I_VECTOR[5]-tmp2*J_VECTOR[5];
      V1=sqrt(v1c[0]*v1c[0]+v1c[1]*v1c[1]+v1c[2]*v1c[2]);
      V2=sqrt(v2c[0]*v2c[0]+v2c[1]*v2c[1]+v2c[2]*v2c[2]);
      ratio=Impact_energy(body[icol],body[jcol],V1,V2);

*/

      /* just set all values of ybnew and body_bnew as the same*/
      for (i=0; i<nbod; i++) {
             y_bnew[6*i+0] = y[6*i+0];
             y_bnew[6*i+1] = y[6*i+1];
             y_bnew[6*i+2] = y[6*i+2];
             y_bnew[6*i+3] = y[6*i+3];
             y_bnew[6*i+4] = y[6*i+4];
             y_bnew[6*i+5] = y[6*i+5];
             body_bnew[i] = body[i];
      }


      /*set the values for y_bnew as the center of mass values*/

      y_bnew[6*icol+0] = y[6*icol+0]*tmp1 + y[6*jcol+0]*tmp2;
      y_bnew[6*icol+1] = y[6*icol+1]*tmp1 + y[6*jcol+1]*tmp2;
      y_bnew[6*icol+2] = y[6*icol+2]*tmp1 + y[6*jcol+2]*tmp2;
      y_bnew[6*icol+3] = y[6*icol+3]*tmp1 + y[6*jcol+3]*tmp2;
      y_bnew[6*icol+4] = y[6*icol+4]*tmp1 + y[6*jcol+4]*tmp2;
      y_bnew[6*icol+5] = y[6*icol+5]*tmp1 + y[6*jcol+5]*tmp2;

      y_bnew[6*jcol+0] = y[6*icol+0]*tmp1 + y[6*jcol+0]*tmp2;
      y_bnew[6*jcol+1] = y[6*icol+1]*tmp1 + y[6*jcol+1]*tmp2;
      y_bnew[6*jcol+2] = y[6*icol+2]*tmp1 + y[6*jcol+2]*tmp2;
      y_bnew[6*jcol+3] = y[6*icol+3]*tmp1 + y[6*jcol+3]*tmp2;
      y_bnew[6*jcol+4] = y[6*icol+4]*tmp1 + y[6*jcol+4]*tmp2;
      y_bnew[6*jcol+5] = y[6*icol+5]*tmp1 + y[6*jcol+5]*tmp2;

      /* the masses of the colliding bodies are simply added */
      body_bnew[icol] = body[icol] + body[jcol];
      body_bnew[jcol] = body[icol] + body[jcol];

}


   
void collision(time,nbod,body,y,index,new_nbod,y_new,body_new,index_new,y_before,radii,aspR,alpha,sigma,temp,beta,M0)
double time; int nbod; double *body, *y, *y_before; int *new_nbod; double *y_new; double *body_new; int *index; int *index_new,M0;
double *radii,*alpha,*sigma,*aspR,*temp,*beta;
{
 int icol, jcol, i, i1;
 double tmp1, tmp2, x[3],v[3],rmin,r,v1c[3],v2c[3],V1,V2,ratio,thard,new_a;
 double s_dens,T,aspR_local,alpha_local,beta_local;

      /* IMIN and JMIN external static integers are given values in N_Body_Force() */
      
      icol = IMIN; 
      jcol = JMIN;
   
      x[0]=I_BVECTOR[0]-J_BVECTOR[0];//y_before[6*icol]-y_before[6*jcol];
      x[1]=I_BVECTOR[1]-J_BVECTOR[1];
      x[2]=I_BVECTOR[2]-J_BVECTOR[2];
      v[0]=I_BVECTOR[3]-J_BVECTOR[3];
      v[1]=I_BVECTOR[4]-J_BVECTOR[4];
      v[2]=I_BVECTOR[5]-J_BVECTOR[5];

      r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      rmin=r_min(x,v,body[icol],body[jcol]);

      tmp1 = body[icol]/(body[icol]+body[jcol]);
      tmp2 = body[jcol]/(body[icol]+body[jcol]);

      v1c[0]=I_BVECTOR[3] - tmp1*I_BVECTOR[3]-tmp2*J_BVECTOR[3];
      v1c[1]=I_BVECTOR[4] - tmp1*I_BVECTOR[4]-tmp2*J_BVECTOR[4];
      v1c[2]=I_BVECTOR[5] - tmp1*I_BVECTOR[5]-tmp2*J_BVECTOR[5];
      v2c[0]=J_BVECTOR[3] - tmp1*I_BVECTOR[3]-tmp2*J_BVECTOR[3];
      v2c[1]=J_BVECTOR[4] - tmp1*I_BVECTOR[4]-tmp2*J_BVECTOR[4];
      v2c[2]=J_BVECTOR[5] - tmp1*I_BVECTOR[5]-tmp2*J_BVECTOR[5];
      V1=sqrt(v1c[0]*v1c[0]+v1c[1]*v1c[1]+v1c[2]*v1c[2]);
      V2=sqrt(v2c[0]*v2c[0]+v2c[1]*v2c[1]+v2c[2]*v2c[2]);
      ratio=Impact_energy(body[icol],body[jcol],&v,r);
      
      // confusing time units, pretty sure time and xki are in years
      // not days?  JMB 6/1/17
      printf("collision between bodies %ld %ld at time %lg days\n",IMIN,JMIN,time); // usedtobe time/365.25

      fprintf(fout3,"collision between bodies %ld %ld at time %lg days\n",index[icol],index[jcol],xki); // same, xki/365.25
      fprintf(fout3,"Current distance %lg \n",r);
      fprintf(fout3,"minimum distance %lg \n",rmin);
      fprintf(fout3,"relative velocity %lg %lg %lg \n",v[0],v[1],v[2]);
      fprintf(fout3,"Kinetic Energy/Binding Energy %lg \n",ratio);
      
      /* up to the colliding body with higher index 'jcol' the components are unaltered, 
         except the case of the colliding body with index 'icol' */

      for (i=0; i<jcol; i++) {

	 if (i!=icol) {
	     y_new[6*i+0] = y[6*i+0];
	     y_new[6*i+1] = y[6*i+1];
	     y_new[6*i+2] = y[6*i+2];
	     y_new[6*i+3] = y[6*i+3];
	     y_new[6*i+4] = y[6*i+4];
	     y_new[6*i+5] = y[6*i+5];
	     body_new[i] = body[i]; 
             index_new[i] = index[i];
	 }
      }      	
	
      /* calculation of the coordinates and velocities the newly formed body
         which inherits the index 'icol'                                     */	
      
      y_new[6*icol+0] = y[6*icol+0]*tmp1 + y[6*jcol+0]*tmp2;
      y_new[6*icol+1] = y[6*icol+1]*tmp1 + y[6*jcol+1]*tmp2;
      y_new[6*icol+2] = y[6*icol+2]*tmp1 + y[6*jcol+2]*tmp2;
      y_new[6*icol+3] = y[6*icol+3]*tmp1 + y[6*jcol+3]*tmp2;
      y_new[6*icol+4] = y[6*icol+4]*tmp1 + y[6*jcol+4]*tmp2;
      y_new[6*icol+5] = y[6*icol+5]*tmp1 + y[6*jcol+5]*tmp2;
      
      /* the masses of the colliding bodies are simply added */ 
      body_new[icol] = body[icol] + body[jcol];
      index_new[icol] = index[icol];    

 
      /* elimination of the 'jcol'th body from the list
         by shifting the indices in the array           */
      for (i=jcol; i<nbod-1; i++) {

	  i1 = i + 1;

	  y_new[6*i+0] = y[6*i1+0];
	  y_new[6*i+1] = y[6*i1+1];
	  y_new[6*i+2] = y[6*i1+2];
	  y_new[6*i+3] = y[6*i1+3];
	  y_new[6*i+4] = y[6*i1+4];
	  y_new[6*i+5] = y[6*i1+5];

	  body_new[i] = body[i1];
	  index_new[i] = index[i1];
      }
      /* the number of bodies is decreased by 1 */
      *new_nbod = nbod - 1;
      
      new_a=sqrt(y_new[6*icol+0]*y_new[6*icol+0]+y_new[6*icol+1]*y_new[6*icol+1]+y_new[6*icol+2]*y_new[6*icol+2]);
      get_local_gas_properties(new_a,radii,aspR,alpha,sigma,temp,beta, &s_dens,&T,&aspR_local,&alpha_local,&beta_local);
      thard=t_hard(T,s_dens,M0,body_new[icol],new_a,r);
      fprintf(fout3,"T_hard %lg \n",thard);

      fflush(fout3);
}



/* here are the numerical integrator procedures, which do not need modifications */

void mmid(central, ntot, bodies, radii, sigma, alpha, temp, beta, y, F_stoc, dydx, xs, htot, nstep, yout)
     double central;
     int ntot;
     double *bodies;
     double *radii, *sigma, *alpha,*temp,*beta;
     double *y,*F_stoc,*dydx;
     double xs, htot;
     long nstep;
     double *yout;
     
{
  long n, i;
  double x, swap, h2, h;
  double dispm[ntot*(ntot-1)*2], dispn[ntot*(ntot-1)*2];
  vektor ym, yn;

  h = htot / nstep;
  Set_Dyn_Disp(ntot);
  Set_temp_Disp(dispm,ntot);
  Set_temp_Disp(dispn,ntot);
  Update_Disp(dispn,dispm,ntot,dydx,h);
  
  for (i = 0; i < nv; i++) {
    ym[i] = y[i];
    yn[i] = y[i] + h * dydx[i];
  }
  x = xs + h;
  N_Body_Force(central,ntot,bodies,radii,x,sigma,alpha,temp,beta,yn,F_stoc,yout);
  h2 = 2 * h;
  for (n = 2; n <= nstep; n++) {
    Update_Disp(dispn,dispm,ntot,yout,h2);

    for (i = 0; i < nv; i++) {
      swap = ym[i] + h2 * yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    N_Body_Force(central,ntot,bodies,radii,x,sigma,alpha,temp,beta,yn,F_stoc,yout);
  }
  for (i = 0; i < nv; i++)
    yout[i] = 0.5 * (ym[i] + yn[i] + h * yout[i]);
}


#define ncol            7

void rzextr(iest, xest, yest, yz, dy, nuse)
long iest;
double xest;
double *yest, *yz, *dy;
long nuse;
{
  long m1, k, j;
  double yy, v, c, b, b1, ddy;
  double fx[ncol];

  glx[iest - 1] = xest;
  if (iest == 1) {
    for (j = 0; j < nv; j++) {
      yz[j] = yest[j];
      gld[j][0] = yest[j];
      dy[j] = yest[j];
    }
    return;
  }
  if (iest < nuse)
    m1 = iest;
  else
    m1 = nuse;
  for (k = 1; k < m1; k++)
    fx[k] = glx[iest - k - 1] / xest;
  for (j = 0; j < nv; j++) {
    yy = yest[j];
    v = gld[j][0];
    c = yy;
    gld[j][0] = yy;
    for (k = 1; k < m1; k++) {
      b1 = fx[k] * v;
      b = b1 - c;
      if (b != 0) {
        b = (c - v) / b;
        ddy = c * b;
        c = b1 * b;
      } else
        ddy = v;
      v = gld[j][k];
      gld[j][k] = ddy;
      yy += ddy;
    }
    dy[j] = ddy;
    yz[j] = yy;

  }
}

#undef ncol
#define imax            11
#define nuse            7

#define one             1.0e0
#define shrink          0.95e0
#define grow            1.2e0


void BSmethod(holdp,period,htemp2print,central,ntot,bodies,radii,sigma,alpha,temp,beta,y, F_stoc,dydx, x, htry, eps, loc_yki, loc_hki, loc_hkovki)
double central;
int ntot;
int period;
double * bodies;
double *radii, *sigma, *alpha, *temp, *beta;
double *y, *F_stoc, *dydx;
double x, htry, eps;
double *loc_yki;
double *loc_hki;
double *loc_hkovki;
int *holdp;
double *htemp2print;
{
  //vektor y, dydx, yscal;
  long j, i, k;
  double xsav, xest, h, errmax;
  vektor ysav, dysav, yseq, yerr;
  long nseq[imax];
  double hdid, hnext, TEMP;
  int hcount=*holdp;
  //memcpy(y, y_, sizeof(vektor));
  //memcpy(dydx, dydx_, sizeof(vektor));
  //memcpy(yscal, yscal_, sizeof(vektor));
  
  nseq[0] = 2;
  nseq[1] = 4;
  nseq[2] = 6;
  nseq[3] = 8;
  nseq[4] = 12;
  nseq[5] = 16;
  nseq[6] = 24;
  nseq[7] = 32;
  nseq[8] = 48;
  nseq[9] = 64;
  nseq[10]= 96;
  h = htry;
  xsav = x;
  for (i = 0; i < nv; i++) {
    ysav[i] = y[i];
    dysav[i] = dydx[i];
  }
  for (;;) {
    for (i = 1; i <= imax; i++) {
      mmid(central,ntot,bodies,radii,sigma,alpha,temp,beta,ysav,F_stoc, dysav, xsav, h, nseq[i - 1], yseq);
      TEMP = h / nseq[i - 1];
      xest = TEMP * TEMP;
      rzextr(i, xest, yseq, y, yerr, (long)nuse);
      errmax = 0.0;
      for (j = 0; j < nv; j++) {
        if (errmax < fabs(yerr[j]))
          errmax = fabs(yerr[j]);
      }
      errmax /= eps;
      if (errmax < one) {
        hdid = h;
        *loc_hki = hdid;
        if (i == nuse)
          hnext = h * shrink;
        else {
          if (i == nuse - 1){
            //printf("grow: %lg\n",hnext);
            hnext = h * grow;
	 }
          else{
		  //printf("other else: %lg\n",hnext);
		  hnext = h * nseq[nuse - 2] / nseq[i - 1];
          }
	}
        *loc_hkovki = hnext;
        for (k=0; k<nv; k++)
             loc_yki[k] = y[k];
        goto _L99;
      }
    }
    h = 0.25 * h;
    if ((imax - nuse) / 2 > 0) {
      for (i = 1; i <= (imax - nuse) / 2; i++)
        h /= 2;
    }
     
    if (x + h == x)
    {
      printf("imax=%d nuse=%d",imax,nuse);
      printf("The BS procedure stopped.\n");
      printf("No convergence for any values of the time step!\n");
      getchar();
    }
  }
_L99: ;

struct hcon my_recordh;
my_recordh.h=h;

fwrite(&my_recordh, sizeof(struct hcon), 1, hconv);

fflush(hconv);
*htemp2print=h;

//saving period
if (hcount>period){
  PRINT_TIME=1;
  hcount=0;
}
else{
  PRINT_TIME=0;
  hcount++;
  //printf("hcount= %d\n",hcount);
}
*holdp=hcount;
}

#undef imax
#undef nuse
#undef one
#undef shrink
#undef grow


void Integrator(holdp,period,htemp2print,cbody,ntot,bodies,radii,sigma,aspR,alpha,temp,beta,y,x,h,sma,yo,xo,hn,F_out,aindexs,asevens,asix5s,asixs,afive5s,afives,afour5s,afours,athree5s,athrees,atwo5s,atwos,aone5s,aones,azero5s,azeros,apzero5s,apones)
     double  cbody;  /*mass of the central body*/
     int     ntot;   /*number of bodies except the central one*/
     int period;
     double *bodies; /*array of the masses*/
     double *radii;
     double *sigma;
     double *aspR;
     double *alpha;
     double *temp;
     double *beta;
     double *y;    /*input vector - initial condtions of the ODE*/
     double  x;    /*input time   */
     double  h;    /*initial time step*/
     double *sma;  //semi-major axes of all bodies
     double *yo;   /*output vector - solution of the ODE*/
     double *xo;   /*time = input time + time step used, this maybe smaller than "h" */
     double *hn;   /*estimated timestep for the next step of the integration*/
     double *F_out; //stochastic force used in this step
     double *aindexs,*asevens,*asix5s,*asixs,*afive5s,*afives,*afour5s,*afours,*athree5s,*athrees,*atwo5s,*atwos,*aone5s,*aones,*azero5s,*azeros,*apzero5s,*apones;
     double *htemp2print;
     int *holdp;
{

 vektor dydx;

 double hki, tol,lhn, ri, s_dens,FS[3],T,aspR_local,alpha_local,beta_local,Torque[2],Damp[3];
 int cv,i,i1,i2,i3,i4,i5,i6,j,i_int,k;

        for (i=0;i<ntot;i++){
             ri=sqrt(y[6*i]*y[6*i]+y[6*i+1]*y[6*i+1]+y[6*i+2]*y[6*i+2]);
             i1=6*i; i2=i1+1; i3=i2+1; i4= i3+1; i5=i4+1; i6=i5+1;
             /*STOCHASTIC FORCES!!!!!*/
             //s_dens = sigma[i_int]+(ri-radii[i_int])*(sigma[i_int+1]-sigma[i_int])/(radii[i_int+1]-radii[i_int]);

             get_local_gas_properties(ri,radii,aspR,alpha,sigma,temp,beta, &s_dens,&T,&aspR_local,&alpha_local,&beta_local);
             // the following will be converted to cgs units
             if (UPDATE_OPACITY==1){
                 //printf("T %g Sigma %g \n",T,s_dens);
                 OPACITY[i]=Opacity(ri,aspR_local,s_dens,T,BODY0,aindexs,asevens,asix5s,asixs,afive5s,afives,afour5s,afours,athree5s,athrees,atwo5s,atwos,aone5s,aones,azero5s,azeros,apzero5s,apones);
                 //printf("OPACITY[i] %g\n",OPACITY[i]);
                }

             get_migration_damping(cbody,bodies[i],y[i1],y[i2],y[i3],y[i4],y[i5],y[i6],ri,s_dens,T,aspR_local,alpha_local,beta_local,Torque,Damp,OPACITY[i]);
             //s_dens=get_sigma(ri,radii,sigma);
             //printf("Torque: %g %g \n",Torque[0],Torque[1]);
             //printf("Damp: %g %g %g \n",Damp[0],Damp[1],Damp[2]);
             get_Fstoc(y[6*i],y[6*i+1],x,ri,cbody,s_dens,FS);
             //printf("FS: %g %g %g \n",FS[0],FS[1],FS[2]);
             F_out[i]=FS[0]+Torque[0]+Damp[0];
             F_out[i+ntot]=FS[1]+Torque[1]+Damp[1];
             F_out[i+2*ntot]=FS[2]+Damp[2];
             // printf("%lg %lg %lg\n",F_out[i],F_out[i+ntot],F_out[i+2*ntot]);
             //  printf("********************** \n");
        }
UPDATE_OPACITY=0;


 tol = 1.0e-13;
 //printf("TIME = %lg\n",x);
 //printf("Integrate \n");
 Set_Disp_Array(ntot,y);
//printf("inthalf \n");
 Set_Dyn_Disp(ntot);
 //printf("Integratetwo \n");
 N_Body_Force(cbody,ntot,bodies,radii,x,sigma,alpha,temp,beta,y,F_out,dydx);
 BSmethod(holdp,period,htemp2print,cbody, ntot, bodies, radii, sigma, alpha, temp,beta, y, F_out, dydx, x, h, tol, yo, &hki, &lhn);

 *xo = x + hki;
 *hn = lhn;
}




/* prints out the results - orbital elements at given period, of gradually */
//now prints every x timesteps
void print_elements(oldhtemp,timeold,time,nbod,cbody,bodies,y_vekt,gradual,index)
double timeold; double time; int nbod; double cbody; 
double *bodies; double *y_vekt; int *index;
int gradual;
double *oldhtemp;
{
  int i,j;
  double x,y,z,vx,vy,vz,mass,period;
  double ax,ex,inc,MAn,omP,OmN,Ex,Ey,vmag,rmag,mu,r,ecc2,inc2;

  // MAn = mean anomaly;  
  // OmN = OmegaN = Longitude of asc. node
  // omP = OmegaP = argument of the perihelion

  struct rec my_record;  
  if (!ptr_myfile)
  {
     printf("Unable to open file!");
  }
 
/*
  switch (gradual) {
          case 0:
          period = peri;
          break;

          case 1:
             //time = (time+EPOCH);
          if (time <= 1000000.0)
              period = 100.0;
          if (time > 1000000.0 && time <= 2000000.0)
              period = 100.0;
          if (time > 2000000.0 && time <= 3000000.0) 
              period = 100.0;
          if (time > 3000000.0 && time <= 4000000.0) 
              period = 100.0;
          break;
  }
*/

//  if (fmod(time,period) <= time-timeold) {
     for (i=0;i<nbod;i++) {
        
           x  = y_vekt[6*i+0];
           y  = y_vekt[6*i+1];
           z  = y_vekt[6*i+2];
           vx = y_vekt[6*i+3];
           vy = y_vekt[6*i+4];
           vz = y_vekt[6*i+5];
	   r=sqrt(x*x+y*y+z*z);
           mass = bodies[i];
	   j= index[i];
	   mu=G*G*(cbody+mass);
	   vmag=sqrt(vx*vx+vy*vy+vz*vz);
	   rmag=sqrt(x*x+y*y+z*z);
	   Ex=(vmag*vmag*x - x*vx*vx - y*vx*vy)/mu -x/rmag;
	   Ey=(vmag*vmag*y- x*vx*vy - y*vy*vy)/mu -y/rmag;

          element(cbody,bodies[i],x,y,z,vx,vy,vz,&ax,&ex,&inc,&omP,&OmN,&MAn); //worried there is something fishy going on here
	
	   omP=acos(Ex/sqrt(Ex*Ex + Ey*Ey));
	   if (omP < 0.0){
		omP=omP+TWOPI;
	    }
	   omP=omP*TODEG;

	/*convert i and j to floats to make binary structure all floats*/
	float ni=(float)i;
	float nj=(float)j;
	
	/*set up structure for binary printing*/
	my_record.i=ni;
	my_record.j=nj;
	my_record.mass=mass;
	my_record.time=time;
	my_record.ax=ax;
	my_record.ex=ex;
	my_record.inc=inc;
	//printf("inc Fifteen = %g \n",inc); *JIA*9/25*
	my_record.MAn=MAn;
	my_record.omP=omP;
	my_record.OmN=OmN;
	my_record.x=x;
	my_record.y=y;
	my_record.z=z;
	my_record.vx=vx;
	my_record.vy=vy;
	my_record.vz=vz;
	my_record.oldhtemp=*oldhtemp;

	//write output binary file
	fwrite(&my_record, sizeof(struct rec), 1, ptr_myfile);

     /* 
	old mehtod for printing time
	
	   if (time <=400.0*365.2)
	 	fprintf(foutb,"%d %d %lg %5.1f %lg %lg %lg %lg %lg %lg %1g %1g %1g %1g %1g %1g %lg\n",i,j,mass,time,ax,ex,inc,MAn,omP,OmN,x,y,z,vx,vy,vz,*oldhtemp);
          if (time > 400.0*365.2 && time <= 800.0*365.2)
	 	fprintf(foutm1,"%d %d %lg %5.1f %lg %lg %lg %lg %lg %lg %1g %1g %1g %1g %1g %1g %lg\n",i,j,mass,time,ax,ex,inc,MAn,omP,OmN,x,y,z,vx,vy,vz,*oldhtemp);
          if (time > 800.0*365.2 && time <=1200.0*365.2)
	 	fprintf(foutm2,"%d %d %lg %5.1f %lg %lg %lg %lg %lg %lg %1g %1g %1g %1g %1g %1g %lg\n",i,j,mass,time,ax,ex,inc,MAn,omP,OmN,x,y,z,vx,vy,vz,*oldhtemp);
          if (time > 1200.0*365.2 && time <= 2000.0*365.2)
	 	fprintf(foute,"%d %d %lg %5.1f %lg %lg %lg %lg %lg %lg %1g %1g %1g %1g %1g %1g %lg\n",i,j,mass,time,ax,ex,inc,MAn,omP,OmN,x,y,z,vx,vy,vz,*oldhtemp);
          // fprintf(fout,"%d %d %lg %5.1f %1g %1g %1g\n",i,j,mass,time,x,y,z);
          
	*/

	   fflush(ptr_myfile);
           fflush(foutb);
           fflush(foutm1);
           fflush(foutm2);
           fflush(foute);
	   
          // printf("time = %lg\n",time);
      // }
       //printf("time = %lg\n",time);
 } 
 *oldhtemp=500;         
}


/* detects if one body leaves the simulation domain between RESCMIN and RESCMAX */

void detect_escape(time,y,nbod,iesc,iflage,index)
double time; double *y; int nbod; int *iesc; int *iflage; int *index;
{
  int i, i1, i2, i3;
  double ri;

   for (i=0;i<nbod;i++) {
   
       i1 = 6*i; i2 = 6*i+1; i3 = 6*i+2;
       ri = sqrt(y[i1]*y[i1] + y[i2]*y[i2] + y[i3]*y[i3]);
       //printf("tavolsag = %lg\n",ri);
       //getchar();
       
       if ((ri < RESCMIN) || (ri > RESCMAX)) {
          *iflage = 1; *iesc = i;
	  printf("escape of body %d at %lg years\n",i,time);
	  fprintf(fout3,"escape of body %d at %lg years\n",index[i],time);
	  fflush(fout3);
	  break;
       }
   }   
           
}

//reads in the file listing all the filenames of the gas profiles--each filename is the year of the profile
void init_time(f_time, tgass)
FILE *f_time;
int *tgass;
{
   int ttt, I;
	for (I=0; I < NTS; I++) {
	 	fscanf(f_time, "%d",&ttt);
		tgass[I]=ttt;
       }
	printf("%d%d",tgass[0],tgass[10]);
}
//reads in first gas profile
void init_radii(timesarr, radarr, sig, alph,temp,beta)
double *radarr, *sig, *alph,*temp,*beta;
int *timesarr;
{
  FILE *fileRad;
  int I;
  double surfdens, r,tt,expal,expbe;
  char file1 [ FILENAME_MAX ];
	sprintf(file1, "%s%d%s",gasdirect, 0, ".dat");
	fileRad = fopen(file1,"r");
	for (I=0; I<NRAD; I++) {
            fscanf(fileRad,"%lg %lg %lg %lg %lg",&r,&expal,&expbe,&surfdens,&tt);
	    radarr[I]=(r/1.49e13); // convert cm to AU 
	    sig[I]=surfdens*SDCONV;  // convert to Msun/AU2
	    alph[I]=expal;
	    temp[I]=tt;
	    beta[I]=expbe;
	}
	fclose(fileRad);
	

/*	for (I = 1; I < NRAD-1; I++) 
       		alph[I]  = -1.0*(radarr[I]/sig[I]*(sig[I+1]-sig[I-1])/(radarr[I+1]-radarr[I-1])); 
       
   	alph[0]  = alph[1];
   	alph[NRAD-1]= alph[NRAD-2]; */
	
	for (I = 0; I < NRAD; I++) {
         //printf("%lg %lg %lg\n",radii[i],sigma[i],alpha[i]);
		fprintf(fout_sd," %lg %lg %lg\n", radarr[I],sig[I],alph[I]);
		fflush(fout_sd);
        }  
   printf("Hello\n");
}
//reads in next gas profile.
void read_density_data(nr,itgas,tgas,sigma_1,sigma_2,alpha_1,alpha_2,temp_1,temp_2,beta_1,beta_2)
int nr, itgas;
int *tgas;
double *sigma_1, *sigma_2,*alpha_1,*alpha_2,*beta_1,*beta_2,*temp_1,*temp_2;
{
 double r,surfdens,expal,expbe,tt;
 char file1 [ FILENAME_MAX ], file2 [FILENAME_MAX ];
 long i;
 int t1, t2;
 FILE *fileSig1, *fileSig2;
   t1=tgas[itgas];
   t2=tgas[itgas+1];
   //   printf("files: %s %s \n",file1,file2);
   sprintf(file1, "%s%d%s",gasdirect, 0, ".dat");
   sprintf(file2, "%s%d%s",gasdirect, 0, ".dat");
   fileSig1 = fopen(file1,"r");
   for (i = 0; i < nr; i++)  {
	    fscanf(fileSig1,"%lg %lg %lg %lg %lg",&r,&expal,&expbe,&surfdens,&tt);
	    sigma_1[i]=surfdens*SDCONV;	   
	    beta_1[i]=expbe;
	    alpha_1[i]=expal;
	    temp_1[i]=tt;
       }
   //   printf("itgas: %d  NTS-1*1-0.1: %g \n ",itgas,(NTS-1)*1.0-0.1);
   if (itgas*1.0 > (NTS-1)*1.0-0.1) { // this is the last time step
				      // JMB
       //       printf("two \n");
	for (i = 0; i < nr; i++){
		sigma_2[i]=sigma_1[i];
		alpha_2[i]=alpha_1[i];
		beta_2[i]=beta_1[i];
		temp_2[i]=temp_1[i];
	}
	//	printf("two and a half \n");
	fileSig2 = fopen(file2,"r");
   } else {
   fileSig2 = fopen(file2,"r");
   for (i = 0; i < nr; i++)  {
            fscanf(fileSig2,"%lg %lg %lg %lg %lg",&r,&expal,&expbe,&surfdens,&tt);
	    sigma_2[i]=surfdens*SDCONV;	  
	    beta_2[i]=expbe;
	    alpha_2[i]=expal;
	    temp_2[i]=tt; 
       }
   }
   //   printf("three \n");
   if (itgas > 0) {
	t1=tgas[itgas-1];
	sprintf(file1, "%d%s", t1, ".dat");
	fclose(fileSig1);
	if (itgas*1.0 < (NTS-1)*1.0-0.1) { // if on last step don't
					   // need to close this file?
	    fclose(fileSig2);
	    }
       }
}

//recalculates current profile by linear interpolation between two nearest timesteps
void update_density(xbe,t1,t2,nr,radii,sigma_1,sigma_2,alpha_1,alpha_2,temp_1,temp_2,beta_1,beta_2,sigma,alpha,temp,beta)
int  t1, t2, nr;
double xbe;
double *radii, *sigma_1, *sigma_2,*alpha_1,*alpha_2,*temp_1,*temp_2,*beta_1,*beta_2, *sigma, *alpha,*beta,*temp;
{
 long i;
  double ratio;
  ratio=(xbe-t1*365.25)/(t2*365.25-t1*365.25);
  for (i = 0; i < nr; i++){
	sigma[i]=(sigma_1[i]+(sigma_2[i]-sigma_1[i])*ratio);
	alpha[i]=(alpha_1[i]+(alpha_2[i]-alpha_1[i])*ratio);
	beta[i]=(beta_1[i]+(beta_2[i]-beta_1[i])*ratio);
	temp[i]=(temp_1[i]+(temp_2[i]-temp_1[i])*ratio);
  } 

 /* for (i = 1; i < nr-1; i++) {
       alpha[i]  = -1.0*(radii[i]/sigma[i]*(sigma[i+1]-sigma[i-1])/(radii[i+1]-radii[i-1])); 
       if(sigma[i] = 1e-7)
	  alpha[i]=0.0;
   }   
   alpha[0]  = alpha[1];
   alpha[nr-1]= alpha[nr-2]; 
   
   for (i = 0; i < nr; i++) {
         //printf("%lg %lg %lg\n",radii[i],sigma[i],alpha[i]);
	fprintf(fout_sd,"%lg %d  %lg %lg %lg\n", xbe, t1, radii[i],sigma[i],alpha[i]);
	fflush(fout_sd);
      }*/
}

/* creates an array of times seperated by 0.4 Myr from time 0.4 to the end of the run (t_final) */
 void body_creation_times(seed,tau,t_final,times)
unsigned int seed;
double tau, t_final, *times;
{
  int i;
  double t_old, t_new;

  for (i = 0; i< 1500; i++)
        times[i] = 0.0;

  i = 0;
  t_old = 0.0;
  srand(seed);
  do {
      
      t_new = t_old +4e5*365.2;
      
      times[i] = t_new;
      printf("%d %lg\n",i,times[i]/365);
      t_old = t_new;
      i++;

  } while (t_new < t_final);
}

/* body creation*/

void create_a_body(n_body,sma_min,inc_max,ecc_max,body_list,y_vekt,x,y,z,vx,vy,vz,m_cr)
int  n_body;
double sma_min, inc_max, ecc_max, *body_list, *y_vekt,*x, *y, *z, *vx, *vy, *vz,*m_cr;

 {

  int i,i1,i2,i3,i4,i5,i6,j;
  double m0,m,sma,ecc,in,om,on,man,r,ri,distance,factor,RHill,
         mnew,lx,ly,lz,lvx,lvy,lvz;
  
  factor = 3.0;
 
  mnew = 3;
_L88:; 
  j++;
  sma = sma_min;
  ecc =  ecc_max* my_rand();
  in  =  inc_max* my_rand();
  //printf("in Sixteen = %g \n",in); *JIA*9/25*
  man =  360.0  * my_rand();
  om  =  360.0  * my_rand();
  on  =  0.0;
  
  coordinate(BODY0,mnew,sma,ecc,in,om,on,man,&lx,&ly,&lz,&lvx,&lvy,&lvz);

  r = sqrt(lx*lx + ly*ly + lz*lz);

  for (i = 0; i < n_body; i++) {
       i1 = 6*i;
       i2 = i1+1;
       i3 = i2+1;
       i4 = i3+1;
       i5 = i4+1;
       i6 = i5+1;
       ri = sqrt(y_vekt[i1]*y_vekt[i1] + y_vekt[i2]*y_vekt[i2] + y_vekt[i3]*y_vekt[i3]);
       RHill = pow((mnew+body_list[i])/(3.0*BODY0),0.33)*(r+ri)/2.0;
       distance = sqrt((lx-y_vekt[i1])*(lx-y_vekt[i1])+(ly-y_vekt[i2])*(ly-y_vekt[i2])+(lz-y_vekt[i3])*(lz-y_vekt[i3]));

       if (distance < factor * RHill) {
           goto _L88;
       }
  }
   *m_cr=mnew; *x = lx; *y = ly; *z = lz; *vx = lvx; *vy = lvy; *vz = lvz;
 }


/* The main program                                                                                       */
/**********************************************************************************************************/
int main(void)

{  
  
  double ar,oldhtemp,xx, yy, zz, vxx, vyy, vzz, m, ax, ex, inc, omP, OmN, MAn, Inc_max, creation_period,
         final_creation_time, dummy,m_cr, Sma_min, Ecc_max, axis, ecce, incl, ompe, omno, mean, r, surfdens;
  double Radii[NRAD], Alpha[NRAD], Alpha_1[NRAD], Alpha_2[NRAD],Sigma[NRAD], Sigma_1[NRAD], Sigma_2[NRAD], Temp[NRAD]; 
  double Temp_1[NRAD], Temp_2[NRAD],Beta[NRAD],Beta_1[NRAD], Beta_2[NRAD],Times[1500],aspR[NRAD], 
    F_stoc_out[NStart*3],x[3],d2,dmin,dcrit; 
  int    I, i, i1, iesc, int_dummy, CONTINUATION, GRADUAL, t, itgas, int_index,GAS_TORQUE,icol,jcol,saving_period;
  int    tgas[NTS], index[MAXNBOD], index_new[MAXNBOD];
  unsigned int SEED1, SEED2, SEED3;
  long nr_of_times;
  char file1 [ FILENAME_MAX ],timefile [ FILENAME_MAX ];
  vektor sma;
  double aindexs[97],asevens[97], asix5s[97], asixs[97],afive5s[97],afives[97],afour5s[97],afours[97],athree5s[97],athrees[97],atwo5s[97],atwos[97],aone5s[97],aones[97],azero5s[97],azeros[97],apzero5s[97],apones[97];  
  /*********************************************************************************************************/
  /* file declarations                                                                                     */
  /*********************************************************************************************************/
  fin1 = fopen("RETRO_for_Poster/elements.dat","r");
  fin2 = fopen("hor.dat","r");
  sprintf(timefile,"%s%s",gasdirect,"time.dat");
  ptr_myfile=fopen("RETRO_for_Poster/bout.bin","wb");
  //surf=fopen("surfR.bin","wb");
  fin_time = fopen(timefile,"r");
  //foutb= fopen("outputb.dat","w");
  //foutm1= fopen("outputm1.dat","w");
  //foutm2= fopen("feb1/outputm2.dat","w");
  //foute= fopen("feb1/outpute.dat","w");
  fout_sd = fopen("RETRO_for_Poster/av_density.dat","w");
  fout3= fopen("RETRO_for_Poster/analysis.dat","w");
  //fout4= fopen("xy.dat","w");
  fout5= fopen("RETRO_for_Poster/everything.bin","wb");
  abkap= fopen("RETRO_for_Poster/abkappa.bin","wb"); 
  hconv=fopen("RETRO_for_Poster/hconv.bin","wb");
  /*********************************************************************************************************/
  
  NBOD = NStart;         /* number of bodies */
  nv   = 6*NBOD;
  SEED1 = 23;        /* random seed for the random number generator to calculate inclinations */
  SEED2 = 26;        /* random seed for the random number generator to calculate the body creation times */
  SEED3 = 54;        /*                                                                   positions      */
  oldhtemp=500; 
  htemp2print=500; 
  Sma_min = 1000.0;    /* NOT IN     */
  Ecc_max = 0.07;   /*  USE !!!   */
  Inc_max = 200; /// *JIA9/21* Inc_max set to 190
  ///Inc_max = 3.25;  /* maximum value of the orbital inclinations. For 2D runs Inc_max = 0.0 setting should be used */
  RHO = 2.0;        /* Bulk density of planetesimals in g/cm3 */
  BODY0 = 1.0e8;      /* Mass of the central body in Solar mass */

  
  ADINDEX = 1.6666667;     /* adiabatic index, if for type I migration an adiabatic disk model is used 
                        only used if the variable "ADIABATIC" is set to 1, see "ADIABATIC".                        */
  
  SIGMA_0 = 7.65e-5; /* surface density at 1 AU in AU*M_star^(-2)
			units in the disk model used for drag or
			migration 
		        JMB 4/16: does not appear to be in use*/ 


  EPOCH = 0.0; //34500.0*365.25; /* beginning time of continuation, is set to 0.0 if the run is a new one */

  htrybe = 0.5;           /* initial time step of integration */
  
  max_time = 50000*365.25; //8000000.0*365.25; /* lenght of the numerical
		    //integration in years                 */

  GRADUAL = 0;                 /* mode of writing out the data                                 */
                               /* USAGE: GRADUAL = 0, data are written periodically, where the
                                  period is specified by the user
                                         GRADUAL = 1, between 0 - 1000 yrs    in every 2 years
                                                           1000 - 10000 yrs   in every 10 years
                                                          10000 - 100000 yrs  in every 100 years 
                                                         100000 - and above   in every 1000 years */

  saving_period = 10.0; //*365.25;       /* period of writing out the data into the output file in years */
  creation_period = 5000.0*365.25;        /* NOT IN USE !!! */
  final_creation_time = max_time;

//  RESCMIN = ;               /* inner boundary of the computation domain */
//  RESCMAX = 30.0;             /* outer boundary of the computation domain */
  
  GASDRAG = 0;                     /* 1 for taking into account gas drag, 0 for no drag !!!NOT IN USE NOW!!!*/

  MIGRATION_TYPE1 = 1;             /* MIGRATION_TYPE1 = 1: type I migration migration, timescales 
                                      are calculated in N_Body_Force()                                      */
  
  ADIABATIC = 1; /* ADIABATIC = 1: type I migration in adiabatic disk */
                 /* ADIABATIC = 0: type I migration in isotherm disk  
                    use together with "MIGRATION_TYPE1 = 1" option    */

  CONTINUATION = 0;    /* if = 0, input file created by inic.c 
                          if = 1, continuation of an integration using an input 
		          created from the outputfile of a previous run           */

  GAS_TORQUE = 1; /*  If 0, don;t read in gas profiles, if 1, gas disk exerts a torque on the planets  */
  UPDATE_OPACITY=1;

  BINFORM=0;
  int holdp=0;

if (GAS_TORQUE ==1){

//printf("Begining\n");
 itgas=0;
 init_time(fin_time, &tgas);
  srand(57);
  init_radii(tgas,&Radii,&Sigma,&Alpha,&Temp,&Beta); // gas disk properties. note radii in
						     // AU now, sigma in Msun/AU2
printf("Radii initiated\n");

int hind;
//read in aspect ratios
for (hind=0;hind<NRAD;hind++) {
           fscanf(fin2, "%lg", &ar);
	   aspR[hind]=ar;
} 

 
//read_density_data(NRAD,itgas,tgas,&Sigma_1,&Sigma_2,&Alpha_1,&Alpha_2,&Temp_1,&Temp_2,&Beta_1,&Beta_2);  // again, sigma in Msun/AU2
 R_inner=Radii[0];
 R_outer=Radii[NRAD-1];
 RESCMIN = R_inner-0.00001;               /* inner boundary of the computation domain */
 RESCMAX = R_outer+0.00001;
 printf("RESCMAX=%f",RESCMAX);
 D_r=(R_outer-R_inner)/(1.0*NRAD);
 read_op_files(&aindexs,&asevens,&asix5s,&asixs,&afive5s,&afives,&afour5s,&afours,&athree5s,&athrees,&atwo5s,&atwos,&aone5s,&aones,&azero5s,&azeros,&apzero5s,&apones);
  } else if (GAS_TORQUE ==0){
	for (I=0;I<NBOD*3;I++) {
 	  Sigma_1[I]=0.0;
	  Sigma_2[I]=0.0;
	  Alpha_1[I]=0.0;
	  Alpha_2[I]=0.0;
	  Temp_1[I]=0.0;
	  Temp_2[I]=0.0;
	  Beta_1[I]=0.0;
	  Beta_1[I]=0.0;
	  Sigma[I]=0.0;
	  Alpha[I]=0.0;
	  Beta[I]=0.0;
	  Temp[I]=0.0;
	}	
  }
///intiating nbody code *JIA9/21*
for (I=0;I<NBOD;I++) {
 index[I]=I;
}
 set_rand_log(); 
Generate_Lambda(Radii,aspR,BODY0);

  if (CONTINUATION == 1) {
  
      printf("This run is a *continuation* of a previous one!\n");
      
      for (I=0;I<NBOD;I++) {
            /* reading the masses and the orbital elements from the given file */
	
	
	    fscanf(fin1,"%d %d %lg %lg %lg %lg %lg %lg %lg %lg",&int_dummy,&int_index,&m,&dummy,&ax,&ex,&inc,&MAn,&omP,&OmN);
	    printf("%lg %lg %lg %lg %lg %lg %lg\n",m,ax,ex,inc,MAn,omP,OmN);
	    // not yet converting migrator masses and SMA to cgs units  JMB 4/12/17
	   

	    EPOCH =dummy*365.25;

	    /* calculation of the star-centered Cartesian coordinates and velocities 
	       from the star-centered orbital elements*/
	    coordinate(BODY0,m,ax,ex,inc,omP,OmN,MAn,&xx,&yy,&zz,&vxx,&vyy,&vzz);
	
	    ybe[6*I]   = xx;
	    ybe[6*I+1] = yy;
	    ybe[6*I+2] = zz;
	    ybe[6*I+3] = vxx;
	    ybe[6*I+4] = vyy;
	    ybe[6*I+5] = vzz;

            yki_after[6*I]   = xx;
            yki_after[6*I+1] = yy;
            yki_after[6*I+2] = zz;
            yki_after[6*I+3] = vxx;
            yki_after[6*I+4] = vyy;
            yki_after[6*I+5] = vzz;
	
	    index[I]=int_index;
	    BODY[I] = m;
	    BODY_BAFTER[I] = m;
     } 
     } else if (CONTINUATION == 0) {

             /* New run from an input file */
	     EPOCH = 0.0;
             for (I=0;I<NBOD;I++) {
	          fscanf(fin1,"%lg %lg %lg %lg %lg %lg %lg",&m,&ax,&ex,&inc,&MAn,&omP,&OmN);
                  //inc = my_rand()*Inc_max;
		  if (i=1)
		  {if (inc < PI); //Q?JIA
		    inc = (114.59*PI) - inc;
}
		       
		  printf("%lg %lg %lg %lg %lg %lg %lg\n",m,ax,ex,inc,MAn,omP,OmN);
		  coordinate(BODY0,m,ax,ex,inc,omP,OmN,MAn,&xx,&yy,&zz,&vxx,&vyy,&vzz);
		    // not yet converting migrator masses and SMA to cgs units  JMB 4/12/17

	          ybe[6*I]   = xx;
	          ybe[6*I+1] = yy;
	          ybe[6*I+2] = zz;
	          ybe[6*I+3] = vxx;
	          ybe[6*I+4] = vyy;
	          ybe[6*I+5] = vzz;

                  yki_after[6*I]   = xx;
                  yki_after[6*I+1] = yy;
                  yki_after[6*I+2] = zz;
                  yki_after[6*I+3] = vxx;
                  yki_after[6*I+4] = vyy;
                  yki_after[6*I+5] = vzz;
	
		  index[I] = I;
	          BODY[I] = m;  
	          BODY_BAFTER[I] = m;  
	    }

	
  }
      

  body_creation_times(SEED2,creation_period,final_creation_time,Times);

  nr_of_times = 0;
  
  xbe = EPOCH;
   do {
         /* collision index, 0 for no collision */
         COLLISION = 0;
	 /* escaping index, 0 for no escape */
	 ESCAPE = 0;
	 /* detection of escapes, but only the escaping body 
	    with the lowest index is detected; 
	    iesc: index of the escaping body */
	 detect_escape(xbe,ybe,NBOD,&iesc,&ESCAPE,index);
	 
	 /* here is eliminated the escaping body  */
	 if (ESCAPE == 1) {
	    
	     NBOD = NBOD - 1;
	     nv   = 6*NBOD;
	     
	     for (i=iesc; i<NBOD; i++) {
	     
	         i1 = i + 1;
	  
	         ybe[6*i+0] = ybe[6*i1+0];
	         ybe[6*i+1] = ybe[6*i1+1];
	         ybe[6*i+2] = ybe[6*i1+2];
	         ybe[6*i+3] = ybe[6*i1+3];
	         ybe[6*i+4] = ybe[6*i1+4];
	         ybe[6*i+5] = ybe[6*i1+5];
	  
	         BODY[i] = BODY[i1];
		 index[i] = index[i1];
	     }   
	 }
       if (GAS_TORQUE ==1){

	//right now this is commented out because Im only reading in one file
	 /*if (itgas < (NTS-1)) {
		 if (xbe > tgas[itgas+1]*365.25) {
			//itgas=itgas+1;
			//read_density_data(NRAD,itgas,tgas,&Sigma_1,&Sigma_2,&Alpha_1,&Alpha_2,&Temp_1,&Temp_2,&Beta_1,&Beta_2);
	 	}
	 
           if (xbe/365.25 > tgas[0]) {
	     update_density(xbe,tgas[itgas],tgas[itgas+1],NRAD,Radii,Sigma_1,Sigma_2,Alpha_1,Alpha_2,Temp_1,Temp_2,Beta_1,Beta_2,&Sigma,&Alpha,&Temp,&Beta);
           }
	 }*/
	}
	 /*performs one step of numerical integration*/
       Lambda_Update(xbe,Radii,aspR,BODY0); // original SS units here
      get_sma(BODY0,NBOD, BODY,Radii,ybe,&sma);

      if (htemp2print<oldhtemp){
        oldhtemp=htemp2print;
      }

      Integrator(&holdp,saving_period,&htemp2print,BODY0,NBOD,BODY,Radii,Sigma,aspR,Alpha,Temp,Beta,ybe,xbe,htrybe,sma,yki,&xki,&hkovki,&F_stoc_out,&aindexs,&asevens,&asix5s,&asixs,&afive5s,&afives,&afour5s,&afours,&athree5s,&athrees,&atwo5s,&atwos,&aone5s,&aones,&azero5s,&azeros,&apzero5s,&apones);


      //commenting out this verion of the integrator that accounts for binary protocal
      //Integrator(&holdp,saving_period,&htemp2print,BODY0,NBOD,BODY,BODY_BAFTER,Radii,Sigma,Alpha,Temp,Beta,ybe,yki_after,xbe,htrybe,sma,yki,&xki,&hkovki,&F_stoc_out,&aindexs,&asevens,&asix5s,&asixs,&afive5s,&afives,&afour5s,&afours,&athree5s,&athrees,&atwo5s,&atwos,&aone5s,&aones,&azero5s,&azeros,&apzero5s,&apones);

      /*if (fmod(xbe,1000)<=1){
          print_positions(NBOD,yki,index,&F_stoc_out,xbe);
	}*/							//x,y,z now printed in output.dat

	
	//commenting out the binary methods I developed for now
	/* if (BINFORM == 0){
           tbmove=xbe;
	   tbmerge=xbe;

          for (i=0; i<NBOD; i++) {
              BODY_BAFTER[i]=BODY[i];
          }

          for (i=0; i<nv; i++)
              yki_bafter[i]=yki[i];
	 }

         if (BINFORM == 1) {

             binary(xki,NBOD,BODY,yki,index,yki_bafter,BODY_BAFTER,ybe); //both indexes of the binary gets COM vector and sum of masses
	     printf("yki_bafter-yki= %g\n",yki_bafter[0]-yki[0]);
	     if ((xbe-tbmove)>thard){
	       printf("It's been thard\n");
	       bpos(yki,yki_bafter,adj_pos,nv);
	       tbmove=xbe;
               for (i=0; i<nv; i++){
                  yki[i] = adj_pos[i];	     
	       }
	    }
	    if ((xbe-tbmerge)>(40*thard)){
		COLLISION=1;
		tbmove=xbe;
		tbmerge=xbe;
	    }
         } */

	 /*re-initialization of the gobal variabes if collision occured*/
	 if (COLLISION == 1) {
	     BINFORM=0; 
	     collision(xki,NBOD,BODY,yki,index,&NEW_NBOD,yki_after,BODY_AFTER,&index_new,ybe,Radii,aspR,Alpha,Sigma,Temp,Beta,BODY0);
	   
	     NBOD = NEW_NBOD;
	     nv   = 6*NEW_NBOD;
	     
	     for (i=0; i<NBOD; i++) {
	          BODY[i] = BODY_AFTER[i];
		  index[i] = index_new[i];
	     }

	     for (i=0; i<nv; i++)   
                  yki[i] = yki_after[i];
	   	  
	 }
	 
         /* creation of a body in a given time */
/*         if (xbe < Times[nr_of_times]  && Times[nr_of_times] <= xki) {
             create_a_body(NBOD,Sma_min,Inc_max,Ecc_max,BODY,yki,&xx,&yy,&zz,&vxx,&vyy,&vzz,&m_cr);
             //NBOD = NBOD+1;
         
             yki[nv] = xx;
             yki[nv+1] = yy;
             yki[nv+2] = zz;
             yki[nv+3] = vxx;
             yki[nv+4] = vyy;
             yki[nv+5] = vzz;
             BODY[NBOD] = m_cr;
             index[NBOD] = NStart+nr_of_times;
	     nr_of_times++;
             printf("Number of bodies: %d\n",NBOD);
	     printf("Next creation time: %lg Number: %d\n",Times[nr_of_times]/365.25,nr_of_times);
	     fprintf(fout3,"Next creation time: %lg Number: %d\n",Times[nr_of_times]/365.25,nr_of_times);
	     fflush(fout3);
             //getchar();
             NBOD ++;
             nv = 6*NBOD;
             
         }
*/
	 
	/*output of the data*/
	//now prints every x timesteps
	//printf("PTIME %d\n",PRINT_TIME);
	if(PRINT_TIME==1){
		print_elements(&oldhtemp,xbe,xki,NBOD,BODY0,BODY,yki,GRADUAL,index); //worrried something is wrong here b/c changing its resolution changes the collision times
	}

	UPDATE_OPACITY=1;

	/*re-initialization of the input data for the next integration step*/
	 
	 xbe = xki;
	 //         printf("xbe %lg  max_time  %g \n",xbe,max_time);
	 htrybe = hkovki;
	 
         for (i=0; i<nv; i++)
             ybe[i] = yki[i];
	 //	 printf(" max_time %g \n",max_time);
   } while (xbe < max_time);
      
 fclose(fin1);
 fclose(fin2);
 fclose(fin_time);
 fclose(ptr_myfile);
 //fclose(surf);
 //fclose(foutb);
 //fclose(foutm1);
 //fclose(foutm2);
 //fclose(foute);
 fclose(fout_sd);
 fclose(fout3);
 //fclose(fout4);
 fclose(fout5);
 fclose(abkap); 
 fclose(hconv); 
}
