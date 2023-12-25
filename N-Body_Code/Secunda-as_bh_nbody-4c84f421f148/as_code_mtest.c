#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#define glimax          11
#define glnmax          1200
#define glncol          7
#define MASSJUP         0.0009548
#define PI		3.141592653589793
#define TWOPI           6.283185307179586
#define	TORAD		PI/180.0
#define TODEG           180.0/PI
#define G               0.01720209895  // the Gaussian grav. constant
#define G_CGS		6.674e-8
#define SUN2GR		1.989e33
#define AU2CM		1.496e13
#define SDCONV		1.12521e-7  // convert Msun/AU^2 to gram/cm^2
#define SPEEDOFLIGHT    3e10  // cm per second
#define UNO		1.0;
#define TWO             2.0;
#define UNO2		0.5;
#define ERROR		1.0e-13;
#define CERO		0.0;
#define MAXNBOD         200              // maximum number of bodies except of the central one
#define BETA            0.3;		//what fraction of a hill radius a merger takes place at
#define NPLMAX          10
#define FLT_MAX         3.4028E+38F      // maximum float value
#define NTS             1000000         //number of Time Steps in time.dat from gas profiles
#define NRAD		5049           //radial resolution in gas simulation: 5049 for SG 1699 for TQM
#define NStart		10           //number of black holes in elements.dat
#define MODES		200          //num of harmonics in turb potential
#define gamma		4.25e-4         //strength of turb potential

char * gasdirect  ="/home/asecunda/as_bh_nbody/SGtimefiles/";//location of gas profiles


// MAXNBOD = glnmax/6


typedef double vektor[glnmax];

 long nv,  NBOD, NEW_NBOD, IMIN,IBMIN, JMIN, JBMIN;
 int CLOSE,nclose[MAXNBOD],COLLISION,BODY_CREATE, PRINT_TIME, ESCAPE, GASDRAG, MIGRATION_TYPE1, ADIABATIC, NPL, NPHI, LamIndex[MODES], UPDATE_OPACITY;
 
double xbe, tbmove,tbmerge, thard, htrybe, hkovki, max_time, xki, BODY0, BODY[MAXNBOD], BODY_AFTER[MAXNBOD],  RHO, Lambda[6*MODES], LamTime[MODES],ProbLog[257],Dyn_Disp[MAXNBOD],Disp_Array[MAXNBOD],I_BVECTOR[6],J_BVECTOR[6],OPACITY[MAXNBOD];
 double RESCMIN, RESCMAX, SDEXP, SIGMA_0, DISP_BEGIN, DISP_DURATION, EPOCH, ADINDEX, EXPSD, DCRITICAL,R_inner,R_outer,D_r;
 
 
 vektor ybe, yki, yki_after, yki_bafter,adj_pos, yki_after_esc;
 
 
 double glx[glimax];
 double gld[glnmax][glncol];
 
FILE *ptr_myfile, *surf, *fin1,*fin2,*fin3,*fin4,*fin_time,*fout3;


/*structure for printing into bout.bin the main output of the code*/
struct rec
        {
        float i,j,mass,time,ax,ex,inc,MAn,omP,OmN,x,y,z,vx,vy,vz;        
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

//calculates the schwarzchild radius of a body, not currently in use
double radius(i,body)
int i; double *body;
{
    double schwarz;
    schwarz = 2.0*body[i]*SUN2GR*G_CGS/(SPEEDOFLIGHT*SPEEDOFLIGHT);
    schwarz = schwarz/AU2CM; // consistent with all AU units
    return schwarz;
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

//calculates the hardening timescale. not currently in use.
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

//calculates eccentricity and inclination
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

//converts from rectangular coordinates to orbit coordinates
void element(tomeg0,tomeg1,x,y,z,vx,vy,vz,ax,exc,ink,omegaP,OmegaN,kozan)
double tomeg0,tomeg1,x,y,z,vx,vy,vz;      // coordinates and velocities of the orbiting body
double *ax;     //semi-major axis
double *exc;    // eccentricity
double *ink;	// inclination          
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
  c = sqrt(cx*cx + cy*cy + cz*cz);		      //absolut value of the angular momentum
  
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

//converts from orbit coordinates to rectangular coordinates
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
	
      INC = INC*TORAD;   // TORAD = to radians
      MEAN = MEAN*TORAD;
      W = W*TORAD;
      OMN = OMN*TORAD;
      
      EXAN = Solve_Kep_eq(ECC,MEAN);
	      
//    DIRECTOR COSINES    
    
      CU  = cos(EXAN);
      SU  = sin(EXAN);
      CUE = CU - ECC;
      PX  = cos(W)*cos(OMN) - sin(W)*sin(OMN)*cos(INC);
      QX  =-sin(W)*cos(OMN) - cos(W)*sin(OMN)*cos(INC);
      PY  = cos(W)*sin(OMN) + sin(W)*cos(OMN)*cos(INC);
      QY  =-sin(W)*sin(OMN) + cos(W)*cos(OMN)*cos(INC);
      PZ  = sin(W)*sin(INC);
      QZ  = cos(W)*sin(INC);
      SMA1 = SMA*sqrt(1.0-ECC*ECC);
      RADI = (1.0-ECC*CU)*SMA;
      UDOT = SMA*VELM/RADI;

//     HELIOCENTRIC COORDINATES AND VELOCITIES X AND XDOT [AU] & [AU/DAY]

      *x = SMA*PX*CUE+SMA1*QX*SU;           // X
      *y = SMA*PY*CUE+SMA1*QY*SU;           // Y
      *z = SMA*PZ*CUE+SMA1*QZ*SU;           // Z
      *vx = UDOT*(-SMA*PX*SU+SMA1*QX*CU);   // velocity components
      *vy = UDOT*(-SMA*PY*SU+SMA1*QY*CU);   
      *vz = UDOT*(-SMA*PZ*SU+SMA1*QZ*CU);   
      
}  

//finds the local value of surface density at r
double get_sigma(r,radii,surfd)
double r, *radii,*surfd;
{
  int index;
  double d_r1,r_remainder,s_dens;

  if (r>R_outer){
     r=R_outer -0.01;
   }
  
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

void get_local_gas_properties(r,radii,aspR,taup,expal,surfd,temp,beta,s_dens,T,aspR_local,taup_local,alpha_local,beta_local)
double r, *radii,*expal,*surfd,*temp,*beta,*aspR,*taup;
double *s_dens,*T,*alpha_local,*beta_local,*aspR_local,*taup_local;
{
  int index;
  double d_r1,r_remainder;

  if (r>R_outer){
     r=R_outer -0.01;
   }

  d_r1=1.0/D_r;
  index=((r-R_inner)*d_r1); // linear interp gives a misleading
			    // (wrong) index
  r_remainder=(r-radii[index])*(r-radii[index])/NRAD;


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
  *taup_local = taup[index]+r_remainder*(taup[index+1]-taup[index]);
}

//derives speed of sound from local temp. This is used to find scale height of disk
double get_speed_of_sound(r,radii,aspR,Mcenter)
double r,*radii,*aspR,Mcenter;
{
  
  int index;
  double d_r1,r_remainder,H,k,mH,cs,adiabat;

  if (r>R_outer){
     r=R_outer -0.01;
   }
  
  d_r1=1.0/D_r;
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
  H = aspR[index]+r_remainder*(aspR[index+1]-aspR[index]);
  cs= H*sqrt(G*G*Mcenter/r);
  return cs;  
}

//calculates the migration torque, eccentricity, and inclination dampening on a black hole in the disk 
void get_migration_damping(body0,bodyi,x,y,z,vx,vy,vz,r,sig,T,aspR,alpha,beta,Torque,Damp,Theta,closei)
double body0,bodyi,x,y,z,vx,vy,vz,r,sig,T,aspR,alpha,beta,*Torque,*Damp,Theta;
int closei;
{
    double ecc, inc, entropy_exp, h, rte,rti,d_torque, torque_rn2,rv,torque_norm;

  get_e_i(body0,bodyi,r,x,y,z,vx,vy,vz,&ecc,&inc);
  h=aspR;


  Theta=((Theta/SUN2GR)*(pow(AU2CM,2)))*(60*60*24); //converting this from cgs to AU, Msun

  rte=h*(r*body0*sqrt(r*body0))*(h*h*h-0.14*ecc*ecc*h+0.06*ecc*ecc*ecc+0.18*ecc*inc*inc)/(0.78*bodyi*sig*G);
  rti=h*(body0*sqrt(body0))*(h*h*h-0.3*inc*inc*h+0.24*inc*inc*inc+0.14*ecc*ecc*inc)/(0.544*bodyi*sig*G*sqrt(r));

  entropy_exp = beta - alpha*(ADINDEX-1.0);
  d_torque = ((1.0/ADINDEX)*(-0.85 - alpha - 1.7*beta + 7.9*entropy_exp/ADINDEX)*Theta*Theta + (-0.85-alpha-0.9*beta))/((Theta+1)*(Theta+1));

  torque_rn2=d_torque*G*G*bodyi*bodyi*sig/(h*h*body0*r);
  rv=x*vx+y*vy+z*vz; //r dot v

  if (closei!=1){
    /* Torque[i] below is actually a Force in x and y components*/
    Torque[0]= -1.0*torque_rn2*y; 
    Torque[1]= 1.0*torque_rn2*x;
  }

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
double Vol_Density(a,h,x,y,z,sig_0)
double a, h, x, y, z, sig_0;
{
   double r;
  
   r=sqrt(x*x + y*y);
   return sig_0*pow(r,a-1.0)/(h*sqrt(2.0*PI)*1.49e13)*exp(z*z/(2.0*h*h*r*r));  
}


// procedure not used in this version
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

//finds gradient of turbulent potential at black hole's location. Lambda is a global variable. 
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
//finds minimum distance between two point particles deemed to be colliding
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
double Impact_energy(m1,m2,v,r)
double m1,m2,*v,r;
{
  double KE, BE,mu;
  mu=(m1*m2)/(m1+m2);
  KE=0.5*mu*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  BE= 0.5*G*G*(m1*m2)/r;
  return (KE/BE);

}


//returns semi-major axis of body

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
      }
}



// procedure not used
double ALPHA(double radius)
{
 double eta = 0.0019*sqrt(radius);
 return sqrt(1.0 - 2.0*eta);
}

// not used 
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

//calculates the local opacity of the disk
double Opacity(r,aspR_local,taup_local,sig,T,BODY0)
     double r,sig,T,BODY0,aspR_local,taup_local;
{
  double omega, kap,taup,taueff,theta, cv;


  sig = sig/SDCONV;  // now in g/cm2
  taup=taup_local; //0.5*kap*sig;
  taueff=0.375*taup+sqrt(3.0)*0.25+0.25/taup;
  cv = 12.5e7; // JMB 4/12/17 updated to monatomic gas 8.314e7*3/2; 
  r = r*AU2CM;
  omega = sqrt(G_CGS*BODY0*SUN2GR/(pow(r,3.0)));// cgs units!

  theta=cv*sig*omega*taueff/(4.0*PI*5.67e-5 *T*T*T);
return theta;
}

//this calculates the distances between bodies.
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
    }
  }
}

// set up an array of distances between all bodies and all other
// bodies and their components
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

/* calculation of the mutual gravitational forces in the SMBH centered 
   coordinate system */
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
      cent_acc[0]+=gsum*xi;    // summing the force on the central body
      cent_acc[1]+=gsum*yi;
      cent_acc[2]+=gsum*zi;

      force[i1] = y[i4];  
      force[i2] = y[i5]; 
      force[i3] = y[i6];

      gcent=grav*rin3*cbody;  // G*M_BH/r3

      force[i4] = - gcent*xi;  // force on orbiting body
      force[i5] = - gcent*yi; 
      force[i6] = - gcent*zi;
      nclose[i]=0;
 }
  for (i=0;i<nbody;i++) {
    force[i*6+3]+=-cent_acc[0];  // cent_acc is a force.
    force[i*6+4]+=-cent_acc[1];
    force[i*6+5]+=-cent_acc[2];
  }
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
      if(rij2<pow((1.2*rhill(body[i],body[j],cbody)),2)){
	    nclose[i]=1;
	    nclose[j]=1;
      }  
    }
  }
	//these are used to check for binary formation/merger
      *mindist = sqrt(dmin2);
      *i_dmin  = imin;
      *j_dmin  = jmin;
}

   

/* adding additional forces: originating from migration, and drag                */
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
    dmin =10000.0;

    /* Gravitational force between all the bodies */
    Gravitational_Force(Ntot,Body0,Body,Y,Force,&dmin,&imin,&jmin);


    /* Migration of the planets: --> IFLAGMIGR = 1 
                    no migration --> IFLAGMIGR = 0 
       using the approach of Cresswell & Nelson (2008)                  */

    if (MIGRATION_TYPE1 == 1 ) {
	
	
	
        for (i=0;i<Ntot;i++) {

             i1 = 6*i;
             i2 = i1+1;
             i3 = i2+1;
             i4 = i3+1;
             i5 = i4+1;
             i6 = i5+1;
             ri = sqrt(Y[i1]*Y[i1] + Y[i2]*Y[i2] + Y[i3]*Y[i3]);

             Force[i4] = Force[i4] +  F_stoc[i];
             Force[i5] = Force[i5]+ F_stoc[i+Ntot];
	     Force[i6] = Force[i6]+ F_stoc[i+2*Ntot];


        }
    }
    /* migration is added */


    /*saves a position and velocity vector for the two nearest bodies*/

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
      bcap=rhill(Body[imin],Body[jmin],Body0)*BETA; //here you can reset the merger 
						   //condition by any factor of rhill

      //lets you print all the timesteps (instead of skipping however many)
      //once two black holes are within 1 rhill and before they merge
      if (dmin<(1.2*rhill(Body[imin],Body[jmin],Body0))){
        CLOSE=1;
      }
      else{
	CLOSE=0;
      }

      if (dmin<bcap){
	 xrel[0]=I_BVECTOR[0]-J_BVECTOR[0];
	 xrel[1]=I_BVECTOR[1]-J_BVECTOR[1];
         xrel[2]=I_BVECTOR[2]-J_BVECTOR[2];
         vrel[0]=I_BVECTOR[3]-J_BVECTOR[3];
         vrel[1]=I_BVECTOR[4]-J_BVECTOR[4];
         vrel[2]=I_BVECTOR[5]-J_BVECTOR[5];
         rrel=sqrt(xrel[0]*xrel[0]+xrel[1]*xrel[1]+xrel[2]*xrel[2]);
	 ratio=Impact_energy(Body[imin],Body[jmin],&vrel,rrel); //checks that kinetic<binding 
	 if (ratio<1){
         	COLLISION=1;
        	IMIN = int_min(imin,jmin);
        	JMIN = int_max(imin,jmin);
	 }
      }

}


/* calculation of the new velocity of the body, elimination of the colliding body
   collisions: always accretional completely non-elastic, the lowest index body will be
   the new one the masses are added
                                            */
void collision(time,nbod,body,y,index,new_nbod,y_new,body_new,index_new,y_before,radii,aspR,taup,alpha,sigma,temp,beta,M0)
double time; int nbod; double *body, *y, *y_before; int *new_nbod; double *y_new; double *body_new; int *index; int *index_new,M0;
double *radii,*alpha,*sigma,*aspR,*taup,*temp,*beta;
{
 int icol, jcol, i, i1;
 double tmp1, tmp2, x[3],v[3],rmin,r,v1c[3],v2c[3],V1,V2,ratio,thard,new_a;
 double s_dens,T,aspR_local,taup_local,alpha_local,beta_local;

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
      
      printf("collision between bodies %ld %ld at time %lg days\n",IMIN,JMIN,time);

      fprintf(fout3,"%ld %ld %lg %lg %lg ",index[icol],index[jcol],body[icol],body[jcol],xki);
      fprintf(fout3,"%lg %lg %lg %lg %lg %lg\n",r,rmin,v[0],v[1],v[2],ratio);
      
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
      body_new[icol] = (body[icol] + body[jcol]);//*0.95;
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
      

      fflush(fout3);
}



/* here are the numerical integrator procedures*/

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


void BSmethod(holdp,period,central,ntot,bodies,radii,sigma,alpha,temp,beta,y, F_stoc,dydx, x, htry, eps, loc_yki, loc_hki, loc_hkovki)
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

//saving period
if (hcount>period){
  PRINT_TIME=1;
  hcount=0;
}
else if (CLOSE==1){
  PRINT_TIME=1;
  hcount=0;
}
else{
  PRINT_TIME=0;
  hcount++;
}
*holdp=hcount;
}

#undef imax
#undef nuse
#undef one
#undef shrink
#undef grow


void Integrator(holdp,period,cbody,ntot,bodies,radii,sigma,aspR,taup,alpha,temp,beta,y,x,h,sma,yo,xo,hn,F_out)
     double  cbody;  /*mass of the central body*/
     int     ntot;   /*number of bodies except the central one*/
     int period;
     double *bodies; /*array of the masses*/
     double *radii;
     double *sigma;
     double *aspR;
     double *taup;
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
     int *holdp;
{

 vektor dydx;

 double hki, tol,lhn, ri, s_dens,FS[3],T,aspR_local,taup_local,alpha_local,beta_local,Torque[2],Damp[3];
 int cv,i,i1,i2,i3,i4,i5,i6,j,i_int,k;

        for (i=0;i<ntot;i++){
             ri=sqrt(y[6*i]*y[6*i]+y[6*i+1]*y[6*i+1]+y[6*i+2]*y[6*i+2]);
             i1=6*i; i2=i1+1; i3=i2+1; i4= i3+1; i5=i4+1; i6=i5+1;
             /*STOCHASTIC FORCES!!!!!*/

             get_local_gas_properties(ri,radii,aspR,taup,alpha,sigma,temp,beta, &s_dens,&T,&aspR_local,&taup_local,&alpha_local,&beta_local);
             if (UPDATE_OPACITY==1){
                 OPACITY[i]=Opacity(ri,aspR_local,taup_local,s_dens,T,BODY0);
                }

             get_migration_damping(cbody,bodies[i],y[i1],y[i2],y[i3],y[i4],y[i5],y[i6],ri,s_dens,T,aspR_local,alpha_local,beta_local,Torque,Damp,OPACITY[i],nclose[i]);
             get_Fstoc(y[6*i],y[6*i+1],x,ri,cbody,s_dens,FS);
             F_out[i]=FS[0]+Torque[0]+Damp[0];
             F_out[i+ntot]=FS[1]+Torque[1]+Damp[1];
             F_out[i+2*ntot]=FS[2]+Damp[2];
        }
UPDATE_OPACITY=0;


 tol = 1.0e-13;
 Set_Disp_Array(ntot,y);
 Set_Dyn_Disp(ntot);
 N_Body_Force(cbody,ntot,bodies,radii,x,sigma,alpha,temp,beta,y,F_out,dydx);
 BSmethod(holdp,period,cbody, ntot, bodies, radii, sigma, alpha, temp,beta, y, F_out, dydx, x, h, tol, yo, &hki, &lhn);

 *xo = x + hki;
 *hn = lhn;
}




/* prints out the results - orbital elements at given period, of gradually */
//now prints every x timesteps
void print_elements(timeold,time,nbod,cbody,bodies,y_vekt,gradual,index)
double timeold; double time; int nbod; double cbody; 
double *bodies; double *y_vekt; int *index;
int gradual;
{
  int i,j;
  double x,y,z,vx,vy,vz,mass,period;
  double ax,ex,inc,MAn,omP,OmN,Ex,Ey,vmag,rmag,mu,r,ecc2,inc2;


  struct rec my_record;  
  if (!ptr_myfile)
  {
     printf("Unable to open file!");
  }
 

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
	my_record.MAn=MAn;
	my_record.omP=omP;
	my_record.OmN=OmN;
	my_record.x=x;
	my_record.y=y;
	my_record.z=z;
	my_record.vx=vx;
	my_record.vy=vy;
	my_record.vz=vz;

	//write output binary file
	fwrite(&my_record, sizeof(struct rec), 1, ptr_myfile);

	   fflush(ptr_myfile);
 } 
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
	  break;
       }
   }   
           
}

//reads in the file listing all the filenames of the gas profiles--each filename is the year of the profile
//this method is for an evolving disk. it is overwritten in this version because the disk does not evolve
void init_time(f_time, tgass)
FILE *f_time;
int *tgass;
{
   int ttt, I;
	for (I=0; I < NTS; I++) {
	 	fscanf(f_time, "%d",&ttt);
		tgass[I]=ttt;
       }
	printf("%d%d",tgass[0],tgass[100]);
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
}

//reads in next gas profile.
//this method is now hard coded to read in a fil called "0.dat" because the disk doesn't evolve
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
   if (itgas*1.0 > (NTS-1)*1.0-0.1) { // this is the last time step
				      // JMB
	for (i = 0; i < nr; i++){
		sigma_2[i]=sigma_1[i];
		alpha_2[i]=alpha_1[i];
		beta_2[i]=beta_1[i];
		temp_2[i]=temp_1[i];
	}
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

}

/* creates an array of times seperated by tau yrs from time tau to the end of the run (t_final) used
to determine when to drop in bodies from the outer disk*/
 void body_creation_times(seed,tau,t_final,times)
unsigned int seed;
double tau, t_final, *times;
{
  int i;
  double t_old, t_new;

  for (i = 0; i< 1500; i++)
        times[i] = 0.0;

  i=0;
  t_old = 0.0;
  srand(seed);
  do {
      
      t_new = t_old +tau;
      
      times[i] = t_new;
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
  
  factor = 10.0;

  /* get mass from file */
  fscanf(fin4, "%lg", &mnew);
 
_L88:; 
  j++;
  sma = sma_min;
  ecc =  ecc_max* my_rand();
  in  =  inc_max* my_rand();
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
  
  double ar,tau,xx, yy, zz, vxx, vyy, vzz, m, ax, ex, inc, omP, OmN, MAn, Inc_max, creation_period,
         final_creation_time, dummy,m_cr, Sma_min, Ecc_max, axis, ecce, incl, ompe, omno, mean, r, surfdens;
  double Radii[NRAD], Alpha[NRAD], Alpha_1[NRAD], Alpha_2[NRAD],Sigma[NRAD], Sigma_1[NRAD], Sigma_2[NRAD], Temp[NRAD]; 
  double Temp_1[NRAD], Temp_2[NRAD],Beta[NRAD],Beta_1[NRAD], Beta_2[NRAD],Times[1500],aspR[NRAD],taup[NRAD], 
    F_stoc_out[MAXNBOD*3],x[3],d2,dmin,dcrit; 
  int    I, i, i1, iesc, int_dummy, CONTINUATION, GRADUAL, t, itgas, int_index,GAS_TORQUE,icol,jcol,saving_period;
  int    tgas[NTS], index[MAXNBOD], index_new[MAXNBOD];
  unsigned int SEED1, SEED2, SEED3;
  long nr_of_times;
  char file1 [ FILENAME_MAX ],hors[FILENAME_MAX],taus[FILENAME_MAX], timefile [ FILENAME_MAX ];
  vektor sma;
  /*********************************************************************************************************/
  /* file declarations                                                                                     */
  /*********************************************************************************************************/
  fin1 = fopen("/scratch/gpfs/asecunda/ligo/test/ex_elements.dat","r"); //this file contains the starting positions and mass of the BHs
  sprintf(hors,"%s%s",gasdirect,"hor.dat");
  fin2 = fopen(hors,"r"); //this file has the aspect ratios for all disk radii
  sprintf(taus,"%s%s",gasdirect,"tau.dat");
  fin3 = fopen(taus,"r"); //this file has the optical depths for all disk radii
  sprintf(timefile,"%s%s",gasdirect,"time.dat");
  ptr_myfile=fopen("/scratch/gpfs/asecunda/ligo/test/bout.bin","wb"); //this is the main output file
  fin_time = fopen(timefile,"r"); //this file has the temperature, surface density and their gradients for all disk radii
  fout3= fopen("/scratch/gpfs/asecunda/ligo/test/analysis.txt","w"); //this file has info on all of the mergers
  /*********************************************************************************************************/
  
  NBOD = NStart;         /* number of bodies */
  nv   = 6*NBOD;
  SEED1 = 23;        /* random seed for the random number generator to calculate inclinations */
  SEED2 = 26;        /* random seed for the random number generator to calculate the body creation times */
  SEED3 = 54;        /*                                                                   positions      */
  saving_period = 500;
  Sma_min = 1000.0;    //where created bodies will appear
  Ecc_max = 0.07;   //maximum value of eccentricities for created bodies
  Inc_max = 0.05;    /* maximum value of the orbital inclinations for created bodies*/
  
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
  
  max_time =  1e6*365.25; /* lenght of the numerical integration in days                 */

  creation_period = 100000.0*365.25;     /*how often black holes from the outer disk will migrate in in days*/
  final_creation_time = max_time;

  
  GASDRAG = 0;                     /* 1 for taking into account gas drag, 0 for no drag !!!NOT IN USE NOW!!!*/

  MIGRATION_TYPE1 = 1;             /* MIGRATION_TYPE1 = 1: type I migration migration, timescales 
                                      are calculated in N_Body_Force()                                      */
  
  ADIABATIC = 1; /* ADIABATIC = 1: type I migration in adiabatic disk */
                 /* ADIABATIC = 0: type I migration in isotherm disk  
                    use together with "MIGRATION_TYPE1 = 1" option    */

  CONTINUATION = 0;    /* if = 0, input file created by inic.c 
                          if = 1, continuation of an integration using an input 
		          created from the outputfile of a previous run           */

  GAS_TORQUE = 1; /*  If 0, don;t read in gas profiles, if 1, gas disk exerts a toque on the planets  */
  UPDATE_OPACITY=1;

  int holdp=0;

  BODY_CREATE=0;  //if BODY_CREATE=1 new BHs will be created every creation_period

  fprintf(fout3, "body1  body2  mass1  mass2  time[days]  current dist  min dist  vrel-x  vrel-y  vrel-z   KE/BE\n\n");
  fflush(fout3);

if (GAS_TORQUE ==1){

 itgas=0;
 init_time(fin_time, &tgas);
  srand(57);
  init_radii(tgas,&Radii,&Sigma,&Alpha,&Temp,&Beta); // gas disk properties. note radii in
						     // AU now, sigma in Msun/AU2

int radp;
for(radp=0;radp<10;radp++){
}

int hind;
//read in aspect ratios
for (hind=0;hind<NRAD;hind++) {
           fscanf(fin2, "%lg", &ar);
	   aspR[hind]=ar;
}

int tind;
//read in aspect ratios
for (tind=0;tind<NRAD;tind++) {
           fscanf(fin3, "%lg", &tau);
           taup[tind]=tau;
} 
 
 R_inner=Radii[0];
 R_outer=Radii[NRAD-1];
 RESCMIN = R_inner-0.00001;               /* inner boundary of the computation domain */
 RESCMAX = R_outer+0.00001;
 D_r=(R_outer-R_inner)/(1.0*NRAD);
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
     } 
     } else if (CONTINUATION == 0) {

             /* New run from an input file */
	     EPOCH = 0.0;
             for (I=0;I<NBOD;I++) {
	          fscanf(fin1,"%lg %lg %lg %lg %lg %lg %lg",&m,&ax,&ex,&inc,&MAn,&omP,&OmN);
	          printf("%lg %lg %lg %lg %lg %lg %lg\n",m,ax,ex,inc,MAn,omP,OmN);
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
	
		  index[I] = I;
	          BODY[I] = m;  
	    }

	
  }
      
  for(i=0;i<MAXNBOD;i++){
    nclose[i]=0;
  }

  if(BODY_CREATE==1){
    fin4 = fopen("test/bcreate_mass.dat","r"); //this file contains the masses of created BHs
    body_creation_times(SEED2,creation_period,final_creation_time,Times);
  }

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
     

      Integrator(&holdp,saving_period,BODY0,NBOD,BODY,Radii,Sigma,aspR,taup,Alpha,Temp,Beta,ybe,xbe,htrybe,sma,yki,&xki,&hkovki,&F_stoc_out);


	 /*re-initialization of the gobal variabes if collision occured*/
	 if (COLLISION == 1) {
	     collision(xki,NBOD,BODY,yki,index,&NEW_NBOD,yki_after,BODY_AFTER,&index_new,ybe,Radii,aspR,taup,Alpha,Sigma,Temp,Beta,BODY0);
	   
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
      if(BODY_CREATE==1){
         if (xbe < Times[nr_of_times]  && Times[nr_of_times] <= xki) {
             create_a_body(NBOD,Sma_min,Inc_max,Ecc_max,BODY,yki,&xx,&yy,&zz,&vxx,&vyy,&vzz,&m_cr);
         
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
             NBOD ++;
             nv = 6*NBOD;
             
         }
       }

	/*output of the data*/
	//now prints every x timesteps set by saving_period
	if(PRINT_TIME==1){
		print_elements(xbe,xki,NBOD,BODY0,BODY,yki,GRADUAL,index); //worrried something is wrong here b/c changing its resolution changes the collision times
	}

	UPDATE_OPACITY=1;

	/*re-initialization of the input data for the next integration step*/
	 
	 xbe = xki;
	 htrybe = hkovki;
	 
         for (i=0; i<nv; i++)
             ybe[i] = yki[i];
   } while (xbe < max_time);
      
 fclose(fin1);
 fclose(fin2);
 fclose(fin3);
 fclose(fin4);
 fclose(fin_time);
 fclose(ptr_myfile);
 fclose(fout3);
}
