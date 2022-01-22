/* use piecewise polytrope EoS */
/* Wolfgang Tichy, Tim Dietrich 5/2014 */

#include "sgrid.h"
#include "EoS_T0.h"

/* functions in this file 
   PwP_polytrope_of_hm1
   PwP_polytrope_rho0_of_hm1
   PwP_polytrope_P_of_hm1
   PwP_polytrope_hm1_of_P
   PwP_init_...
 */

#define Power pow

/* global vars the describe polytrope pieces */
/* we switch from piece i to piece i+1 if rho0>=PwP_rho0[i+1] */
int     PwP_n_pieces = 0; /* number of pieces, i.e. how many different n */
double *PwP_rho0 = NULL;  /* PwP_rho0[0]=0 always,  PwP_rho0[i] is rho0 where we switch */
double *PwP_kappa = NULL; /* PwP_kappa[i] is kappa in piece i */
double *PwP_n = NULL;     /* PwP_n[i] is n in piece i */
double *PwP_k = NULL;     /* PwP_k[i] is k in piece i */
double *PwP_q = NULL;     /* PwP_q[i] = q(PwP_rho0[i]),  q=h-1 */
double *PwP_P = NULL;     /* PwP_P[i] = P(PwP_rho0[i]) */

/* allocate memory for global PwP variables. */
void alloc_PwP_globals(int pieces)
{
  int maxpieces = pieces + 3;
  PwP_rho0  = calloc(maxpieces, sizeof(*PwP_rho0));
  PwP_kappa = calloc(maxpieces, sizeof(*PwP_kappa));
  PwP_n     = calloc(maxpieces, sizeof(*PwP_n));
  PwP_k     = calloc(maxpieces, sizeof(*PwP_k));
  PwP_q     = calloc(maxpieces, sizeof(*PwP_q));
  PwP_P     = calloc(maxpieces, sizeof(*PwP_P));
}
/* free  memory for global PwP variables. */
void free_PwP_globals(void)
{
  free(PwP_rho0);
  free(PwP_kappa);
  free(PwP_n);
  free(PwP_k);
  free(PwP_q);
  free(PwP_P);
  PwP_rho0 = PwP_kappa = PwP_n = PwP_k = PwP_q = PwP_P = NULL;
}


/**************************************************************************/
/* functions that deal with pure piecewise polytropes (i.e. cold EoS) */
/**************************************************************************/

/* for pwp */
/* find piece in which q=hm1 lies
   it is inside piece m, if PwP_q[m]<hm1<PwP_q[m+1] */
void PwP_select_polytrope_n_kappa_k_of_hm1(double hm1,
                                           double *n, double *kappa, double *k)
{
  int m;

  /* find m such that PwP_q[m]<hm1<PwP_q[m+1] */
  for(m=0; m<PwP_n_pieces-1; m++)
    if(hm1 < PwP_q[m+1]) break;
  /* Note, if PwP_n_pieces=1 we always get m=0 */

  *n     = PwP_n[m];
  *kappa = PwP_kappa[m];
  *k     = PwP_k[m];
}

/* for pwp select based on P */
void PwP_select_polytrope_n_kappa_k_of_P(double P,
                                         double *n, double *kappa, double *k)
{
  int m;

  /* find m such that PwP_P[m]<P<PwP_P[m+1] */
  for(m=0; m<PwP_n_pieces-1; m++)
    if(P < PwP_P[m+1]) break;
  /* Note, if PwP_n_pieces=1 we always get m=0 */

  *n     = PwP_n[m];
  *kappa = PwP_kappa[m];
  *k     = PwP_k[m];
}

/* PwP */
/* inputs: hm1 := h-1 , n=1/(Gamma-1) , kappa, constant k
   here: P = kappa rho0^Gamma 
         rho0 espilon = P/(Gamma-1) + k P^(1/Gamma)
   outputs: rest mass density rho0, pressure P, total energy density rhoE and polytropic index n*/
void PwP_polytrope_of_hm1(double hm1,
                          double *rho0, double *P, double *rhoE,
                          double *drho0dhm1)
{
  double hmk, n,kappa,k;

  PwP_select_polytrope_n_kappa_k_of_hm1(hm1, &n, &kappa, &k);
 	 
  hmk = hm1 + (1.0 - k);  /* hmk = hm1 + 1 - k */
  *rho0 = pow(hmk/(kappa*(n+1.0)), n);
  *P    = (*rho0) * hmk/(n+1.0);
//  *rhoE = (hm1+1.0) * (*rho0) - (*P);
  *rhoE = n*(*P) + k*(*rho0); 
  /* drho0/dhm1 = (n/(kappa*(n+1))) [hmk/(kappa*(n+1))]^{n-1} */
  if(*P>0.0)
    *drho0dhm1 = (*rho0)*(*rho0)*n/((*P)*(n+1.0));
  else /* limit for hmk->0 */
  {
    if(n<1.0) *drho0dhm1 = (n/(kappa*(n+1.0)))*pow(hmk/(kappa*(n+1.0)), n-1.0);
    else if(n==1.0) *drho0dhm1 = (n/(kappa*(n+1.0)));
    else            *drho0dhm1 = 0.0;
  }
}

/* get rho0 only from PwP_polytrope_of_hm1 */
double PwP_polytrope_rho0_of_hm1(double hm1)
{
  double rho0, dummy;
  PwP_polytrope_of_hm1(hm1, &rho0, &dummy, &dummy, &dummy);
  return rho0;
}
/* get P only from PwP_polytrope_of_hm1 */
double PwP_polytrope_P_of_hm1(double hm1)
{
  double rho0, P, dummy;
  PwP_polytrope_of_hm1(hm1, &rho0, &P, &dummy, &dummy);
  return P;
}

/* get hm1=h-1 from P */
/* P = kappa^{-n} [hmk/(n+1)]^{n+1}
   h-k = hmk = (n+1) kappa (P/kappa)^{1/(n+1)}  */
double PwP_polytrope_hm1_of_P(double P)
{
	
  double n, k, kappa, hmk; 

  PwP_select_polytrope_n_kappa_k_of_P(P, &n, &kappa, &k);

  hmk = (n+1.0) * kappa * pow(P/kappa, 1.0/(n+1.0));
  return hmk + (k-1.0);
}

/* get rho0, rhoE from P */
/* rho0 = pow(P/kappa, n/(n+1) 
   rhoE = n*P + k*rho0          */
void PwP_polytrope_rho0_rhoE_of_P(double P, double *rho0, double *rhoE)
{
  double n, k, kappa; 
  int m;

  PwP_select_polytrope_n_kappa_k_of_P(P, &n, &kappa, &k);

  *rho0 = pow(P/kappa, n/(n+1));
  *rhoE = n*P + k*(*rho0);
}

/* test if we get the correct n,kappa,k for different q and P */
void PwP_test_piecewise_n_kappa_k()
{
  int m;
  double hm1, P, n,kappa,k;

  printf("values of n,kappa,k for different test q and P:\n");
  for(m=-1; m<PwP_n_pieces; m++)
  {
    if(m<0) hm1 = P = -0.1;
    else
    {
      hm1 = 0.5*(PwP_q[m]+PwP_q[m+1]);
      P   = 0.5*(PwP_P[m]+PwP_P[m+1]);
    }
    PwP_select_polytrope_n_kappa_k_of_hm1(hm1, &n, &kappa, &k);
    printf("q=%e -> n=%e kappa=%e k=%e\n", hm1, n,kappa,k);
    PwP_select_polytrope_n_kappa_k_of_P(P, &n, &kappa, &k);
    printf("P=%e -> n=%e kappa=%e k=%e\n", P, n,kappa,k);
  }
}


/* initialize piecewise polytropes from files in Tim's format */
int PwP_init_file(void)
{
  int i;
  char *fname = Gets("DNSdata_EoS_file");
  char sline[1024];
  int counter = 0;

  printf(" reading PWP pars from file:\n    %s\n", fname); 

  // open file and read the data in
  FILE *ifp = fopen(fname, "r");  
  if (!ifp) errorexits(" could not open file %s", fname); 

  while (fgets(sline,256,ifp)) 
  {
    if (sline[0]=='#') {
      // skip comment
    }
    else if(counter==0){
      sscanf(sline,"%d", &PwP_n_pieces);
      free_PwP_globals();
      alloc_PwP_globals(PwP_n_pieces); /*get mem and set everything to zero */
      counter = counter + 1;
    }
    else if(!(counter==0)){
        sscanf(sline,"%le %le %le",&PwP_rho0[counter-1],&PwP_kappa[counter-1], &PwP_n[counter-1]); 
        counter = counter + 1;
    }
  }
  fclose (ifp);
 
  if(!((counter-1)==PwP_n_pieces)) errorexit("input EoS_file wrong");
  /* done with reading*/

  /* activate PwP_pwp*/
  //PwP_pwp  = 1;
  PwP_k[0] = 1.;

  printf("PwP_init_file: using %d polytrope pieces \n", PwP_n_pieces);
  printf("rho0           q            P            kappa        n            k \n");

  /* compute the constant k, k = a+1 according to arXiv:0812.2163*/
  for(i=0;i<PwP_n_pieces;i++){
    if( i > 0) {    
      PwP_k[i] = ((PwP_k[i-1]*PwP_rho0[i]+PwP_n[i-1]
               * PwP_kappa[i-1]*pow(PwP_rho0[i],1.+1./PwP_n[i-1]))/PwP_rho0[i])
               - PwP_kappa[i]*PwP_n[i]*pow(PwP_rho0[i],1./PwP_n[i]);
    }
    PwP_q[i]  = PwP_k[i] + (PwP_n[i]+1)*PwP_kappa[i]*pow(PwP_rho0[i],1./PwP_n[i])-1.;
    PwP_P[i]  = PwP_kappa[i]*pow(PwP_rho0[i],1.+1./(PwP_n[i]));
    
    /* print rho, kappa, n, k in log-file*/ 
    printf("%e %e %e %e %e %e \n", 
    PwP_rho0[i], PwP_q[i], PwP_P[i], PwP_kappa[i], PwP_n[i], PwP_k[i]);
  }

  /* copy last column and set large end-value*/  
  PwP_rho0[PwP_n_pieces]  = PwP_rho0[(PwP_n_pieces-1)]  + 1.e10;
  PwP_q[PwP_n_pieces]     = PwP_q[(PwP_n_pieces-1)]     + 1.e10;
  PwP_P[PwP_n_pieces]     = PwP_P[(PwP_n_pieces-1)]     + 1.e10;
  PwP_kappa[PwP_n_pieces] = PwP_kappa[(PwP_n_pieces-1)];
  PwP_n[PwP_n_pieces]     = PwP_n[(PwP_n_pieces-1)];
  PwP_k[PwP_n_pieces]     = PwP_k[(PwP_n_pieces-1)];
  
  printf("%e %e %e %e %e %e \n", 
  PwP_rho0[PwP_n_pieces], PwP_q[PwP_n_pieces], PwP_P[PwP_n_pieces], 
  PwP_kappa[PwP_n_pieces], PwP_n[PwP_n_pieces], PwP_k[PwP_n_pieces]);
    

  /* Set parameters, 
     be careful that kappa is set correct, we do not recompute it, 
     when we use a file, but the EoS uses only kappa[0] 
     best: use notebook in /DNSdata/pwpfits */
  char par[1000];

  Setd("EoS_PwP_kappa",PwP_kappa[0]);

  sprintf(par,"%15.14f",PwP_n[0]);
  for(i=1;i<PwP_n_pieces;i++) sprintf(par,"%s %15.14f",par, PwP_n[i]);
  Sets("EoS_PwP_n",par);

  /* we do not use:  sprintf(par,"%15.14f", PwP_rho0[0]); 
     it's redundant. */
  sprintf(par,"%15.14f", PwP_rho0[1]);
  for(i=2;i<PwP_n_pieces;i++) sprintf(par,"%s %15.14f",par, PwP_rho0[i]);
  Sets("EoS_PwP_rho0",par);

  return 0;
}

/* initialize from parameters in sgrid parfile */
int PwP_init_parameter(void)
{
  int i;
  int rho0count, ncount;
  char *str, *s;

  /* count words in par EoS_PwP_n */
  str = strdup(Gets("EoS_PwP_n"));
  ncount = 0;
  for(s=strtok(str, " "); s!=NULL; s=strtok(NULL, " "))
    ncount++;
  free(str);

  /* initialize mem, and set everything to zero */
  free_PwP_globals();
  alloc_PwP_globals(ncount);

  /* activate PwP_pwp*/
  PwP_k[0] = 1.;

  /* loop over pars we read */
  for (i=1; i<3; i++)
  {
    int count, start;
    char str[200];
    char *par;
    double val;

    start = 0;
    count = 0;
    if(i==1) par = Gets("EoS_PwP_n");
    if(i==2) par = Gets("EoS_PwP_rho0");
    while(sscanf(par+start, "%s", str)==1)
    {
      start += strlen(str);
      if(par[start]==' ') start++;
      val = atof(str);
      if(i==1) PwP_n[count] = val;
      if(i==2)
      {
        /* if there is no leading 0 in EoS_PwP_rho0, skip to next PwP_rho0 */
        if(count==0 && val!=0.0) count++;
        PwP_rho0[count] = val;
      }
      count++;
    }
    if(i==1)  ncount    = count;
    if(i==2)  rho0count = count;
  }

  PwP_n_pieces = ncount;
  PwP_kappa[0] = Getd("EoS_PwP_kappa");

  printf("PwP_init_parameter: using %d polytrope pieces \n", PwP_n_pieces);
  printf("rho0           q            P            kappa        n            k \n");
 
  /* compute the constant k, k = a+1 according to arXiv:0812.2163*/
  for(i=0;i<PwP_n_pieces;i++)
  {
    if(i > 0)
    {  
      PwP_kappa[i] = PwP_kappa[i-1]*pow(PwP_rho0[i],(PwP_n[i]-PwP_n[i-1])/
                     (PwP_n[i-1]*PwP_n[i]));

      PwP_k[i] = ((PwP_k[i-1]*PwP_rho0[i]+PwP_n[i-1]
               * PwP_kappa[i-1]*pow(PwP_rho0[i],1.+1./PwP_n[i-1]))/PwP_rho0[i])
               - PwP_kappa[i]*PwP_n[i]*pow(PwP_rho0[i],1./PwP_n[i]);
    }
    PwP_q[i]  = PwP_k[i] + (PwP_n[i]+1)*PwP_kappa[i]*pow(PwP_rho0[i],1./PwP_n[i])-1.;
    PwP_P[i]  = PwP_kappa[i]*pow(PwP_rho0[i],1.+1./(PwP_n[i]));
    
    /* print rho, kappa, n, k in log-file*/ 
    printf("%e %e %e %e %e %e \n", 
    PwP_rho0[i], PwP_q[i], PwP_P[i], PwP_kappa[i], PwP_n[i], PwP_k[i]);
  }

  /* Copy last column and set large end-value.
     This is probably not needed at all. */  
  PwP_rho0[PwP_n_pieces]  = PwP_rho0[(PwP_n_pieces-1)]  + 1.e10;
  PwP_q[PwP_n_pieces]     = PwP_q[(PwP_n_pieces-1)]     + 1.e10;
  PwP_P[PwP_n_pieces]     = PwP_P[(PwP_n_pieces-1)]     + 1.e10;
  PwP_kappa[PwP_n_pieces] = PwP_kappa[(PwP_n_pieces-1)];
  PwP_n[PwP_n_pieces]     = PwP_n[(PwP_n_pieces-1)];
  PwP_k[PwP_n_pieces]     = PwP_k[(PwP_n_pieces-1)];
  /*
  printf("after last piece:\n");
  printf("%e %e %e %e %e %e\n",
  PwP_rho0[PwP_n_pieces], PwP_q[PwP_n_pieces], PwP_P[PwP_n_pieces], 
  PwP_kappa[PwP_n_pieces], PwP_n[PwP_n_pieces], PwP_k[PwP_n_pieces]);
  */
  PwP_test_piecewise_n_kappa_k();

  printf("\n");
  return 0;
}

/* free all PwP stuff */
int PwP_finalize(tGrid *grid)
{
  free_PwP_globals();
  return 0;
}
