/* EoS_T0.c */
/* Wolfgang Tichy 2022 */


#include "sgrid.h"
#include "EoS_T0.h"


/* global struct with EoS info */
tEoS_T0 EoS_T0[1];



/**************************************************************************/
/* some general functions in that hold for any EoS */
/**************************************************************************/

/* return h-1 = epsl + P/rho0 */
double hm1_of_rho0_epsl_P(double rho0, double epsl, double P)
{
  if(rho0==0.) return 0.;
  else         return epsl + P/rho0;
}

/* return h = 1 + epsl + P/rho0 */
double h_of_rho0_epsl_P(double rho0, double epsl, double P)
{
  return hm1_of_rho0_epsl_P(rho0, epsl, P) + 1.;
}

/* return epsl = rhoE/rho0 - 1,  since: rhoE = rho0 * (1 + epsl) */
double epsl_of_rho0_rhoE(double rho0, double rhoE)
{
  if(rho0==0.) return 0.;
  else         return rhoE/rho0 - 1.;
}

/**************************************************************************/
/* wrappers for some funcs in tab1d_Of_rho0_AtT0.c */
/**************************************************************************/

/* for EoS->vars_from_hm1 */
void tab1d_rho0_P_rhoE_drho0dhm1_from_hm1(double hm1, double *rho0,
                                          double *P, double *rhoE,
                                          double *drho0dhm1)
{
  double epsl;

  tab1d_rho0_epsl_P_drho0dhm1_Of_hm1_AtT0(hm1, rho0, &epsl, P,
                                          drho0dhm1);
  /* total energy density */
  *rhoE = *rho0 + (*rho0) * epsl;
}

/* for EoS->rho0_of_hm1 */
double tab1d_rho0_of_hm1(double hm1)
{
  double rho0, epsl, P, dPdrho0, dPdepsl;

  tab1d_Of_hm1_AtT0(hm1, &rho0, &epsl, &P, &dPdrho0, &dPdepsl);
  return rho0;
}

/* for EoS->hm1_of_P */
double tab1d_hm1_of_P(double P)
{
  double rho0, epsl, dPdrho0, dPdepsl;

  tab1d_Of_P_AtT0(P, &rho0, &epsl, &dPdrho0, &dPdepsl);

  return hm1_of_rho0_epsl_P(rho0, epsl, P);
}

/* for EoS->rho0_rhoE_from_P */
void tab1d_rho0_rhoE_from_P(double P, double *rho0, double *rhoE)
{
  double epsl, dPdrho0, dPdepsl;

  tab1d_Of_P_AtT0(P, rho0, &epsl, &dPdrho0, &dPdepsl);

  /* total energy density */
  *rhoE = *rho0 + (*rho0) * epsl;
}


/**************************************************************************/
/* fill EoS structure depending on the EoS we use */
/**************************************************************************/
int EoS_T0_init_from_pars(tGrid *grid)
{
  /* for TOV we need an EoS so we need to init EoS here */
  if(Getv("EoS_type", "PwP"))
  {
    PwP_init_parameter();
    EoS_T0->vars_from_hm1    = PwP_polytrope_of_hm1;
    EoS_T0->rho0_of_hm1      = PwP_polytrope_rho0_of_hm1;
    EoS_T0->hm1_of_P         = PwP_polytrope_hm1_of_P;
    EoS_T0->rho0_rhoE_from_P = PwP_polytrope_rho0_rhoE_of_P;
  }
  else if(Getv("EoS_type", "tab1d_AtT0"))
  {
    EoS_tab1d_load_rho0_epsl_P_AtT0(Gets("EoS_tab1d_load_file"));
    EoS_T0->vars_from_hm1    = tab1d_rho0_P_rhoE_drho0dhm1_from_hm1;
    EoS_T0->rho0_of_hm1      = tab1d_rho0_of_hm1;
    EoS_T0->hm1_of_P         = tab1d_hm1_of_P;
    EoS_T0->rho0_rhoE_from_P = tab1d_rho0_rhoE_from_P;
  }
  else
  {
    errorexit("unkown EoS_type");
  }

  return 0;
}

/* free all EoS stuff */
int EoS_T0_finalize(tGrid *grid)
{
  PwP_finalize(grid);
  EoS_free_tab1d_rho0();
  return 0;
}


/**************************************************************************/
/* simple interface funcs ro access EoS structure */
/**************************************************************************/
void EoS_T0_rho0_P_rhoE_from_hm1(double hm1,
                                 double *rho0, double *P, double *rhoE)
{
  double drho0dhm1[1]; /* not used outside this func */

  EoS_T0->vars_from_hm1(hm1, rho0, P, rhoE, drho0dhm1);
}
