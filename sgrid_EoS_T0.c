/* sgrid_EoS_T0.c */
/* Wolfgang Tichy, Jan 2022 */

#include "sgrid.h"
#include "EoS_T0.h"


int sgrid_EoS_T0() 
{
  if (!Getv("physics", "EoS_T0")) return 0;
  printf("Adding EoS_T0\n");

  /* functions */
  AddFun(POST_PARAMETERS, EoS_T0_init_from_pars, "init EoS_T0 structure");
  AddFun(POST_FINALIZE_GRID, EoS_T0_finalize, "free mem. in EoS_T0 and PwP");

  /* variables */
  //AddVar("EoS_T0_u", "",  "wave function");
  
  /* parameters */
  AddPar("EoS_type", "PwP", "equation of state we use [PwP,tab1d_AtT0]");
  AddPar("EoS_PwP_n", "1", "polytropic index n, Gamma = 1 + 1/n");
  AddPar("EoS_PwP_kappa",  "1",  "kappa in EOS: P = kappa rho0^Gamma");
  AddPar("EoS_PwP_rho0", "2.36701096e-04 8.11322219e-04 1.61880065e-03",
         "densities where we switch between PwP pieces (used only if "
         "EoS_PwP_n has more than 1 entry)");
  AddPar("EoS_tab1d_load_file", "", "file to load for 1d EoS");

  return 0;
}
