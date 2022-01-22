/* sgrid_EoS_T0.h */
/* Wolfgang Tichy, Jan 2022 */


/* struct that contains important EoS info and function pointers */
typedef struct tEOS_T0 {
  /* function pointers to compute EoS vars from different inputs */
  void (*vars_from_hm1)(double hm1,
                        double *rho0, double *P, double *rhoE,
                        double *drho0dhm1);
  double (*rho0_of_hm1)(double hm1);
  double (*hm1_of_P)(double P);
  void (*rho0_rhoE_from_P)(double P, double *rho0, double *rhoE);
} tEoS_T0;
