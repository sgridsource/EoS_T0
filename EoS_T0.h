/* EoS_T0.h */
/* Wolfgang Tichy, Jan 2022 */


/* from PwP.c */
void PwP_polytrope_of_hm1(double hm1,
                          double *rho0, double *P, double *rhoE,
                          double *drho0dhm1);
double PwP_polytrope_rho0_of_hm1(double hm1);
double PwP_polytrope_P_of_hm1(double hm1);
double PwP_polytrope_hm1_of_P(double P);
void PwP_polytrope_rho0_rhoE_of_P(double P, double *rho0, double *rhoE);
int PwP_init_file(void);
int PwP_init_parameter(void);
int PwP_poly_init(void);
int PwP_finalize(tGrid *grid);


/* from tab1d_Of_rho0_AtT0.c */
void EoS_tab1d_load_rho0_epsl_P_AtT0(char *fname);
int tab1d_Of_P_AtT0(double P, double *rho0, double *epsl,
                    double *dPdrho0, double *dPdepsl);
int tab1d_Of_hm1_AtT0(double hm1, double *rho0, double *epsl,
                      double *P, double *dPdrho0, double *dPdepsl);
void tab1d_rho0_epsl_P_drho0dhm1_Of_hm1_AtT0(double hm1,
                                             double *rho0, double *epsl,
                                             double *P, double *drho0dhm1);
void EoS_free_tab1d_rho0(void);

/* from EoS_T0.c */
double hm1_of_rho0_epsl_P(double rho0, double epsl, double P);
void tab1d_rho0_P_rhoE_drho0dhm1_from_hm1(double hm1, double *rho0,
                                          double *P, double *rhoE,
                                          double *drho0dhm1);
double tab1d_rho0_of_hm1(double hm1);
double tab1d_hm1_of_P(double P);
void tab1d_rho0_rhoE_from_P(double P, double *rho0, double *rhoE);
int EoS_T0_init_from_pars(tGrid *grid);
int EoS_T0_finalize(tGrid *grid);
