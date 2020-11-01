/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3550283937265723398);
void inv_err_fun(double *nom_x, double *true_x, double *out_4208484015671583759);
void H_mod_fun(double *state, double *out_299537516037088952);
void f_fun(double *state, double dt, double *out_5631650506502684714);
void F_fun(double *state, double dt, double *out_3059307807673509865);
void h_25(double *state, double *unused, double *out_2009596230628797542);
void H_25(double *state, double *unused, double *out_3504271223961127907);
void h_24(double *state, double *unused, double *out_4561368411239011603);
void H_24(double *state, double *unused, double *out_3502069605190963543);
void h_30(double *state, double *unused, double *out_5949752251454604008);
void H_30(double *state, double *unused, double *out_930216555814872301);
void h_26(double *state, double *unused, double *out_6893262104356266443);
void H_26(double *state, double *unused, double *out_2198035452006286762);
void h_27(double *state, double *unused, double *out_7011269498289643187);
void H_27(double *state, double *unused, double *out_4040991950962615117);
void h_29(double *state, double *unused, double *out_2912966200010271323);
void H_29(double *state, double *unused, double *out_8550804458305265301);
void h_28(double *state, double *unused, double *out_90904213451541964);
void H_28(double *state, double *unused, double *out_3528484607297725079);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
