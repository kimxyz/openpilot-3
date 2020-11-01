/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6250306268315751845);
void inv_err_fun(double *nom_x, double *true_x, double *out_715827553028028078);
void H_mod_fun(double *state, double *out_2703875983461333353);
void f_fun(double *state, double dt, double *out_7072834655337394657);
void F_fun(double *state, double dt, double *out_2271556917003262232);
void h_3(double *state, double *unused, double *out_6215483252332514156);
void H_3(double *state, double *unused, double *out_7952728626426121366);
void h_4(double *state, double *unused, double *out_7108172694207513842);
void H_4(double *state, double *unused, double *out_5994334220325727439);
void h_9(double *state, double *unused, double *out_5300014547476410143);
void H_9(double *state, double *unused, double *out_3374200612179985356);
void h_10(double *state, double *unused, double *out_6161794153708182530);
void H_10(double *state, double *unused, double *out_4290626965068938111);
void h_12(double *state, double *unused, double *out_7783598547602621082);
void H_12(double *state, double *unused, double *out_3887778708586517870);
void h_31(double *state, double *unused, double *out_1348777693866794466);
void H_31(double *state, double *unused, double *out_6562868483729937691);
void h_32(double *state, double *unused, double *out_5327306578188609016);
void H_32(double *state, double *unused, double *out_8513370205377689338);
void h_13(double *state, double *unused, double *out_8530893049992882354);
void H_13(double *state, double *unused, double *out_7957139106131299503);
void h_14(double *state, double *unused, double *out_5300014547476410143);
void H_14(double *state, double *unused, double *out_3374200612179985356);
void h_19(double *state, double *unused, double *out_5523785801753562172);
void H_19(double *state, double *unused, double *out_7806880463282416537);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);