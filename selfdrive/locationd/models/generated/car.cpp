
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3550283937265723398) {
   out_3550283937265723398[0] = delta_x[0] + nom_x[0];
   out_3550283937265723398[1] = delta_x[1] + nom_x[1];
   out_3550283937265723398[2] = delta_x[2] + nom_x[2];
   out_3550283937265723398[3] = delta_x[3] + nom_x[3];
   out_3550283937265723398[4] = delta_x[4] + nom_x[4];
   out_3550283937265723398[5] = delta_x[5] + nom_x[5];
   out_3550283937265723398[6] = delta_x[6] + nom_x[6];
   out_3550283937265723398[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4208484015671583759) {
   out_4208484015671583759[0] = -nom_x[0] + true_x[0];
   out_4208484015671583759[1] = -nom_x[1] + true_x[1];
   out_4208484015671583759[2] = -nom_x[2] + true_x[2];
   out_4208484015671583759[3] = -nom_x[3] + true_x[3];
   out_4208484015671583759[4] = -nom_x[4] + true_x[4];
   out_4208484015671583759[5] = -nom_x[5] + true_x[5];
   out_4208484015671583759[6] = -nom_x[6] + true_x[6];
   out_4208484015671583759[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_299537516037088952) {
   out_299537516037088952[0] = 1.0;
   out_299537516037088952[1] = 0.0;
   out_299537516037088952[2] = 0.0;
   out_299537516037088952[3] = 0.0;
   out_299537516037088952[4] = 0.0;
   out_299537516037088952[5] = 0.0;
   out_299537516037088952[6] = 0.0;
   out_299537516037088952[7] = 0.0;
   out_299537516037088952[8] = 0.0;
   out_299537516037088952[9] = 1.0;
   out_299537516037088952[10] = 0.0;
   out_299537516037088952[11] = 0.0;
   out_299537516037088952[12] = 0.0;
   out_299537516037088952[13] = 0.0;
   out_299537516037088952[14] = 0.0;
   out_299537516037088952[15] = 0.0;
   out_299537516037088952[16] = 0.0;
   out_299537516037088952[17] = 0.0;
   out_299537516037088952[18] = 1.0;
   out_299537516037088952[19] = 0.0;
   out_299537516037088952[20] = 0.0;
   out_299537516037088952[21] = 0.0;
   out_299537516037088952[22] = 0.0;
   out_299537516037088952[23] = 0.0;
   out_299537516037088952[24] = 0.0;
   out_299537516037088952[25] = 0.0;
   out_299537516037088952[26] = 0.0;
   out_299537516037088952[27] = 1.0;
   out_299537516037088952[28] = 0.0;
   out_299537516037088952[29] = 0.0;
   out_299537516037088952[30] = 0.0;
   out_299537516037088952[31] = 0.0;
   out_299537516037088952[32] = 0.0;
   out_299537516037088952[33] = 0.0;
   out_299537516037088952[34] = 0.0;
   out_299537516037088952[35] = 0.0;
   out_299537516037088952[36] = 1.0;
   out_299537516037088952[37] = 0.0;
   out_299537516037088952[38] = 0.0;
   out_299537516037088952[39] = 0.0;
   out_299537516037088952[40] = 0.0;
   out_299537516037088952[41] = 0.0;
   out_299537516037088952[42] = 0.0;
   out_299537516037088952[43] = 0.0;
   out_299537516037088952[44] = 0.0;
   out_299537516037088952[45] = 1.0;
   out_299537516037088952[46] = 0.0;
   out_299537516037088952[47] = 0.0;
   out_299537516037088952[48] = 0.0;
   out_299537516037088952[49] = 0.0;
   out_299537516037088952[50] = 0.0;
   out_299537516037088952[51] = 0.0;
   out_299537516037088952[52] = 0.0;
   out_299537516037088952[53] = 0.0;
   out_299537516037088952[54] = 1.0;
   out_299537516037088952[55] = 0.0;
   out_299537516037088952[56] = 0.0;
   out_299537516037088952[57] = 0.0;
   out_299537516037088952[58] = 0.0;
   out_299537516037088952[59] = 0.0;
   out_299537516037088952[60] = 0.0;
   out_299537516037088952[61] = 0.0;
   out_299537516037088952[62] = 0.0;
   out_299537516037088952[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5631650506502684714) {
   out_5631650506502684714[0] = state[0];
   out_5631650506502684714[1] = state[1];
   out_5631650506502684714[2] = state[2];
   out_5631650506502684714[3] = state[3];
   out_5631650506502684714[4] = state[4];
   out_5631650506502684714[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5631650506502684714[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5631650506502684714[7] = state[7];
}
void F_fun(double *state, double dt, double *out_3059307807673509865) {
   out_3059307807673509865[0] = 1;
   out_3059307807673509865[1] = 0;
   out_3059307807673509865[2] = 0;
   out_3059307807673509865[3] = 0;
   out_3059307807673509865[4] = 0;
   out_3059307807673509865[5] = 0;
   out_3059307807673509865[6] = 0;
   out_3059307807673509865[7] = 0;
   out_3059307807673509865[8] = 0;
   out_3059307807673509865[9] = 1;
   out_3059307807673509865[10] = 0;
   out_3059307807673509865[11] = 0;
   out_3059307807673509865[12] = 0;
   out_3059307807673509865[13] = 0;
   out_3059307807673509865[14] = 0;
   out_3059307807673509865[15] = 0;
   out_3059307807673509865[16] = 0;
   out_3059307807673509865[17] = 0;
   out_3059307807673509865[18] = 1;
   out_3059307807673509865[19] = 0;
   out_3059307807673509865[20] = 0;
   out_3059307807673509865[21] = 0;
   out_3059307807673509865[22] = 0;
   out_3059307807673509865[23] = 0;
   out_3059307807673509865[24] = 0;
   out_3059307807673509865[25] = 0;
   out_3059307807673509865[26] = 0;
   out_3059307807673509865[27] = 1;
   out_3059307807673509865[28] = 0;
   out_3059307807673509865[29] = 0;
   out_3059307807673509865[30] = 0;
   out_3059307807673509865[31] = 0;
   out_3059307807673509865[32] = 0;
   out_3059307807673509865[33] = 0;
   out_3059307807673509865[34] = 0;
   out_3059307807673509865[35] = 0;
   out_3059307807673509865[36] = 1;
   out_3059307807673509865[37] = 0;
   out_3059307807673509865[38] = 0;
   out_3059307807673509865[39] = 0;
   out_3059307807673509865[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3059307807673509865[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3059307807673509865[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3059307807673509865[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3059307807673509865[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3059307807673509865[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3059307807673509865[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3059307807673509865[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3059307807673509865[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3059307807673509865[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3059307807673509865[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3059307807673509865[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3059307807673509865[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3059307807673509865[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3059307807673509865[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3059307807673509865[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3059307807673509865[56] = 0;
   out_3059307807673509865[57] = 0;
   out_3059307807673509865[58] = 0;
   out_3059307807673509865[59] = 0;
   out_3059307807673509865[60] = 0;
   out_3059307807673509865[61] = 0;
   out_3059307807673509865[62] = 0;
   out_3059307807673509865[63] = 1;
}
void h_25(double *state, double *unused, double *out_2009596230628797542) {
   out_2009596230628797542[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3504271223961127907) {
   out_3504271223961127907[0] = 0;
   out_3504271223961127907[1] = 0;
   out_3504271223961127907[2] = 0;
   out_3504271223961127907[3] = 0;
   out_3504271223961127907[4] = 0;
   out_3504271223961127907[5] = 0;
   out_3504271223961127907[6] = 1;
   out_3504271223961127907[7] = 0;
}
void h_24(double *state, double *unused, double *out_4561368411239011603) {
   out_4561368411239011603[0] = state[4];
   out_4561368411239011603[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3502069605190963543) {
   out_3502069605190963543[0] = 0;
   out_3502069605190963543[1] = 0;
   out_3502069605190963543[2] = 0;
   out_3502069605190963543[3] = 0;
   out_3502069605190963543[4] = 1;
   out_3502069605190963543[5] = 0;
   out_3502069605190963543[6] = 0;
   out_3502069605190963543[7] = 0;
   out_3502069605190963543[8] = 0;
   out_3502069605190963543[9] = 0;
   out_3502069605190963543[10] = 0;
   out_3502069605190963543[11] = 0;
   out_3502069605190963543[12] = 0;
   out_3502069605190963543[13] = 1;
   out_3502069605190963543[14] = 0;
   out_3502069605190963543[15] = 0;
}
void h_30(double *state, double *unused, double *out_5949752251454604008) {
   out_5949752251454604008[0] = state[4];
}
void H_30(double *state, double *unused, double *out_930216555814872301) {
   out_930216555814872301[0] = 0;
   out_930216555814872301[1] = 0;
   out_930216555814872301[2] = 0;
   out_930216555814872301[3] = 0;
   out_930216555814872301[4] = 1;
   out_930216555814872301[5] = 0;
   out_930216555814872301[6] = 0;
   out_930216555814872301[7] = 0;
}
void h_26(double *state, double *unused, double *out_6893262104356266443) {
   out_6893262104356266443[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2198035452006286762) {
   out_2198035452006286762[0] = 0;
   out_2198035452006286762[1] = 0;
   out_2198035452006286762[2] = 0;
   out_2198035452006286762[3] = 0;
   out_2198035452006286762[4] = 0;
   out_2198035452006286762[5] = 0;
   out_2198035452006286762[6] = 0;
   out_2198035452006286762[7] = 1;
}
void h_27(double *state, double *unused, double *out_7011269498289643187) {
   out_7011269498289643187[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4040991950962615117) {
   out_4040991950962615117[0] = 0;
   out_4040991950962615117[1] = 0;
   out_4040991950962615117[2] = 0;
   out_4040991950962615117[3] = 1;
   out_4040991950962615117[4] = 0;
   out_4040991950962615117[5] = 0;
   out_4040991950962615117[6] = 0;
   out_4040991950962615117[7] = 0;
}
void h_29(double *state, double *unused, double *out_2912966200010271323) {
   out_2912966200010271323[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8550804458305265301) {
   out_8550804458305265301[0] = 0;
   out_8550804458305265301[1] = 1;
   out_8550804458305265301[2] = 0;
   out_8550804458305265301[3] = 0;
   out_8550804458305265301[4] = 0;
   out_8550804458305265301[5] = 0;
   out_8550804458305265301[6] = 0;
   out_8550804458305265301[7] = 0;
}
void h_28(double *state, double *unused, double *out_90904213451541964) {
   out_90904213451541964[0] = state[5];
   out_90904213451541964[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3528484607297725079) {
   out_3528484607297725079[0] = 0;
   out_3528484607297725079[1] = 0;
   out_3528484607297725079[2] = 0;
   out_3528484607297725079[3] = 0;
   out_3528484607297725079[4] = 0;
   out_3528484607297725079[5] = 1;
   out_3528484607297725079[6] = 0;
   out_3528484607297725079[7] = 0;
   out_3528484607297725079[8] = 0;
   out_3528484607297725079[9] = 0;
   out_3528484607297725079[10] = 0;
   out_3528484607297725079[11] = 0;
   out_3528484607297725079[12] = 0;
   out_3528484607297725079[13] = 0;
   out_3528484607297725079[14] = 1;
   out_3528484607297725079[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
