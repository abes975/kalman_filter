#include <iostream>
#include "tools.h"
#include "math.h"
#include <cfloat>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  if (!estimations.size() || (estimations.size() != ground_truth.size()))
    return rmse;

  for (int i = 0; i < estimations.size(); ++i) {
      VectorXd res = estimations.at(i) - ground_truth.at(i);
      res = res.array() * res.array();
      rmse += res;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px*px + py*py;
  if(fabs(c1) <= FLT_EPSILON) {
      //std::cout << "Division by zero" << std::endl;
      return Hj;
   }

  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py * (vx * py - vy * px)/ c3, px * ( vy* px - vx * py) / c3, px/c2, py/c2;

  return Hj;
}

VectorXd Tools::ConvertPolarToCartesian(const VectorXd& polar) {
    VectorXd cartesian(4);

    double rho = polar(0);
    double phi = polar(1);
    // velocity
    double rho_v = polar(2);

    double px = rho * cos(phi);
    double py = rho * sin(phi);
    double vx = rho_v * cos(phi);
    double vy = rho_v * sin(phi);

    cartesian << px, py, vx, vy;
    return cartesian;
}

VectorXd Tools::ConvertCartesianToPolar(const VectorXd& cartesian) {
    float px = cartesian(0);
    float py = cartesian(1);
    float vx = cartesian(2);
    float vy = cartesian(3);
    double rho = sqrt(px*px+py*py);

    VectorXd polar(3);
    polar << 0,0,0;

    if (!px) {
        px = FLT_EPSILON;
        rho = sqrt(px*px+py*py);
    }
    if (!py) {
        py = FLT_EPSILON;
        rho = sqrt(px*px+py*py);
    }

    if (!rho) {
        rho = FLT_EPSILON;
    }

    double phi = atan2(py, px);
    double rho_v = (px*vx+py*vy)/rho;
    polar << rho,
             phi,
             rho_v;

    return polar;
}
