#include <iostream>
#include "tools.h"
#include "math.h"

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
  int acc = 0;
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

  if(!px && !py) {
    std::cout << "Division by zero" << std::endl;
    return Hj;
  }

  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py * (vx * py - vy * px)/ c3, px * ( vy* px - vx * py) / c3, px/c2, py/c2;

  return Hj;
}

VectorXd Tools::ConvertPolarToCartesian(const VectorXd& polar) {
    VectorXd cartesian(4);

    float rho = polar(0);
    float phi = polar(1);
    // velocity
    float rho_v = polar(2);

    float px = rho * cos(phi);
    float py = rho * sin(phi);
    float vx = rho_v * cos(phi);
    float vy = rho_v * sin(phi);

    cartesian << px, py, vx, vy;
    return cartesian;
}

VectorXd Tools::ConvertCartesianToPolar(const Eigen::VectorXd& cartesian) {
    VectorXd polar(3);

    float px = cartesian(0);
    float py = cartesian(1);
    float vx = cartesian(1);
    float vy = cartesian(2);

    float rho = sqrt(px * px + py * py);
    float phi = atan2(py, px);
    float rho_v = (px * vx + py * vy)/rho;

    polar << rho, phi, rho_v;
    return polar;
}
