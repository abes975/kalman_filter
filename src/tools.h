#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);


  /*
  * A helper method to convert polar coordinates to cartesian coordinates.
  * Used when the source of measurement is RADAR
  */
  Eigen::VectorXd ConvertPolarToCartesian(const Eigen::VectorXd& polar);

  /*
  * A helper method to convert cartesian coordinates to polar.
  * Used when the source of measurement is  RADAR in UpdateEKF
  */
  Eigen::VectorXd ConvertCartesianToPolar(const Eigen::VectorXd& cartesian);

};

#endif /* TOOLS_H_ */
