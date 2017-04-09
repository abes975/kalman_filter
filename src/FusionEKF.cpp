#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //state covariance matrix P
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;


  //the initial transition matrix F_
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  Q_ = MatrixXd(4,4);
  Q_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0 ,0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    cout << "EKF: " << endl;
    VectorXd x_(4);
    x_ << 1, 1, 1, 1;
    cout << "FINE EKF: " << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ = tools.ConvertPolarToCartesian(measurement_pack.raw_measurements_);
      Hj_ = tools.CalculateJacobian(x_);

      ekf_.Init(x_, P_, F_, Hj_, R_radar_, Q_);
      cout << "INIT RADAR EKF: " << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "BEFORE INIT LASER EKF: " << endl;
      x_ << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1], 0, 0;
      ekf_.Init(x_, P_, F_, H_laser_, R_laser_,Q_);
      cout << "INIT LASER EKF: " << endl;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //compute the time elapsed between the current and previous measurements
  //dt - expressed in seconds
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  cout << "AFTER INIT BEFORE DT: " << endl;
	previous_timestamp_ = measurement_pack.timestamp_;
  float noise_ax = 9;
	float noise_ay = 9;
  cout << "CREATING F: " << endl;
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0 , dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  cout << "CREATING Q: " << endl;


   ekf_.Q_ << pow(dt,4)/4 * noise_ax, 0, pow(dt,3)/2 * noise_ax, 0,
         0, pow(dt,4)/4 * noise_ay, 0 , pow(dt,3)/2 * noise_ay,
         pow(dt,3)/2 * noise_ax, 0, pow(dt,2) * noise_ax, 0,
         0, pow(dt,3)/2 * noise_ay, 0, pow(dt,2) * noise_ay;

  cout << "BEFORE PREDICT: " << endl;
  ekf_.Predict();
  cout << "END PREDICT: " << endl;
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    //VectorXd z = tools.ConvertPolarToCartesian(measurement_pack.raw_measurements_);
    cout << "MEASURE RADAR EKF: " << endl;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    cout << "FINISH MEASURE RADAR EKF: " << endl;
  } else {
    cout << "MEASURE LASER EKF: " << endl;
    // Laser updates
    ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
    cout << "FINISH MEASURE LASER EKF: " << endl;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
