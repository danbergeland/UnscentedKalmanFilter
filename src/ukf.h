#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  
  ///* state vector: predicted z: rho, phi, rho_dot
  VectorXd z_;
  
  ///* Measurement space covariance matrix
  MatrixXd S_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Q matrix for updates
  MatrixXd Q_;
  
  ///*R matrix for radar updates
  MatrixXd R_;

  ///* Weights of sigma points
  VectorXd weights_;
  
  ///*Augmented sigma matrix
  MatrixXd Xsig_aug_;
  
  ///*Zsig radar prediction
  MatrixXd Zsig_;

  ///* State dimension
  int n_x_;
  
  ///* Radar dimension
  int n_z_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
  
   /**
   * Generates Sigma Points, given current x_, P_ and lamda
   * @param Provide a pointer to a MatrixXd.  The Sigma Points matrix will update Xsig_out
   */
  void GenerateSigmaPoints(MatrixXd* Xsig_out);
     /**
   * Generates Augmented Sigma Points, given current x_, P_, lamda and Q_laser

   */
  void AugmentedSigmaPoints(void);
  
   /**
   * Updates class Xsig_pred matrix, based on in augmented sigma matrix.
   * @param Xsig_aug is the augmented laser matrix, dt is the time step in seconds
   */
  void SigmaPointPrediction(float dt);
  
  /**
   * Updates x_ and P_ based on current augmented sigma points matrix
   * 
   */
  void PredictMeanAndCovariance(void);
  
   /**
   * Updates class Xsig_pred matrix, based on in augmented sigma matrix.
   * @param z_out is the predicted radar value, S_out is the predicted covariance
   */
  void PredictRadarMeasurement(void);
};

#endif /* UKF_H */
