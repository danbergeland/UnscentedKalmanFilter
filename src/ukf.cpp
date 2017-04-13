#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::iostream;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0,    0,  0,
        0, 1, 0,    0,  0,
        0, 0, 100,  0,  0,
        0, 0, 0,  100,  0,
        0, 0, 0,    0,100;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //30 in template
  std_a_ = .01;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //30 in template
  std_yawdd_ = 0.0;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.05;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.05;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.1;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.0175;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.1;


  Q_ = MatrixXd(2,2);
  Q_ << std_a_*std_a_,0,
              0,  std_yawdd_*std_yawdd_;
              
  R_ = MatrixXd(3,3);
  R_ << std_radr_*std_radr_,0,0,
              0,std_radphi_*std_radphi_,0,
              0,0, std_radrd_*std_radrd_;
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  

  time_us_ = 0;
  
    ///* State dimension
  n_x_ = 5;
  n_z_ = 3;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_,1+2*n_aug_);

  weights_ = VectorXd(1+2*n_aug_);
  //Augmented laser matrix
  Xsig_aug_ = MatrixXd(n_aug_,1+2*n_aug_);
  
  z_ = VectorXd(n_z_);
  
  Zsig_ = MatrixXd(n_z_,1+2*n_aug_);
  
  S_ = MatrixXd(n_z_,n_z_);
  ///* the current NIS for radar
  NIS_radar_ = 0;

  ///* the current NIS for laser
  NIS_laser_ = 0;
  
    //set weights, they are based on unchanging parameters, so we only need this once.
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) 
  {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  bool predict = false;
  
  if(is_initialized_ == false)
  {

    
    time_us_ = meas_package.timestamp_;
    
    //initialize with laser
    if(use_laser_ && meas_package.sensor_type_ ==  MeasurementPackage::SensorType::LASER)
    {
      std::cout<<"Initialize with Laser\r\n";
      x_ << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),0,0,0;
      std::cout<<x_<<"\r\n";
      predict = true;
      is_initialized_ = true;
    }
    //initialize with radar
    else if(use_radar_ && meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
    {
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);
      float px = cos(phi)*rho;
      float py = sin(phi)*rho;
      std::cout<<"Initialize with Radar\r\n";
      x_ << px,py,0,0,0;
      std::cout<<x_<<"\r\n";
      predict = true;
      is_initialized_ = true;
    }

  }
  //If initialized, run normal updates
  else
  {
    if((use_laser_==true) && (meas_package.sensor_type_ ==  MeasurementPackage::SensorType::LASER))
    {
      std::cout<<"State before laser update\n\r"<<x_<<"\r\n";
      UpdateLidar(meas_package);
      predict = true;
    }
    else if((use_radar_==true) && (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR))
    {
      std::cout<<"State before radar update\n\r"<<x_<<"\r\n";
      UpdateRadar(meas_package);
      predict = true;
    }
  }
  //only run prediction if the data was added (to allow for enable of radar/laser)
  if(predict==true)
  {
    float dt = (meas_package.timestamp_-time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;
    std::cout<<"State before predict\n\r"<<x_<<"\r\n";
    Prediction(dt);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //only predict if initialized
  if(is_initialized_)
  {
    AugmentedSigmaPoints();
    SigmaPointPrediction(delta_t);
    PredictMeanAndCovariance();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
  //Update z and S matrices for x and Zsig
  PredictRadarMeasurement();
  
  VectorXd z = VectorXd(3);
  z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),meas_package.raw_measurements_(2);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_;
    //angle normalization
    z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  std::cout<<"K:\n\r"<<K<<"\n\r";
  
  //residual
  VectorXd z_diff = z - z_;

  //angle normalization
  z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

  std::cout << "Z_diff:\n\r"<<z_diff<<"\n\r";
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();
  
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out){
      //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  
  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }
  
    *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(){
   //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug << x_, 0,0;
  Xsig_aug_.col(0) = x_aug;
  //create augmented covariance matrix, using Q 
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5)=std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  //create square root matrix
  MatrixXd pSqrt = P_aug.llt().matrixL();
  double scaler = sqrt(lambda_+n_aug_);
  //create augmented sigma points
  for(int i=0;i<n_aug_;i++)
  {
      Xsig_aug_.col(i+1) =     x_aug + scaler*pSqrt.col(i);
      Xsig_aug_.col(i+1+n_aug_) = x_aug - scaler*pSqrt.col(i);
  }
}

void UKF::SigmaPointPrediction(float dt){
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*dt) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*dt) );
    }
    else {
        px_p = p_x + v*dt*cos(yaw);
        py_p = p_y + v*dt*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*dt;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*dt*dt * cos(yaw);
    py_p = py_p + 0.5*nu_a*dt*dt * sin(yaw);
    v_p = v_p + nu_a*dt;

    yaw_p = yaw_p + 0.5*nu_yawdd*dt*dt;
    yawd_p = yawd_p + nu_yawdd*dt;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

}

void UKF::PredictMeanAndCovariance(){
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  
  //predict state mean
  for(int i=0; i< (2*n_aug_+1);i++)
  {
      x = x + weights_(i)*Xsig_pred_.col(i);
  }
  
  x_ = x;
  
  MatrixXd P = MatrixXd(n_x_, n_x_);
    //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));

    P = P + weights_(i) * x_diff * x_diff.transpose();
  }
  P_ = P;
}

void UKF::PredictRadarMeasurement() {

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    if(p_x>.001 || p_x<-.001)
    {
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
      Zsig(1,i) = atan2(p_y,p_x);                                 //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

//angle normalization
    z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_;
  //write result
  z_ = z_pred;
  std::cout<<"z pred:\r\n"<<z_<<"\r\n";
  S_ = S;
  Zsig_ = Zsig;
}

