# Matrix-Fisher-Gaussian
Distribution and filtering on SO(3) x Euclidean space

## List of parameters

* angular velocity noise parameters
  * **parameters.omegaNoise.randomWalk**: random walk noise for angular velocity
  * **parameters.omegaNoise.biasInstability**: bias random walk noise for angular velocity

* measurement noise parameters
  * **parameters.meaNoise**
    * A 3-by-3 matrix if **meaIsVec**==false, representing the covariance for attitude uncertainty (**GaussMea**==true), or S for attitude uncertainty (**GaussMea**==false).
    * A **nVecRef**-by-1 vector if **meaIsVec**==true, representing the variance for isotropic vector measurement uncertainty (**GaussMea**==true), or kappa for unit vector uncertainty (**GaussMea**==false).
    
* other settings
  * **parameters.setting.omegaLocal**: true if angular velocity is measured in body-fixed frame; false if angular velocity is measured in inertial frame.
  * **parameters.setting.GaussMea**: true if the measurement noise is Gaussian; false if the measurement noise is matrix Fisher (**meaIsVec**==false) or von Mises Fisher (**meaIsVec**==true).
  * **parameters.setting.meaIsVec**: true if measurements are vectors; false if measurement is attitude.
  
* vector measurement settings
  * **parameters.setting.vecRefInertial**: true if reference vector is in inertial frame; false if reference vector is in body-fixed frame.
  * **parameters.setting.nVecRef**: number of reference vectors.
  * **parameters.setting.vRef**: reference vectors, a 3***nVecRef**-by-1 vector.
  
* attitude measurement settings
  * **parameters.setting.attMeaLocal**: true if the noise is in the body-fixed frame; false if the noise is in inertial frame.
  
* initialization settings
  * **parameters.initValue.Miu**: initial bias, which is to be substracted.
  * **parameters.initValue.U**
  * **parameters.initValue.V**
  * **parameters.initValue.xNoise**: covariance for the initial bias
  * **parameters.initValue.RNoise**: covariance for initial attitude (**GaussMea**==true), or S for initial attitude (**GaussMea**==false).
