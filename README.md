# Time-varying frequency Estimation of narrow-band signals via Extended Kalman Filter and Unscented Kalman Filter

## Files list
* **main_audio.m** Performs frequency tracking of 3 audio files at 3 different frequencies. Generates a plot of the estimation. 
* **main_pi.m** Performs analysis of the Performance Index and generates plots and data in the `data/pi` folder. Iterates over several values of `q` and `sigma_noise`
* **main_ri.m** performs analysis of the Robustness Index and generates plots and data in the `data/ri` folder.

## Flags
### main_pi.m
* **LOGSPACE** determines if the array of sigmas is spaced linearly (`false`) or logarithmically (`true`)
### ekf/ekf.m
* **SAVE_EKF_PLOT** if true, generates a video in `generated_video/` of the estimation process.
* **LOG** if true, generates a complete log of the values during the estimation
* **VERBOSE** if true, prints to console the iteration number. Useful in conjunction with `SAVE_EKF_PLOT` to check on the progress, which makes the process really slow.
### ukf/ukf.m
* **VERBOSE** if true, prints to console the iteration number. Useful in conjunction with `SAVE_EKF_PLOT` to check on the progress, which makes the process really slow.
* **LOG** if true, generates a complete log of the values during the estimation
* **SAVE_UKF_PLOT** if true, generates a video in `generated_video/` of the estimation process.
* **SAVE_SP_PLOT** if true, generates a video in `generated_video/` showing the motion of sigma points in the state space.
### pi/pi_plot.m
* **PLOT_FIT** if false, plots raw data cloud. If true, plots the Least Square fit of a hyperbola of the data
### main_ri.m
* **PLOT_PREDICTION** if true, generates a video in `generated_video/` with an iteration for every frame
* **PLOT_PROFILE** if true, opens a plot with the frequency profile of the generated signal
* **COMPUTE** if false, just plots without computing, loading previously generated data.
