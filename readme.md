# EKFPNP

We compute camera pose parameters from a sequence of images using a sequential estimation procedure.
Taking into account both the camera motion model and the noisy observation model, the results
tend to be more accurate and more robust compared to when the observations are considered alone.
To combine the camera pose history and the uncertainty of image point measurements, the Extended
Kalman Filter (EKF) is employed.

## PnP Toolbox

This toolbox is an extension of the toolbox provided by the authors of
MLPnP, CEPPnP and OPnP
We extended it to show the use of EKFPnP.

## Scripts to produce the synthetic results:
main_ordinary_3d.m        
main_planar.m          
main_ordinary_3d_sigma.m   
main_planar_sigma.m        
main_ordinary_3d_time.m    
main_planar_time.m	   


## Acknowledgements
This toolbox is based on the toolbox provided in:

MLPnP: 
Urban, S.; Leitloff, J.; Hinz, S..
MLPNP - A REAL-TIME MAXIMUM LIKELIHOOD SOLUTION TO THE PERSPECTIVE-N-POINT PROBLEM.
ISPRS Annals of Photogrammetry,
https://github.com/urbste/MLPnP_matlab

CEPPnP:
Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer.
Leveraging Feature Uncertainty in the PnP Problem. 
In Proceedings of BMVC, 2014. 

REPPnP: 
Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer.
Very Fast Solution to the PnP Problem with Algebraic Outlier Rejection. 
In Proceedings of CVPR, 2014.

OPnP: 
Y. Zheng, Y. Kuang, S. Sugimoto, K. Astro ?m, and M. Okutomi. Revisiting the pnp problem: A fast, general and opti-
mal solution. In ICCV, pages 4321?4328, 2013.

ASPnP: 
Y. Zheng, S. Sugimoto, and M. Okutomi. Aspnp: An accurate and scalable solution to the perspective-n-point problem. 
Trans. on Information and Systems, 96(7):1525?1535, 2013.



