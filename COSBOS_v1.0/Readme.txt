/* ********************************************************************* */

Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
Signal Analysis and Machine Perception Laboratory, 
Department of Electrical, Computer, and Systems Engineering, 
Rensselaer Polytechnic Institute, Troy, NY 12180, USA

You are free to use this software, but we would appreciate it if you can 
cite our papers: 

[1] Quan Wang, Xinchi Zhang, Kim L. Boyer, "Occupancy distribution 
estimation for smart light delivery with perturbation-modulated light 
sensing", Journal of Solid State Lighting 2014 1:17, ISSN 2196-1107, 
doi:10.1186/s40539-014-0017-2. 

[2] Quan Wang, Xinchi Zhang, Meng Wang, Kim L. Boyer, "Learning Room 
Occupancy Patterns from Sparsely Recovered Light Transport Models", 22nd 
International Conference on Pattern Recognition (ICPR), 2014. 

[3] Quan Wang, Xinchi Zhang, Kim L. Boyer, "3D Scene Estimation with
Perturbation-Modulated Light and Distributed Sensors", 10th IEEE Workshop 
on Perception Beyond the Visible Spectrum (PBVS). 

[4] Xinchi Zhang, Quan Wang, Kim L. Boyer, "Illumination Adaptation with 
Rapid-Response Color Sensors", SPIE Optical Engineering + Applications, 2014. 

/* ********************************************************************* */

This package includes: 

1. The code for Light Transport Model (LTM) recovery, both overdetermined 
and underdetermined. See 'LTM_Recovery/demo_LTM.m' for a demo. This work 
is described in [2]. 

2. The code for 3D scene estimation with light blockage model and wall-mounted
sensors. See 'BlockageModel/demo_Blockage.m' for a demo. This work is 
described in [1] and [3]. 

3. The code for floor-plane occupancy mapping with light reflection model 
and ceiling-mounted sensors. See 'ReflectionModel/demo_Reflection.m' for 
a demo. This work is described in [1]. 

* [4] is not directly related to this package. It is about the new color 
sensors that we have built for occupancy sensing. 

* In each demo, we included example data. But the code for data collection 
has too many dependencies: platform, hardware, driver, and other software 
packages. Thus we are not including the code for data collection here. 

* More information on this work can be found here: 
https://sites.google.com/site/cosboswiki/

/* ********************************************************************* */

This work was supported primarily by the Engineering Research Centers 
Program (ERC) of the National Science Foundation under NSF Cooperative 
Agreement No. EEC-0812056 and in part by New York State under NYSTAR 
contract C090145. 

/* ********************************************************************* */


