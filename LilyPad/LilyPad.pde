/*********************************************************
2D+T Simulation approximations of flow past slender bodies
at an angle of attack.

James Schulmeister
MIT, 2016
*********************************************************/

Ellipsoid run;
CirculationFinder cf;

// --input parameters-------
int p = 8;
int n=(int)pow(2, p);  // number of grid points along a side of domain
float aoa = 10;  // angle of attack (degrees)
float fineness = 5;  // major/minor axis
float d = 0.25;  // diameter at beginning and end of sim (<1)
float Re = 1e4;  // Reynolds number by ellipsoid length
String name = "/Volumes/Macintosh HD/Users/jamesschulmeister/Dropbox (MIT)/2D+T_ellipsoid/simOutput/spheroid/D" + str(p) + "_AOA" + str(round(aoa));
boolean recording = true;
boolean recordingFields = false;
// -------------------------

void setup(){
  run = new Ellipsoid(aoa, fineness, d, Re, name, recording);
  size(600, 600);
}

void draw(){
  run.update();
}
