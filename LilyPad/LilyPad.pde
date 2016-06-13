/*********************************************************
                  Main Window!

Click the "Run" button to Run the simulation.

Change the geometry, flow conditions, numercial parameters
visualizations and measurments from this window.

This screen has an example. Other examples are found at
the top of each tab. Copy/paste them here to run, but you
can only have one setup & run at a time.

*********************************************************/


Ellipsoid run;
CirculationFinder cf;

// --input parameters-------
int n=(int)pow(2, 8);  // number of grid points along a side of domain
float aoa = 30;  // angle of attack (degrees)
float fineness = 5;  // major/minor axis
float d = 0.25;  // diameter at beginning and end of sim (<1)
float Re = 1e4;  // Reynolds number by ellipsoid length
String name = "/Volumes/Macintosh HD/Users/jamesschulmeister/Dropbox (MIT)/2D+T_ellipsoid/D" + str(round(n/10)) + "_AOA" + str(round(aoa));
boolean recording = false;
// -------------------------

void setup(){
  run = new Ellipsoid(aoa, fineness, d, Re, name, recording);
  size(600, 600);
}

void draw(){
  run.update();
}
