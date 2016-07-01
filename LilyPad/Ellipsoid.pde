/*************************
Ellipsoid Class
Circle that expands and shrinks, modeling the 2D+T flow past an ellipsoid

Example code:

Ellipsoid run;

// --input parameters-------
int n=(int)pow(2, 9)+2;  // number of grid points along a side of domain
float aoa = 15;  // angle of attack (degrees)
float fineness = 5;  // major/minor axis
float d = 0.25;  // diameter at beginning and end of sim (<1)
float Re = 1e4;  // Reynolds number by ellipsoid length
String name = "AOA_"+str(aoa);
boolean saveData = false;
// -------------------------

void setup(){
  run = new Ellipsoid(aoa, fineness, d, Re, name, saveData);
}

void draw(){
  run.update();
}

***********************/

class Ellipsoid{

  ExpandingBDIM flow;
  EllipsoidBody body;
  FloodPlot flood;
  float t=0, D, aoa, d, L_D;
  SaveDataJ u;
  SaveDataJ v;
  SaveDataJ w;
  SaveDataJ p;
  SaveDataJ psi;
  SaveDataJ pProfile;
  SaveDataJ drag;
  SaveDataJ dist;
  SaveDataJ diameter;
  int uin;

  Ellipsoid(float a, float fine, float ds, float Re, String name, boolean recording){
    this.D = n/10.;  // blockage ratio of 10%
    this.L_D = fine;
    this.d = ds;
    this.aoa = a*PI/180;  // convert angle of attack from degrees to radians

    // set window view area and zoom
    float zoom = 3;
    Window view = new Window(int((n-n/zoom)/2), int((n-n/zoom)/2), int(n/zoom), int(n/zoom));

    body = new EllipsoidBody(n/2, n/2, d*D, 1, view); // define geom

    this.uin = 1;  // free stream
    flow = new ExpandingBDIM(n, n, float(0), body, L_D*D/Re, true, uin);   // QUICK with adaptive dt
    flow.dt = .01; // initial time step

    cf = new CirculationFinder(flow,body,view);
    cf.setAnnotate(true,1.0/d);

    flood = new FloodPlot(view);
    flood.range = new Scale(-.5, .5);
    flood.setLegend("vorticity");

    // initialize output files
    if (recording) {
      pProfile = new SaveDataJ(name + "_pProfile.txt");
      drag = new SaveDataJ(name + "_drag.txt");
      diameter = new SaveDataJ(name + "_diameter.txt");
    }
    if (recordingFields) {
      u = new SaveDataJ(name + "_u.txt");
      v = new SaveDataJ(name + "_v.txt");
      w = new SaveDataJ(name + "_w.txt");
      p = new SaveDataJ(name + "_p.txt");
      psi = new SaveDataJ(name + "_psi.txt");
      dist = new SaveDataJ(name + "_dist.txt");
    }
  }

  void update() {
    float dt = flow.dt;
    t += dt;
    float aspect = 1;

    float trans = .5*d;  // time to slow body to rest
    if (t/D < trans) { // smoothly slow body to rest
      body.translate(.5*dt*(cos(PI*(t/D)/trans)+1), 0);
    }

    // increase and decrease smoothly
    float T = L_D*tan(aoa);
    float Ds = d;
    float De = sqrt(1 - sq(t/D - T/2 - trans)/sq(T/2));
    float mergeTime = 0.2*sin(aoa);
    float DedtMerge = -4*(trans+mergeTime - T/2 - trans)/(sq(T)*sqrt(1 - 4*sq(trans+mergeTime - T/2 - trans)/sq(T)));
    float beta = atan(1/DedtMerge);
    float Rmerge = (sqrt(1 - sq(trans+mergeTime - T/2 - trans)/sq(T/2)) - d)/(1 - sin(beta));
    float k = Rmerge + d;
    float h = trans + mergeTime - Rmerge*cos(beta);
    float hEnd = trans + T - mergeTime + Rmerge*cos(beta);
    float circ = k - sqrt(sq(Rmerge) - sq(t/D - h));
    float circEnd = k - sqrt(sq(Rmerge) - sq(t/D - hEnd));

    // Merging section
    if ((t/D > h) && (t/D < (trans+mergeTime))) {
      Ds = circ;
    }

    // Ellipse section
    if((t/D > trans+mergeTime) && (t/D < trans + T - mergeTime)){
      Ds = De;
    }

    // End Merging section
    if((t/D > trans + T - mergeTime) && (t/D < trans + T - mergeTime + Rmerge*cos(beta))){
      Ds = circEnd;
    }

    // Final Sting section
    if(t/D > trans + T - mergeTime + Rmerge*cos(beta)){
      Ds = d;
    }

    body.update(Ds*D,aspect);  // pass the new diameter and aspect ratio to the body
    flow.update(body);
    flow.update2();

    flood.display(flow.u.vorticity());
    body.display();
    flood.displayTime(t/D);

    if (recording) {
      PVector force = body.pressForce(flow.p);
      drag.addFloat(t/D,force.x);
      diameter.addFloat(t/D,Ds);
    }

    if (((t)%(0.1*D)) <= dt) {
      // add data to output files
      if (recording) {
        pProfile.addProfileData(t/D,body.coords,flow.p);
      }
      if (recordingFields) {
        u.addField(t/D, flow.u.x);
        v.addField(t/D, flow.u.y);
        w.addField(t/D, flow.u.vorticity());
        p.addField(t/D, flow.p);
        psi.addField(t/D, flow.u.streamFunc());
        dist.addField(t/D, flow.bodyNull);
      }
      if (t/D > 1.25*(trans+T)) {
        // close the output data files
        if (recording) {
          pProfile.finish();
          drag.finish();
          diameter.finish();
        }
        if (recordingFields) {
          u.finish();
          v.finish();
          w.finish();
          p.finish();
          psi.finish();
          dist.finish();
        }
        cf.update();
        cf.display();
        // saveFrame("AOA20_ellipse.png");
        exit();
      }
    }
  }
}

// EllipsoidBody is a circle that expands and contracts according to the
// 2D+T approximation for the flow past an ellipsoid at and angle of attack.
class EllipsoidBody extends EllipseBody {
    float dh = 0;
    PVector dcen[];

    EllipsoidBody( float x, float y, float _h, float _a, Window window ) {
      super(x, y, _h, _a, window);
      dcen = new PVector[m];
      for ( int i=0; i<m; i++ ) dcen[i] = new PVector(0, 0);
    }

    void update(float d, float as) {
      // update the circle's diameter
      h = d;
      a = as;

      // update body coordinates
      EllipseBody nbod = new EllipseBody(xc.x, xc.y, h, a, window);
      for ( int i=0; i<m; i++ ) {
        dcen[i] = PVector.sub(nbod.orth[i].cen, orth[i].cen);
      }
      coords = nbod.coords;
      end();
    }

    // get "nearest" section body velocity for collapse
    float velocity( int d, float dt, float x, float y ) {
      float dis = -1e10;
      PVector v = new PVector(0, 0);
      for ( int i=0; i<m; i++ ) {
        float d2 = orth[i].distance(x, y);
        if ( d2 > dis ) {
          dis = d2;
          v = dcen[i];
        }
      }

      PVector r = new PVector(x, y);
      r.sub(xc);
      if (d==1) return dxc.x/dt - r.y*dphi/dt + v.x/dt;
      else     return dxc.y/dt + r.x*dphi/dt + v.y/dt;
    }

    // calculate flux
    float get_flux ( VectorField p ) {
      float pv = 0;
      for ( OrthoNormal o: orth ) {
        float pdlx = p.x.linear( o.cen.x, o.cen.y )*o.nx*o.l;
        float pdly = p.y.linear( o.cen.x, o.cen.y )*o.ny*o.l;
        pv += pdlx + pdly;
      }
      return pv;
    }

     void display( color C, Window window ) { // note: can display while adding
       noFill();
       stroke(bodyOutline);
       strokeWeight(1);
       beginShape();
       for ( PVector x: coords ) vertex(window.px(x.x), window.py(x.y));
       endShape(CLOSE);
     }
  }

class ExpandingBDIM extends BDIM {
  float flux;
  Field delc;
  Field bodyNull;

  ExpandingBDIM( int n, int m, float dt, EllipsoidBody body, float nu, boolean QUICK) {
    super(n, m, dt, body, nu, QUICK);
  }

  ExpandingBDIM( int n, int m, float dt, EllipsoidBody body, float nu, boolean QUICK, int uInf) {
    super(n, m, dt, body, nu, QUICK, uInf);
  }

  void update( Body body ) {
    get_coeffs(body);
    update();
  }

  void updateUP( VectorField R, VectorField coeff ) {
    run.flow.u.x.gradientExit = false;
    R.plusEq(PVector.mult(g, dt));
    u.eq(del.times(R).minus(ub.times(del.plus(-1))));
    if (mu1) u.plusEq(del1.times((R.minus(ub)).normalGrad(wnx, wny)));
    u.setBC();

    // update cell center delta
    delc = new Field(n, m, 0, 1);
    for ( int i=1 ; i<n-1 ; i++ ) {
      for ( int j=1 ; j<m-1; j++ ) {
        float dis = run.body.distance( (float)(i), (float)(j) );
        delc.a[i][j] = delta0(dis);
      }
    }

    // update bodyNull field to remove extraneous vorticity
    bodyNull = new Field(n, m, 0, 1);
    for ( int i=1 ; i<n-1 ; i++ ) {
      for ( int j=1 ; j<m-1; j++ ) {
        float dis = run.body.distance( (float)(i), (float)(j) );
        bodyNull.a[i][j] = delta0(dis - 2*eps);
      }
    }

    // new exit treatment
    float sink = run.body.get_flux(run.flow.ub)/(2*PI);

    // downstream
    for ( int j=0 ; j<m ; j++ ) {
      float x = n/2, y = j - m/2;
      float r = sqrt(sq(x) + sq(y));
      float theta = atan(abs(y/x));
      // u.x.a[n-1][j] += sink/r*cos(theta);
      u.x.a[n-1][j] += sink/m*2*PI;
    }

    // new source term
    Field s = u.divergence().plus(delc.plus(-1).times(ub.divergence()));
    p = u.project(coeff, p, s);
  }
}
