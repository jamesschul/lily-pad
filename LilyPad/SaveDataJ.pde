/**********************************
 SaveDataJ class

 Saves data to a text file with customizable header--James special:)
***********************************/

class SaveDataJ{
  ArrayList<PVector> coords;
  PrintWriter output;

  SaveDataJ(String name) {
    output = createWriter(name);
  }

  void addFloat(float t, float f) {
    output.print(t + " ");
    output.print(f + " ");
    output.println(";");
  }

  void addPVector(float t,  PVector v) {
    output.print(t + " ");
    output.print(v.x + " ");
    output.print(v.y + " ");
    output.println(";");
  }

  void addProfileData(float t, ArrayList<PVector> coords, Field a){
    output.print(t + " ");
    int n = coords.size();
    for(int i=0; i<n; i++){
      output.print(a.linear( coords.get(i).x, coords.get(i).y ) +" ");
    }
    output.println(";");
  }

  void addField(float t, Field a){ //save the entire field data (square)
    for (int j=1; j<a.n-1; j++) {
    for (int i=1; i<a.n-1; i++) {
      output.print(a.a[i][j] +" ");
    }
    output.println(";");
    }
  }

  void addField_5(float t, Field a){ //save every fifth data point of the field (square)
    for (int j=1; j<a.n/5; j++) {
    for (int i=1; i<a.n/5; i++) {
      output.print(a.a[i*5][j*5] +" ");
    }
    output.println(";");
    }
  }

  void finish(){
    output.flush(); // Writes the remaining data to the file
    output.close(); // Finishes the file
  }
}


/**********************************
Write and Resume BDIM functions

For saving and then restarting where you left off.
***********************************/

void writeBDIM( String name, float time, BDIM flow ) {
  PrintWriter output;
  output = createWriter(name);
  output.println(time);
  output.println(flow.dt);
  int n = flow.p.n;
  for ( int i=0; i<n; i++ ) {
    for ( int j=0; j<n; j++ ) {
      output.println(""+flow.u.x.a[i][j]+", "+flow.u.y.a[i][j]+", "+flow.p.a[i][j]);
    }
  }
  output.flush();
  output.close();
}

float resumeBDIM( String name, BDIM flow) {
  float[] data;
  String[] stuff = loadStrings(name);
  float t = float(stuff[0]);
  flow.dt = float(stuff[1]);
  int n = flow.p.n;
  for ( int i=0; i<n; i++ ) {
    for ( int j=0; j<n; j++ ) {
      data = float(split(stuff[2+i*n+j], ','));
      flow.u.x.a[i][j] = data[0];
      flow.u.y.a[i][j] = data[1];
      flow.p.a[i][j] = data[2];
    }
  }
  return t;
}
