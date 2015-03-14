# Introduction #

In [Lily Pad](LilyPadDocumentation.md), the solid and fluid mechanical governing equations are used in their non-dimensional form. Because the dimensions (size, speed, and time scales) of a problem are defined by the user, these cannot be "hard-wired" into the governing equations. Instead, the grid size _h_ and uniform flow velocity _U_ are used. This makes certain [grid-based parameters](Nondimensionalization#Numerical_Implications.md) easy to compute, but means any [engineering coefficients](Nondimensionalization#Engineering_Implications.md) will require some scaling.

This document reviews the basics of nondimensionalization for the [Navier-Stokes](Nondimensionalization#Navier_Stokes_Equation.md) equation, and gives a few [examples](Nondimensionalization#Example.md). For a more detailed discussion of nondimensionalization in physics, you can check out wikipedia (  [general](http://en.wikipedia.org/wiki/Nondimensionalization) & [Navier-Stokes](http://en.wikipedia.org/wiki/Non-dimensionalization_and_scaling_of_the_Navier%E2%80%93Stokes_equations) ) or classics texts such as [Bachelor 1967](http://books.google.co.uk/books?id=Rla7OihRvUgC&lpg=PP1&pg=PP1#v=onepage&q&f=false) or [Newman 1977](http://books.google.co.uk/books?id=nj-k_lAmaBYC&lpg=PP1&dq=newman%201977&pg=PP1#v=onepage&q=newman%201977&f=false).

# Navier Stokes Equation #

The Navier-Stokes equations governing the dynamics of a simple fluid with uniform material properties is:

> _rho D/Dt(u) = -d/dx(p)+mu d<sup>2</sup>/dx<sup>2</sup>(u)_

where _u_ is the local fluid velocity vector, and _p_ is the local pressure of the fluid. The properties _rho_ and _mu_ are the fluid density and viscosity. The operators _D/Dt_, _d/dx_, and _d<sup>2</sup>/dx<sup>2</sup>_ are the material derivative, gradient and laplacian, respectively.

These variables have dimensions of mass, length, and time, and need to be scaled by characteristic values in the flow. As discussed above, the values typically chosen by engineers (such as the length or mass of a body) are defined by the user, and therefore cannot be hard-wired into the solver. Additionally, because Lily Pad allows for multi-body problems there may be any number of length and velocity scales to choose from for a given simulations.

Instead the code is scaled by the fixed numerical parameters; the grid cell dimension _h_, the uniform flow velocity _U_, and the fluid density _rho_. The non-dimensional variables are then:

  * _**u** = u / U_
  * _**p** = p / ( rho U<sup>2</sup> )_
  * _**x** = x / h_
  * _**t** = t U / h_

where _**x**_ and _**t**_ are non-dimensional distances and times used in operators. The Navier-Stokes equation is then:

> _**D/Dt**(**u**) = -**d/dx**(**p**)+**nu** **d<sup>2</sup>/dx<sup>2</sup>**(**u**)_

where the non-dimensional kinematic viscosity is

  * _**nu** = mu / ( rho U h )_


# Example #

Because the equations have been nondimensionalized using _h_ and _U_, it is most appropriate to use these variables when setting up the engineering problem of interest.

For example, consider simulating the flow around a cylinder moving with imposed harmonic motion. The size and location of the body, and amplitude of the motion is expressed in units of _h_. The frequency is expressed in units of _U/h_ and the velocity in units of _U_. The code to place this body is then:
```
CircleBody body;
float t=0;                // initialize non-dim time at zero
void setup(){
  size(400,400);          // display size in pixels

  int n=64+2;             // number of grid cells in each direction
  float x  = 32, y = 32;  // circle center is located at grid center
  float diam = 20;        // diameter is 20h
  
  // define the geometry (the `Window' sizes the circle to the display) 
  body = new CircleBody( x, y, diam, new Window(n,n) );
}
void draw(){
  t += 0.1;               // increment time by 0.1 h/U
  float amp = 20, cen=32; // amplitude of motion is 20h
  float omega = PI/6.;    // frequency is pi/6 U/h

  // displacement is desired location - current location
  float dy = cen+amp*sin(omega*t)-body.xc.y;

  // apply translation
  body.translate(0,dy);   
  body.update();

  // draw background and current position of body
  background(0);
  body.display();

  // write the velocity to screen (displacement/time per step)
  println("velocity: "+body.dxc.y/0.1);
}
```

Copy and paste this code into Lily Pad to try it out. What is the magnitude of the velocity printed to the output? Does it match your expectations given the defined frequency and amplitude?

# Engineering Implications #

The above discussion and example state that _h_,_U_, and _rho_ are the fundamental scales in Lily Pad. This has implication for working with the code:

  * As in the previous section, the problem of interest should be set up in terms of these units.
  * In addition to size and velocities, the viscosity must also be scaled. To achieve a Reynolds number based on length of, _**Re<sub>L</sub>** = rho U L / mu_, then the non-dimensional kinematic viscosity should be set to _**nu**=**L Re<sub>L</sub>**<sup>-1</sup>_. Where _**L**=L / h_ as above.
  * In general the output will need to be converted into engineering scales as well. For instance to a force _**F**_ from Lily Pad corresponds to a two-dimensional drag coefficient of _**C<sub>D</sub>**=2F / (rho U<sup>2</sup> L ) = 2 **F / L**_.

# Numerical Implications #

However, many numerical parameters of interest are easier to work with in Lily Pad's non-dimensional form:

  * The non-dimensional length _**L** = L/h_ is the body resolution, giving the number of grid points along the body. This is a heuristic for how well resolved the simulation will be.
  * The non-dimensional time step _**dt** = dt U / h_ is the global Courant number, specifying what fraction of a grid cell is traveled by the flow in a time step. This is important for the numerical stability of the method.
  * The non-dimensional viscosity _**nu** = mu / ( rho U h ) = **Re<sub>h</sub>**<sup>-1</sup>_ is the inverse grid-based Reynolds number. This indicates how well resolved the viscous stresses will be. Setting _**nu** << 1_ will under-resolve those forces. Increasing the body resolution  _**L**_ will allow you to match a given _**Re<sub>L</sub>**_ without an excessively small _**nu**_.