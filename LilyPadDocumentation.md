# Introduction #

Lily Pad is a Computational Fluid Dynamics (CFD) solver for fluid-structure interactions with a set of built-in interactive visualization tools. The goal of Lily Pad is to lower the barrier to CFD by adopting simple high-speed methods, and giving immediate visual feed-back to the user.

It would be imposible to genuinely achieve this goal unless Lily Pad was genuine CFD solver. Therefore the full two-dimensional Navier-Stokes equations are solved and the exact body boundary conditions are applied. However, most of the complications plaguing CFD are avoided by using the Boundary Data Immersion Method (BDIM) to immerse solid bodies into the fluid domain. This allows simulations to be set up and run extremely efficiently. Despite this simplicity (or perhaps because of it), BDIM is very accurate, it has some nice analytic properties, and has been extensively validated.

This documentation covers the [numerical algorithms](LilyPadDocumentation#Methodology.md) used in Lily Pad and then goes over the [code organization](LilyPadDocumentation#Code_Organization.md). Please scan the [bibliography](LilyPadDocumentation#Bibliography.md) for the references with the numerical details and interesting applications. To use Lily Pad, please see the GettingStarted and [Nondimensionalization](https://code.google.com/p/lily-pad/wiki/Nondimensionalization) page.

# Numerical Methodology #

The Boundary Data Immersion Method `[4]` is implemented in Lily Pad to perform simulations of fluid-structure interactions. The basic idea of BDIM (and all Immersed Boundary methods) is to adjust the equations of motion to account for the interaction of the fluid and solid and solve those equations on a simple numerical grid. In Lily Pad this is a uniform Cartesian grid. There is also a simplified version of the conservative Volume of Fluid (cVOF `[3]`) scheme implemented but this is in development.

Lily Pad uses the implicit Large Eddy Simulation (ILES) for subgrid scale modeling, which avoids all of the standard turbulence model selection and parameterization issues. ILES is a validated research method, but as with all numerical modeling of fluid flow, the Reynolds number set in the code will not be reflected in the results if the flow is under-resolved. Additionally, you should keep in mind that there is no such thing as physical two-dimensional turbulence.

The pressure equation is solved using a very efficient Multi-Grid method. Either QUICK or Semi-Lagrangian convection schemes can be used (SL enables much larger time steps but reduces the accuracy near the body). Advancement in time is done using a two-step explicit scheme. If both steps are done, this is second-order accurate using either the QUICK or SL schemes. Single precision floats are used throughout.

BDIM imposes no-slip boundary conditions on the solid/fluid interface. The far field boundary conditions are adjustable, but the default conditions are reflection conditions on the top and bottom plane, a uniform flow inlet condition, and a zero gradient exit condition. The numerical details of the simulation method follow those used in `[1-4]` but the upshot is that Lily Pad is a serious two-dimensional solver.  The algorithms in the solver are state-of-the-art and the code is being continually updated by researchers and capabilities added.

# Code Organization #

Lily Pad is written in [Processing](http://www.processing.org) in an [object-oriented](http://en.wikipedia.org/wiki/Object-oriented_programming) programming style. The code is organized into a set of files, each of which defines a class of objects such as **Body** and **VectorField**. While there are many files in Lily Pad you shouldn't get overwhelmed: only a few are fundamental to the solver, a few handle the visualizations, and the majority are special extension classes and test cases.

## Solver Classes ##

The core classes for the solver are **Field**, **VectorField**, **Body**, and **BDIM**. As with all the classes, there is a simple example code snippet at the top of each file which can be used for testing.

  * **Field** contains the data structures and routines for handling scalar fields such as a pressure or distance field. This includes routines to; initialize, add, and multiply fields; return the gradient or laplacian of a field; interpolate the field value at an input point; as well as setting boundary conditions and advecting the field given a velocity field.

  * **VectorField** contains the data structures and routines for handling vector fields such as the fluid velocity. The vector is made up of two shifted **Field** objects, and most of those routines are reused. Additional routines; take the divergence, curl, and inner product of the vector field; apply the advection/diffusion equations; and project a velocity field to be divergence-free, returning the required pressure field. The **PoissonMatrix** and **MGSolver** (Multi-Grid) classes set-up and invert the matrices for this projection.

  * **Body** is the parent class of all the solid geometries which are used in Lily Pad. This general geometry can be created, translated, rotated, or otherwise updated using prescribed or free-motion or interactively using the mouse. The **OrthoNormal** class is used in body query routines (such as distance?, velocity?) for use in the solver.

  * **BDIM** contains the data and routines to set-up and solve the general solid-fluid interaction problem. The grid size, time step, **Body**, initial velocity field, and fluid viscosity can be set, as well as a flag to use QUICK or Semi-Lagrangian advection. Routines update() and update2() integrate the equations in time, and adaptive time stepping is used if the time step is set to zero.

## Visualization and IO Classes ##

Most classes have a `Class`.display() routine which can be used to visualize an object. This is sufficient for **Body** and **VectorField** objects. **FloodPlot** makes nice customizable flood plots to extend Field.display().
The **Window** class maps from the numerical location in the grid to the screen pixels. Future visualization classes are planned to take advantage of the **Particle** class.

The basic data IO classes are **ReadData** and **SaveData** which read and save text file data line by line. Examples of these routines and more advanced ones such as **SaveVectorField** are found in the code examples on the top of each file.

## Extensions and Test Classes ##

The remainder of the Classes are simply extensions of the classes above and test cases:

  * **SharpField** is an extension of **Field** when there is a sharp discontinuity using the cVOF method. This is used in the **FreeInterface** and **TwoPhase** classes which are under development.

  * **EllipseBody**, **CircleBody**, **NACA**, **FlexNACA**, are extensions of the general **Body** class for specific shapes. These help initialize objects, and may have specialized versions of the display and query routines.

  * **BodyUnion** enables complex bodies to be created by taking the union (as in [set theory](http://en.wikipedia.org/wiki/Union_(set_theory))) of two simple bodies. For instance, a multi-section foil, a foil passing a circle, spinning circles embedded in a larger circular structure, movable flaps on a foil section, etc.

Special test classes such as **AudreyTest** and **InlineFoilTest** are used to set up and run non-trivial tests. These include multi-object tests, tests with complex prescribed or reactive motion, and controllers. Trying out these tests is a great way to see what Lily Pad is capable of and inspire new projects.

## Bibliography ##

  1. G. D. Weymouth, D. G. Dommermuth, K. Hendrickson, and D. K.-P. Yue. Advancements in cartesian-grid methods for computational ship hydrodynamics. _26th Symposium on Naval Hydrodynamics_, Rome, Italy, 17-22 September 2006.
  1. G. D. Weymouth. Physics and Learning Based Computational Models for Breaking Bow Waves Based on New Boundary Immersion Approaches, _MIT PhD Thesis_, 2008.
  1. G. D. Weymouth and D. K.-P. Yue. Boundary data immersion method for cartesian-grid simulations of fluid-body interaction problems. _Journal of Computational Physics_, 2011.
  1. G. D. Weymouth and D. K.-P. Yue. Conservative volume-of-fluid method for free-surface simulations on cartesian-grids. _Journal of Computational Physics_, 229(8):2853 – 2865, 2010.
  1. G. D. Weymouth and M. Triantafyllou. Global vorticity shedding for a shrinking cylinder. _Journal of Fluid Mechanics_, 702(1):470–487, 2012.
  1. J. Izraelevitz, G. D. Weymouth, and M. Triantafyllou. Inline motion in flapping foils for improved force vectoring performance. _Bulletin of the American Physical Society_, 57, 2012.
  1. A. Maertens, G. D. Weymouth, and M. Triantafyllou. Limits of the potential flow model for obstacle detection using a lateral line. _Bulletin of the American Physical Society_, 57, 2012.
  1. J. Schulmeister, J. Dahl, G. D. Weymouth, and M. Triantafyllou. Flow separation control with rotating cylinders. _Bulletin of the American Physical Society_, 57, 2012.
  1. G. D. Weymouth and M. Triantafyllou. Viscous flow around a rapidly collapsing cylinder as a model of animal locomotion. _Bulletin of the American Physical Society_, 57, 2012.
  1. G. D. Weymouth and M. Triantafyllou. Ultra-fast escape of a deformable jet-propelled body. _Journal of Fluid Mechanics_, 721:367–385, 2013.