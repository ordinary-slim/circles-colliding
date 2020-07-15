# Circles colliding

Circles doing ellastic collisions in a 2D box

### Prerequisites

In order to run this code you need Julia and the following Julia packages:

1. Plots
2. LinearAlgebra

### Explanation

The main function of the code is _many_circles_in_box_. It generates N random circles in a square of dimensions _L_ and _H_, initializes their velocity randomly and their acceleration to zero, and simulates their motion. The only modelled forces are contact forces.

The motivation was to implement a basic version of DEM and learn some more Julia. I decided to do it the same way as Cundall in his first paper in the sense that I restricted myself to 2D and circles. In order to not preoccupy myself with force modelling, this code only supports elastic collisions. With the assumption that the forces are normal to the contact planes, the impulses are solved for and the motion is updated:
1. Update motion from Newton's second law, assuming constant acceleration during a time-step.
2. Model forces and obtain "next" acceleration.
    1. Contact forces: Contact detection and contact forces modelling.

Although the elements are considered to be rigid solids, the _soft sphere_ approach is very common: up to a small overlap is allowed, from which elastic forces associated to deformation are calculated. This code allows small overlap although deformation forces are not taken into account: collisions are only checked for **after** updating the position.

This code plots in real time and can handle up to 25 circles relatively smoothly. It can also save _gifs_, like this one:

![](24FPS_example.gif)

### Improvements
1. The main aspect that can be improved on is the collision check: right now, they are done without taking into account much information. The common approach is to set a grid with a grid size sufficiently big, and only test for neighbouring grid cells.
2. Add ellipses as shape.
3. Add rotational DOF
