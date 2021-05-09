# Circles colliding

Circles colliding in a 2D box

### Prerequisites

In order to run this code you need Julia and the following Julia packages:

1. Plots
2. LinearAlgebra

### Explanation

The code is organized into multiple files "circles_in_box_/contact model name/.jl". e.g. "circles_in_box_hertzian_spring_dashpot.jl". It generates N random circles in a square of dimensions _L_ and _H_, initializes their velocities randomly and their acceleration to zero, and simulates their motion. The only modelled forces are contact forces.

The motivation was to implement a basic version of DEM and learn some more Julia. I decided to do it the same way as Cundall in his first paper in the sense that I restricted myself to 2D and circles. In order to not preoccupy myself with force modelling, this code only supports elastic collisions. With the assumption that the forces are normal to the contact planes, the impulses are solved for and the motion is updated:
1. Update motion from Newton's second law, assuming constant acceleration during a time-step.
2. Model forces and obtain "next" acceleration.
    1. Contact forces: Contact detection and contact forces modelling. Contact detection is generally done in two steps, a Broad Collision Check (BCP) and a Narrow Collision Check (NCP). Latest version of this code uses _Sweep And Prune_ (SAP) as BCP.

Although the elements are considered to be rigid solids, the _soft sphere_ approach is very common: up to a small overlap is allowed, from which elastic forces associated to deformation are calculated. This code allows small overlap although deformation forces are not taken into account: collisions are only checked for **after** updating the position.

This code plots in real time and can handle up to 50 circles relatively smoothly. It can also save _gifs_, like this one:

![](30LinearSpringDashpot24fps_spedupncropped.gif)

### Improvements
1. Consider other BCP methods (space oriented like grids, or object oriented which builds a tree structure)
2. Add ellipses as shape.
