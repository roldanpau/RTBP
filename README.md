# RTBP
Restricted Three Body Problem library

##Introduction
  This is a library of functions for different computations related to the
  Restricted Three Body Problem, or RTBP for short. The software is available
  upon request (just email me), so that I can keep track of who's using it
  and why. 
  This software has been developed for my own research, in particular for the
  following paper:
  "Diffusion along mean motion resonance in the restricted planar three-body
  problem" by J.Fejoz, M.Guardia, V.Kaloshin and P.Roldan.
  (<a href="http://arxiv.org/abs/1109.2892">link</a>)
  The kind of functions that you will find here are
  - Numerical integration of the RTBP equations (in Cartesian variables)
  - Numerical integration of the RTBP equations (in Delaunay variables)
  - Cartesian to Delaunay change of coordinates
  - Numerical integration of the RTBP variational equations 
  - Poincare section and associated Poincare map (in Cartesian)
  - Poincare section and associated Poincare map (in Delaunay)
  - Numerical computation of (some family of) periodic orbits
  - Hyperbolic splitting (eigenvalues/vects) associated to a fixed point of
  the Poincare map
  - Numerical computation of hyperbolic invariant manifolds associated to
  fixed point
  - Homoclinic intersection of invariant manifolds
  - Numerical computation of splitting angle of the manifolds at the
  homoclinic intersection
  
##Acknowledgements
  A big thanks goes to Angel Jorba for altruistically sharing his ideas about
  RTBP computations with me. Moreover, he provided me with the high-precision 
  Taylor method integrator which is at the base of some of my computations.
  See 
  <a href="http://www.maia.ub.edu/~angel/soft.html">http://www.maia.ub.edu/~angel/soft.html</a>.
