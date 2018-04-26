/*! \file section.h
    \brief Definition of Poincare sections.

    $Author: roldan $
    $Date: 2013-03-11 11:36:04 $
*/

#ifndef SECTION_H_INCLUDED
#define SECTION_H_INCLUDED

extern const double TWOPI;

/// Different Poincare sections in Cartesian coordinates.
/// SEC1 is {y=0, p_y>0}, SEC2 is {y=0, p_y<0}.
/// 
/// We reuse this type to distinguish different Poincare sections
/// in Delaunay coordinates. 
/// SEC1 is {l=0}, SEC2 is {l=\pi}.
///
/// The 3:1 periodic orbit should be a fixed point for the section {g=0}. 
/// Thus we define a new section SECg that corresponds to {g=0}.
//
/// The 3:1 periodic orbit should be a fixed point for the section {g=\pi}. 
/// Thus we define a new section SECg2 that corresponds to {g=\pi}.
typedef enum {SEC1, SEC2, SECg, SECg2} section_t;    

#endif // SECTION_H_INCLUDED
