/*! \file section.h
    \brief Definition of Poincare sections.

    $Author: roldan $
    $Date: 2013-03-11 11:36:04 $
*/

#ifndef SECTION_H_INCLUDED
#define SECTION_H_INCLUDED

/// Different Poincare sections in Cartesian coordinates.
/// SEC1 is {y=0, p_y>0}, SEC2 is {y=0, p_y<0}.
/// 
/// We reuse this type to distinguish different Poincare sections
/// in Delaunay coordinates. 
/// SEC1 is {l=0}, SEC2 is {l=\pi}.
typedef enum {SEC1, SEC2} section_t;    

#endif // SECTION_H_INCLUDED
