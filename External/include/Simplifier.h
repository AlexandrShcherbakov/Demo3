#ifndef DEMO3_SIMPLIFIER
#define DEMO3_SIMPLIFIER

#include <vector>
#include <utility>
#include <tuple>
#include <algorithm>
#include <set>
#include <map>
#include <list>
#include <cfloat>
#include <iostream>

#include "VectorMath.h"
#include "Geometry.h"

///Convenient typedefs for long typenames
typedef std::tuple<std::vector<vec3>, std::vector<uint> > geometry;

///Convenient small functions
geometry make_geometry(const std::vector<vec3> points, const std::vector<uint> indices);

///*****************************************************************************
///Geometry class
///*****************************************************************************



///*****************************************************************************
///QuadricErrorMetricsSimplifier
///*****************************************************************************

geometry QuadricErrorMetricsSimplifier(const geometry& baseGeometry, const uint polygonLimit);
Geometry NewQuadricErrorMetricsSimplifier(const Geometry baseGeometry, const uint polygonLimit);

///*****************************************************************************
///SOMmapSimplifier
///*****************************************************************************

geometry SOMmapSimplifier(const geometry& baseGeometry, const uint polygonLimit);


///*****************************************************************************
///SmallEdgesSimplifier
///*****************************************************************************

Geometry SmallEdgesSimplifier(const Geometry &inGeometry, const uint polygonLimit);

#endif
