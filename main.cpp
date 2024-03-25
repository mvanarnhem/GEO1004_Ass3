#include <iostream>
#include <fstream>
#include <string>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Plane_3.h>
#include "../include/json.hpp"
#include <cstdlib>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Cartesian<double>     K;
typedef K::Point_2                  Point_polygon;
typedef CGAL::Polygon_2<K>          Polygon;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3                      Plane;
typedef Kernel::Point_3                      Point_3;