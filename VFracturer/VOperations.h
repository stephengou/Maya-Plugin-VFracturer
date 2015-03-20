/** A collection of utility methods for generating and manipulating geometries
**  by Zhen Gou, Zhenghan Mei  2014
**/
#ifndef VOPERATIONS_H
#define VOPERATIONS_H
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_inventor_ostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Vector_3.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "../glm/glm.hpp"
#include "../glm/gtc/matrix_transform.hpp"

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq RT;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::Quotient<CGAL::MP_Float> RT;
#endif

typedef CGAL::Cartesian<RT>                        K;//it was Cartesian<RT>
typedef K::Point_3                                 Point_3;
typedef CGAL::Polyhedron_3< K>                     Polyhedron_3;

typedef CGAL::Convex_hull_d_traits_3<K>            Hull_traits_3;
typedef CGAL::Convex_hull_d< Hull_traits_3 >       Convex_hull_3;
typedef CGAL::Creator_uniform_3<double, Point_3>   Creator;

Polyhedron_3 * sphereConvexHull (int n, int r);
Polyhedron_3 * testOperation();
Polyhedron_3 * intersection(Polyhedron_3 &A,Polyhedron_3& B);
Polyhedron_3 * genCube();
void genRegularVoroCell(std::vector<Polyhedron_3 *> &cells, int num, glm::vec3 max, glm::vec3 min);
void genVoroCell(std::vector<Polyhedron_3 *> &cells, std::vector<Polyhedron_3 *> &otherCells, Polyhedron_3 * &impactCell,int num, glm::vec3 max, glm::vec3 min, glm::vec3 pos, float r, int type);
void genVoroPatternCell(std::vector<Polyhedron_3 *> &cells, int num, glm::vec3 max, glm::vec3 min, glm::vec3 pos,float percentage, int type);
void translate(Polyhedron_3 & poly, float x, float y, float z);
void scale(Polyhedron_3 & poly, float x, float y, float z);

void dimensionPolyhedron_3( Polyhedron_3 & poly, float & minX,  float & maxX,  float & minY,  float & maxY,  float & minZ,  float & maxZ );
bool withinDimension(glm::vec3 pos, glm::vec3 Max, glm::vec3 Min);
void clipDimension(glm::vec3 & pos, glm::vec3 & max, glm::vec3 & min);

#endif