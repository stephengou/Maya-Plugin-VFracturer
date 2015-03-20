/** a clipper class that breaks the input model according to pattern and output a list of models.
**  by Zhen Gou, Zhenghan Mei  March 2014
**/
#ifndef VCLIPPER_H
#define VCLIPPER_H
#include "VOperations.h"
#include "VConverter.h"
//model and pattern cells are assumed to be transformed to desired positions and scales already
class VClipper
{
public:
	VClipper(Polyhedron_3 & model);
	bool loadDecompositionCell(Polyhedron_3 & cell); //EXCEPT THE IMPACT CELL
	bool loadImpactCell(Polyhedron_3 & cell);
	bool loadPatternCell(Polyhedron_3 & cell);
	void output(vector<Polyhedron_3*> & list);
	void decomp(vector<Polyhedron_3*> & list);//for debug use only
private:
	CGAL::Nef_polyhedron_3<K> m_model;
	vector<CGAL::Nef_polyhedron_3<K>> m_decomp;
	vector<CGAL::Nef_polyhedron_3<K>> m_pattern;
	CGAL::Nef_polyhedron_3<K> m_impact_cell;

	VClipper(); //do not allow default construction
};
#endif