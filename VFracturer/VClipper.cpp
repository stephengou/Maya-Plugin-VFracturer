/** 
**  by Zhen Gou, Zhenghan Mei  March 2014
**/
#include "VClipper.h"

bool VClipper:: loadDecompositionCell(Polyhedron_3 & cell)
{
	if(!cell.is_closed()) 
	{
		cout<<"INVALID Decomp CELL: cell mesh is not closed!"<<endl;
		return false;
	}
		
	CGAL::Nef_polyhedron_3<K> nef_poly(cell);
	m_decomp.push_back(nef_poly);
	return true;
}

bool VClipper::loadPatternCell(Polyhedron_3 & cell)
{
	if(!cell.is_closed()) 
	{
		cout<<"INVALID pattern CELL: cell mesh is not closed!"<<endl;
		return false;
	}
		
	CGAL::Nef_polyhedron_3<K> nef_poly(cell);
	m_pattern.push_back(nef_poly);
	return true;
}

bool VClipper::loadImpactCell(Polyhedron_3 & cell)
{
	if(!cell.is_closed()) 
	{
		cout<<"INVALID impact CELL: cell mesh is not closed!"<<endl;
		return false;
	}

	CGAL::Nef_polyhedron_3<K> nef_poly(cell);
	m_impact_cell=nef_poly;
	return true;
}

void VClipper::decomp(vector<Polyhedron_3*> & list)
{
	//calculating the broken parts
	for(int i=0;i<m_decomp.size();++i)
	{
	  cout<<"calculating intersection..."<<endl;
	  CGAL::Nef_polyhedron_3<K> intersect=m_model.intersection(m_decomp.at(i));
	  if(intersect.is_simple()) 
	  {
		  Polyhedron_3 * poly_intersection=new Polyhedron_3();
		  intersect.convert_to_polyhedron(*poly_intersection);
		  list.push_back(poly_intersection);
		  cout<<"intersection completed..."<<endl;
	  }
	  else cout<<"INVALID INTERSECTION: nef_polyhedron is not simple! hence totally removed to avoid exception!"<<endl;
     
	}

}


void VClipper::output(vector<Polyhedron_3*> & list)
{
	//calculating decomposition of whole model
	for(int i=0;i<m_decomp.size();++i)
	{
	  cout<<"calculating intersection..."<<endl;
	  CGAL::Nef_polyhedron_3<K> intersect=m_model.intersection(m_decomp.at(i));
	  if(intersect.is_simple()) 
	  {
		  Polyhedron_3 * poly_intersection=new Polyhedron_3();
		  intersect.convert_to_polyhedron(*poly_intersection);
		  list.push_back(poly_intersection);
		  cout<<"intersection completed..."<<endl;
	  }
	  else cout<<"INVALID INTERSECTION: nef_polyhedron is not simple! hence totally removed to avoid exception!"<<endl;
     
	}

	//decomposing impact cell using pattern
	CGAL::Nef_polyhedron_3<K> impact=m_model.intersection(m_impact_cell);
	for(int i=0;i<m_pattern.size();++i)
	{
		cout<<"calculating intersection..."<<endl;
		CGAL::Nef_polyhedron_3<K> intersect=impact.intersection(m_pattern.at(i));
		if(intersect.is_simple()) 
	    {
		    Polyhedron_3 * poly_intersection=new Polyhedron_3();
		    intersect.convert_to_polyhedron(*poly_intersection);
		    list.push_back(poly_intersection);
		    cout<<"intersection completed..."<<endl;
	    }
	    else cout<<"INVALID INTERSECTION: nef_polyhedron is not simple! hence totally removed to avoid exception!"<<endl;
	}


}

VClipper::VClipper( Polyhedron_3 & model):m_model(model){}
