/**
**  by Zhen Gou, Zhenghan Mei  March 2014
**/
#include "VOperations.h"
#include "VConverter.h"
#include "voro++.hh"
using namespace voro;

double rnd() {return double(rand())/RAND_MAX;}
Polyhedron_3 * sphereConvexHull (int n, int r)
{
	 
	  Convex_hull_3 hullA(3); 
	  CGAL::Random_points_in_sphere_3<Point_3, Creator> gen(r);

	  for (int i = 0; i < n ; i++, ++gen)
		 hullA.insert(*gen);

	  assert(hullA.is_valid());

	  Polyhedron_3 * P=new Polyhedron_3(); 
	  CGAL::convex_hull_d_to_polyhedron_3(hullA,*P);
	  std::cout << "The convex hull has " << P->size_of_vertices() 	<< " vertices" << std::endl;

	  return P;
}
//
Polyhedron_3 * intersection(Polyhedron_3 &A,Polyhedron_3& B)
{  
	  if(!A.is_closed()) cout<<"A is no valid because it's not closed!"<<endl;
	  if(!B.is_closed()) cout<<"B is no valid because it's not closed!"<<endl;
	  Polyhedron_3 * poly_intersection=new Polyhedron_3();
	  CGAL::Nef_polyhedron_3<K> NA(A);
	  CGAL::Nef_polyhedron_3<K> NB(B);
	  cout<<"calculating intersection..."<<endl;
	  CGAL::Nef_polyhedron_3<K> intersect=NA.intersection(NB);
	  if(intersect.is_simple()) intersect.convert_to_polyhedron(*poly_intersection);
	  else cout<<"nef_polyhedron is not simple! hence totally removed to avoid exception!"<<endl;
      cout<<"intersection completed..."<<endl;
	  return poly_intersection;
}

Polyhedron_3 * genCube()
{
	Polyhedron_3 * cube=new Polyhedron_3();
	Mesh m;
	m.addVertex(-1,1,-1);m.addVertex(1,1,-1);m.addVertex(1,1,1);m.addVertex(-1,1,1);m.addVertex(-1,-1,-1);m.addVertex(1,-1,-1);m.addVertex(1,-1,1);m.addVertex(-1,-1,1);
	m.init();
	vector<int> f;
	f.push_back(0);f.push_back(1);f.push_back(2);f.push_back(3); m.addFace(f);f.clear();
	f.push_back(4);f.push_back(7);f.push_back(6);f.push_back(5); m.addFace(f);f.clear();
	f.push_back(0);f.push_back(4);f.push_back(5);f.push_back(1); m.addFace(f);f.clear();
	f.push_back(1);f.push_back(5);f.push_back(6);f.push_back(2); m.addFace(f);f.clear();
	f.push_back(2);f.push_back(6);f.push_back(7);f.push_back(3); m.addFace(f);f.clear();
	f.push_back(0);f.push_back(3);f.push_back(7);f.push_back(4); m.addFace(f);f.clear();
	VConverter::mesh_to_CGAL_Polyhedron_3(m, *cube);
	return cube;
}

void translate(Polyhedron_3 & poly, float x, float y, float z)
{
	CGAL::Aff_transformation_3<K> A(1,0,0,x,0,1,0,y,0,0,1,z,1);
	std::transform(poly.points_begin(),poly.points_end(),poly.points_begin(),A);
}

void scale(Polyhedron_3 & poly, float x, float y, float z)
{
	CGAL::Aff_transformation_3<K> A(x,0,0,0,0,y,0,0,0,0,z,0,1);
	std::transform(poly.points_begin(),poly.points_end(),poly.points_begin(),A);
}
Polyhedron_3 * testOperation()
{
	 Convex_hull_3 hullA(3);  // create instance of the class with dimension == 3
	  // generate 250 points randomly on a sphere of radius 100
	  // and insert them into the convex hull
	  Convex_hull_3 hullB(3);
	  CGAL::Random_points_in_sphere_3<Point_3, Creator> gen(2.5);
	  CGAL::Random_points_in_sphere_3<Point_3, Creator> genB(2.5);

	  for (int i = 0; i < 30 ; i++, ++gen)
		 hullA.insert(*gen);
	  for (int i = 0; i < 12 ; i++, ++genB)
		 hullB.insert(*genB);

	  assert(hullA.is_valid());
	  assert(hullB.is_valid());

	  // define polyhedron to hold convex hull and create it
	  Polyhedron_3 P; Polyhedron_3 PB;Polyhedron_3 * poly_intersection=new Polyhedron_3();
	  CGAL::convex_hull_d_to_polyhedron_3(hullA,P);
	  CGAL::convex_hull_d_to_polyhedron_3(hullB,PB);



	  CGAL::Nef_polyhedron_3<K> NP(P);
	  CGAL::Nef_polyhedron_3<K> NPB(PB);
	  CGAL::Nef_polyhedron_3<K> intersect=NP.intersection(NPB);
	  intersect.convert_to_polyhedron(*poly_intersection);

	  std::cout << "The convex hull A has " << P.size_of_vertices() 
				<< " vertices" << std::endl;
	  std::cout << "The convex hull B has " << PB.size_of_vertices() 
				<< " vertices" << std::endl;
	  std::cout << "intersection of A,B has" << poly_intersection->size_of_vertices()
				<< " vertices" << std::endl;

	  return poly_intersection;
}


void genRegularVoroCell(std::vector<Polyhedron_3 *> &cells, int num, glm::vec3 max, glm::vec3 min)
{
	const double x_min(min[0]),x_max(max[0]);
	const double y_min(min[1]),y_max(max[1]);
	const double z_min(min[2]),z_max(max[2]);
	const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

	// Set up the number of blocks that the container is divided into
	const int n_x=6,n_y=6,n_z=6;

	// Set the number of particles that are going to be randomly introduced
	int particles=num;
	int counter = 0;
	double dx = max[0] - min[0];
	double dy = max[1] - min[1];
	double dz = max[2] - min[2];
	double x,y,z, cx, cy, cz, px, py, pz;

	particle_order po;
	voronoicell c;
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	// Randomly add particles into the container
	int it = 0;

	for(int i=0;i<particles;i++) 
	{
		x=x_min+rnd()*dx;
		y=y_min+rnd()*dy;
		z=z_min+rnd()*dz;
		con.put(po,i,x,y,z);
	}


	c_loop_order clo(con,po);
	vector<int> face_vertices;
	vector<int> face_orders;

	if(clo.start()) 
	{
		do if(con.compute_cell(c,clo))
		{
			Mesh m;
			Polyhedron_3 * cell=new Polyhedron_3();
			clo.pos(px,py,pz);
			vector<double> vertices_global_position;
			c.vertices(px,py,pz,vertices_global_position);
			for(int i = 0; i < c.p; i++)
			{
				cx = vertices_global_position[3*i];
				cy = vertices_global_position[3*i+1];
				cz = vertices_global_position[3*i+2];
				m.addVertex(cx,cy,cz);
			}
			
			c.face_vertices(face_vertices);
			c.face_orders(face_orders);

			int num_vertices_in_face = face_vertices[0];
			int start = 1;
			int end = 0;
			vector<int> f,tmp;
			for(int i = 0; i<c.number_of_faces(); i++)
			{
				num_vertices_in_face = face_vertices[end];
				end = start + num_vertices_in_face;
				for(int j = start; j < end; j++)
				{	
					f.push_back(face_vertices[j]);
				}
				m.init();
				reverse(f.begin(),f.end());
				m.addFace(f);
				f.clear();
				start = end + 1;
			}
			VConverter::mesh_to_CGAL_Polyhedron_3(m, *cell);
			cells.push_back(cell);
			counter++;
		}while(clo.inc());
	}
}


void genVoroPatternCell(std::vector<Polyhedron_3 *> &cells, int num, glm::vec3 max, glm::vec3 min, glm::vec3 pos, float percentage, int type)
{

	const double x_min(min[0]),x_max(max[0]);
	const double y_min(min[1]),y_max(max[1]);
	const double z_min(min[2]),z_max(max[2]);
	const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

	// Set up the number of blocks that the container is divided into
	const int n_x=6,n_y=6,n_z=6;

	// Set the number of particles that are going to be randomly introduced
	const int particles=num;
	int i,it;
	double x,y,z, cx, cy, cz, px, py, pz;
	particle_order po;
	voronoicell c;
	double dx = x_max-x_min;
	double dy = y_max-y_min;
	double dz = z_max-z_min;
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	
	// Randomly add particles into the container
	it = 0;
	for(i=0;i<particles;i++) {

		do 
		{
			x=x_min+rnd()*(dx);
			y=y_min+rnd()*(dy);
			z=z_min+rnd()*(dz);
			it++;
		}
		while( (abs(x - pos[0]) > dx * percentage || abs(y - pos[1]) > dy * percentage || abs(z - pos[2]) > dz * percentage) && it < 1000);

		con.put(po,i,x,y,z);
	}
	c_loop_order clo(con,po);
	vector<int> face_vertices;
	vector<int> face_orders;

	if(clo.start()) 
	{
		do if(con.compute_cell(c,clo))
		{
			Mesh m;
			Polyhedron_3 * cell=new Polyhedron_3();
			clo.pos(px,py,pz);
			vector<double> vertices_global_position;
			c.vertices(px,py,pz,vertices_global_position);
			for(int i = 0; i < c.p; i++)
			{
				cx = vertices_global_position[3*i];
				cy = vertices_global_position[3*i+1];
				cz = vertices_global_position[3*i+2];
				m.addVertex(cx,cy,cz);
			}
			
			c.face_vertices(face_vertices);
			c.face_orders(face_orders);

			int num_vertices_in_face = face_vertices[0];
			int start = 1;
			int end = 0;
			vector<int> f,tmp;
			for(int i = 0; i<c.number_of_faces(); i++)
			{
				num_vertices_in_face = face_vertices[end];
				end = start + num_vertices_in_face;
				for(int j = start; j < end; j++)
				{	
					f.push_back(face_vertices[j]);
				}
				m.init();
				reverse(f.begin(),f.end());
				m.addFace(f);
				f.clear();
				start = end + 1;
			}
			VConverter::mesh_to_CGAL_Polyhedron_3(m, *cell);
			cells.push_back(cell);
		}while(clo.inc());
	}

	// Output the particle positions in gnuplot format
	//con.draw_particles("random_points_p.gnu");

	// Output the Voronoi cells in gnuplot format
	//con.draw_cells_gnuplot("random_points_v.gnu");

}

void genVoroCell(std::vector<Polyhedron_3 *> &cells, std::vector<Polyhedron_3 *> &otherCells, Polyhedron_3 * &impactCell,int num, glm::vec3 max, glm::vec3 min, glm::vec3 pos, float r, int type)
{
	const double x_min(min[0]),x_max(max[0]);
	const double y_min(min[1]),y_max(max[1]);
	const double z_min(min[2]),z_max(max[2]);
	const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);
	double minr = 1000000;

	// Set up the number of blocks that the container is divided into
	const int n_x=6,n_y=6,n_z=6;

	// Set the number of particles that are going to be randomly introduced
	int particles=num;
	int counter = 0;
	int impactID;
	
	double dx = max[0] - min[0];
	double dy = max[1] - min[1];
	double dz = max[2] - min[2];
	if(dx < minr)
		minr = dx;
	if(dy < minr)
		minr = dy;
	if(dz < minr)
		minr = dz;
	double radius = r * minr;
	double x,y,z, cx, cy, cz, px, py, pz;
	double impactX, impactY, impactZ;
	int div = 8;
	double k = 0.007;

	particle_order po;
	voronoicell c;
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	if( r >= 1)
	{
		particles = 1;
	}
	// Randomly add particles into the container
	int it = 0;
	if(type == 1)
	{
		for(int i=0;i<particles;i++) 
		{
			if( i == 0)
			{
				x = pos[0];
				y = pos[1];
				z = pos[2];
			}
			else
			{
				do 
				{
					x=pos[0] - radius + rnd()*(radius*2.0f);
					y=pos[1] - radius + rnd()*(radius*2.0f);
					z=pos[2] - radius + rnd()*(radius*2.0f);
					it++;
				}
				while( (x >= x_max  || y  >= y_max||  z>=z_max || x <=x_min || y<=y_min || z <= z_min) && it < 1000);
			}
		con.put(po,i,x,y,z);
		}
	}
	else if(type == 2)                        // spider web pattern
	{
		double divx = dx/(div*2+1.0f);
		double divy = dy/(div*2+1.0f);
		double divz = dz/(div*2+1.0f);
		int cp = 0;
		for(int i = 0; i < div ; i++)
		{
			x = pos[0] + (i+1)*divx + (rand()%3-1.0f)*dx*k;
			y = pos[1] + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;

			x = pos[0] - (i+1)*divx + (rand()%3-1.0f)*dx*k;
			y = pos[1] + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;
		}

		for(int i = 0; i < div ; i++)
		{
			x = pos[0] + (rand()%3-1.0f)*dx*k ;
			y = pos[1] + (i+1)*divy + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;

			x = pos[0] + (rand()%3-1.0f)*dx*k ;
			y = pos[1] - (i+1)*divy + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;
		}

		for(int i = 0; i < div ; i++)
		{
			x = pos[0] + (rand()%3-1.0f)*dx*k;
			y = pos[1] + (rand()%3-1.0f)*dy*k;
			z = pos[2] +(i+1)*divz + (rand()%3-1.0f)*dy*k;
			con.put(po,cp,x,y,z);
			cp++;

			x = pos[0] + (rand()%3-1.0f)*dx*k;
			y = pos[1] + (rand()%3-1.0f)*dy*k;
			z = pos[2] - (i+1)*divz + (rand()%3-1.0f)*dy*k;
			con.put(po,cp,x,y,z);
			cp++;
		}

		for(int i = 0; i < div ; i++)
		{
			x = pos[0] + (i+1)*divx + (rand()%3-1.0f)*dx*k;
			y = pos[1] + (i+1)*divy + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;

			x = pos[0] - (i+1)*divx + (rand()%3-1.0f)*dx*k;
			y = pos[1] - (i+1)*divy + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;

		}

		for(int i = 0; i < div ; i++)
		{
			x = pos[0] - (i+1)*divx + (rand()%3-1.0f)*dx*k;
			y = pos[1] + (i+1)*divy + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;

			x = pos[0] + (i+1)*divx + (rand()%3-1.0f)*dx*k;
			y = pos[1] - (i+1)*divy + (rand()%3-1.0f)*dy*k;
			z = pos[2] + (rand()%3-1.0f)*dz*k;
			con.put(po,cp,x,y,z);
			cp++;
		}
	}
	c_loop_order clo(con,po);
	vector<int> face_vertices;
	vector<int> face_orders;

	if(clo.start()) 
	{
		do if(con.compute_cell(c,clo))
		{
			Mesh m;
			Polyhedron_3 * cell=new Polyhedron_3();
			clo.pos(px,py,pz);
			vector<double> vertices_global_position;
			c.vertices(px,py,pz,vertices_global_position);
			for(int i = 0; i < c.p; i++)
			{
				cx = vertices_global_position[3*i];
				cy = vertices_global_position[3*i+1];
				cz = vertices_global_position[3*i+2];
				m.addVertex(cx,cy,cz);
			}
			
			c.face_vertices(face_vertices);
			c.face_orders(face_orders);

			int num_vertices_in_face = face_vertices[0];
			int start = 1;
			int end = 0;
			vector<int> f,tmp;
			for(int i = 0; i<c.number_of_faces(); i++)
			{
				num_vertices_in_face = face_vertices[end];
				end = start + num_vertices_in_face;
				for(int j = start; j < end; j++)
				{	
					f.push_back(face_vertices[j]);
				}
				m.init();
				reverse(f.begin(),f.end());
				m.addFace(f);
				f.clear();
				start = end + 1;
			}
			VConverter::mesh_to_CGAL_Polyhedron_3(m, *cell);
			cells.push_back(cell);
			if (counter == 0)
				impactCell = cell;
			else
			{
				otherCells.push_back(cell);
				cout << "Push into other cell" << endl;
			}
			counter++;
		}while(clo.inc());
	}
}

void dimensionPolyhedron_3( Polyhedron_3 & poly, float & minX,  float & maxX,  float & minY,  float & maxY,  float & minZ,  float & maxZ )
{
	float x(0.0f),X(0.0f),y(0.0f),Y(0.0f),z(0.0f),Z(0.0f);
	CGAL::Polyhedron_3<K>::Point_iterator iter;
	for(iter=poly.points_begin();iter!=poly.points_end();++iter)
	{
		if(iter->x().to_double()<x) x=iter->x().to_double();
		if(iter->x().to_double()>X) X=iter->x().to_double();

		if(iter->y().to_double()<y) y=iter->y().to_double();
		if(iter->y().to_double()>Y) Y=iter->y().to_double();

		if(iter->z().to_double()<z) z=iter->z().to_double();
		if(iter->z().to_double()>Z) Z=iter->z().to_double();
	}
	minX=x;maxX=X;minY=y;maxY=Y;minZ=z;maxZ=Z;
	
}

bool withinDimension(glm::vec3 pos, glm::vec3 max, glm::vec3 min)
{
	if(pos[0]<=min[0] || pos[0]>=max[0]) return false;
	if(pos[1]<=min[1] || pos[1]>=max[1]) return false;
	if(pos[2]<=min[2] || pos[2]>=max[2]) return false;

	return true;
}

void clipDimension(glm::vec3 & pos, glm::vec3 & max, glm::vec3 & min)
{
	float eX(0.0),eY(0),eZ(0);
	eX=0.005f*abs(min[0]-max[0]);eY=0.005f*abs(min[1]-max[1]);eZ=0.005f*abs(min[2]-max[2]);
	if(pos[0]<=min[0]) pos[0]=min[0]+eX;
	if(pos[1]<=min[1]) pos[1]=min[1]+eY;
	if(pos[2]<=min[2]) pos[2]=min[2]+eZ;

	if(pos[0]>=max[0]) pos[0]=max[0]-eX;
	if(pos[1]>=max[1]) pos[1]=max[1]-eY;
	if(pos[2]>=max[2]) pos[2]=max[2]-eZ;
}