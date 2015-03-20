/**
**  by Zhen Gou, Zhenghan Mei 2014
**/
#ifndef MESH_H
#define MESH_H
#include <vector>
#include "SceneGraph.h"
using namespace std;
using namespace glm;

class Mesh
{
public:
	int id;
	float r,g,b;
	//refers to indices of vertex, a face comprises of 3 int, SHOULD TO GO IBO
	vector<int> faceList;
	//indicate num of vertices of each face
	vector<int> faceCount;
	//a vertex consists of 4 floats  SHOULD GO TO VBO
	vector<float> vertexList;
	
	//a normal consists of 4 floats, and it's normal for each vertex  SHOULD GO TO NBO
	vector<float> normalList;
	
	//a color consists of 3 floats   SHOULD GO TO CBO
	vector<float> colorList;
	//whether init() has been called
	bool initialized;

	void scale(float sx,float sy,float sz)
	{

		glm::mat4 M=glm::mat4(1.0);
		M=glm::scale(M,vec3(sx,sy,sz));
		for (int i=0;i<numVertices();++i)
		{
			vec4 v(getVertex(i),1.0);
			v=M*v;
			setVertex(i,v[0],v[1],v[2]);
		}

		
		
	}
	int numVertices()
	{
		return (vertexList.size()/4);
	}
	int numFaces()
	{
		return (faceList.size()/3);
	}
	void addVertex(float x, float y, float z)
	{
		vertexList.push_back(x);
		vertexList.push_back(y);
		vertexList.push_back(z);
		vertexList.push_back(1.0);
	}

	//only call after loaded all vertices
	void init()
	{
		if(initialized) return;
		for(int i=0;i<vertexList.size()/4;++i)
		{
			normalList.push_back(0.0);normalList.push_back(0.0);normalList.push_back(0.0);normalList.push_back(1.0);
		}
		for(int i=0;i<vertexList.size()/4;++i)
		{
			colorList.push_back(1.0);colorList.push_back(0.0);colorList.push_back(0.0);
		}
		initialized=true;
	}

	//should only be called after vertices have been loaded
	void addFace(vector<int> & list)
	{
		if(!initialized) cout<<"ERROR: Mesh, not initialized!"<<endl;
		int num=list.size();
		if(num<3) { cout<<"VConverter::VMesh_toMesh::loading face and face has less than 3 sides"<<endl; } //pop msg when face has < 3 sides
		for(int j=1;j<num-1;++j)
		{
			int i1,i2,i3; //index of the 3 vertices
			vec3 v1,v2,v3,u,v,normal; //vec3 of the 3 vertices
			i1=list.at(0);i2=list.at(j);i3=list.at(j+1);	
			v1=getVertex(i1);v2=getVertex(i2);v3=getVertex(i3);

			u=v2-v1;v=v3-v2;
			normal=glm::normalize(glm::cross(u,v));
			setNormal(i1,normal); setNormal(i2,normal); setNormal(i3,normal);
			addFace(i1,i2,i3);
		}
	}
	//wrong method, do not use
	void addFaceWithJagginess(vector<int> & list, float epsilon, int num)
	{
		//calculate center of mass
		vec3 com(0,0,0);
		for(int i=0;i<list.size();++i)
			com+=getVertex(list.at(i));
		com*=1/list.size();
		////////////////////////////////////////////////////////

		//foreach 2 vertex in the original list
		
		
		vector<int> result;result.clear();
		result.push_back(list.at(0));
			
		for(int i=0;i<list.size()-1;++i)
		{
			vec3 A=	getVertex(list.at(i));
			vec3 B=	getVertex(list.at(i+1));
			vec3 C=0.5f*(B-A)+A;
		//	C=(1.0f-epsilon)*(C-com)+com;

			int index=numVertices();
			
		addVertex(C[0],C[1],C[2]);
		result.push_back(index);




			result.push_back(list.at(i+1));
		}
		
		/*
		for(int i=0;i<list.size()-1;++i)
		{
			result.push_back(list.at(i));
		}*/
		addFace(list);
	}

	vec3 getVertex(int i)
	{
		vec3 v=vec3(vertexList.at(4*i),vertexList.at(4*i+1),vertexList.at(4*i+2));
		return v;
	}

	void setVertex(int i,float x, float y, float z)
	{
		vertexList.at(4*i)=x;vertexList.at(4*i+1)=y;vertexList.at(4*i+2)=z;vertexList.at(4*i+3)=1.0;
	}
	


	Mesh():id(-1),r(0.5),g(0.5),b(0.5),initialized(false){}

private:
	void addFace(int a, int b, int c)
	{
		faceList.push_back(a);
		faceList.push_back(b);
		faceList.push_back(c);
	}

	void addNormal(float x, float y, float z)
	{
		normalList.push_back(x);
		normalList.push_back(y);
		normalList.push_back(z);
		normalList.push_back(1.0);
	}

	void setNormal(int i,vec3 norm)
	{
		normalList.at(4*i)+=norm[0];
		normalList.at(4*i+1)+=norm[1];
		normalList.at(4*i+2)+=norm[2];
		normalList.at(4*i+3)=1.0;
	}
};
#endif