/**
**  by Zhen Gou, Zhenghan Mei 2014
**/
#include "SceneGraph.h"


Node::Node():parent(NULL),furniture(NULL),children(NULL)
{	
	selected=false;
	mult=false;
	red=1.0;
	green=0.0;
	blue=0.0;
	geometry=-1;
	material=1;
	mesh=-1;
	meshHeight=1.0;

	transX=0.0;
	transZ=0.0;
	transY=0.0;
	rotationY=0.0;
	scaleX=1.0;
	scaleY=1.0;
	scaleZ=1.0;
	multX=1;
	multZ=1;
}

void Node::addChild(Node n)
{
	
	children.push_back(n);
	children[children.size()-1].parent=this;

}
void Node::select()
{
	selected=true;
	if(furniture!=NULL)
	{
		furniture->select();
	}
	for (int i=0;i<children.size();++i)
	{
		children[i].select();
	}
}

void Node::deselect()
{
	selected=false;
	if(furniture!=NULL)
	{
		furniture->deselect();
	}
	for (int i=0;i<children.size();++i)
	{
		children[i].deselect();
	}
}

float Node::getHeight()
{
	if (children.size()==0||mult)
	{
		float H=0.0;
		
		if(furniture!=NULL)
		{
			float furH=furniture->getHeight();
			H=std::max(H,furH);
		}
		if(geometry>=0)
		{
			float primitiveH=0.5;
			H=std::max(primitiveH,H);
		}
		if(mesh>=0)
		{
			H=std::max(meshHeight,H);
		}

		return transY+scaleY*H;
	}
	else
	{
		float H=0.0;
		for (int i=0;i<children.size();++i)
		{
			float h=children[i].getHeight();
			H=std::max(h,H);
		}
		if (geometry>=0)
		{
			float primitiveH=0.5;
			H=std::max(H,primitiveH);
		}

		if(mesh>=0)
		{
			H=std::max(meshHeight,H);
		}

		if (furniture!=NULL)
		{
			float furH=furniture->getHeight();
			H=std::max(furH,H);
		}

		return transY+scaleY*H;
		
	}

}

void Node::floodChildMaterial()
{
	for (int i=0;i<children.size();++i)
	{
		children.at(i).material=material;
	}
}
