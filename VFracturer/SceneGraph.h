/**
**  by Zhen Gou, Zhenghan Mei 2014
**/
#ifndef SG_H
#define SG_H
#include <vector>
#include "../glm/glm.hpp"
#include "../glm/gtc/matrix_transform.hpp"
using namespace std;
class Node
{	
public:
	
	glm::mat4 transform;
	Node* parent;
	Node* furniture;
	vector<Node> children;

	float red;
	float green;
	float blue;
	float meshHeight;
	int material;

	int geometry;
	int mesh;

	float transX;
	float transZ;
	float transY;
	float rotationY;
	float scaleX;
	float scaleY;
	float scaleZ;

	int multX;
	int multZ;

	bool selected;
	bool mult;

	Node();
	float getHeight();
	void select();
	void deselect();
	void addChild(Node n);
	void floodChildMaterial();

};
#endif
