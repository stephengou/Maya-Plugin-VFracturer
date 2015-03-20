/** Program for real-time view of the fracture simulation plug-in, mainly for debugging purposes.
**  This code is based on my raytracer code from CIS 560 final project.
**  
**  Zhen Gou
**/
#include "glew.h"
#include "../freeglut/include/GL/glut.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>

#include "VOperations.h"
#include "VConverter.h"
#include "VClipper.h"
using namespace std;
using namespace glm;


//mesh buffer
unsigned int m_vbo=0;
unsigned int m_cbo=1;
unsigned int m_nbo=2;
unsigned int m_ibo=3;

//attributes
unsigned int positionLocation = 0;
unsigned int colorLocation = 1;
unsigned int normalLocation = 2;

//uniforms
unsigned int u_modelMatrixLocation;
unsigned int u_projMatrixLocation;
unsigned int u_lightPos;
unsigned int u_lightColor;
vec3 lightPos=vec3(0.0,0.0,5.0);
vec3 lightColor=vec3(1.0,1.0,1.0);
	
//needed to compile and link and use the shaders
unsigned int vertexShader;
unsigned int fragmentShader;
unsigned int shaderProgram;

//Animation/transformation stuff
clock_t old;
float rotation = 0.0f;

//helper function to read shader source and put it in a char array
//thanks to Swiftless
char* textFileRead(const char*);

//some other helper functions from CIS 565
void printLinkInfoLog(int);
void printShaderInfoLog(int);

//standard glut-based program functions
void init(void);
void resize(int, int);
void display(void);
void keypress(unsigned char, int, int);
void cleanup(void);

void createMesh(mat4 modelView, int id,float r,float g, float b);

void initRoom(string);

void traverseSceneGraph(Node&, mat4);
void generatePreorderNodeList(Node* root);
vector<string> splitStringOnSpace(string);
vector<string> readLineSplitBySpace();
void loadItem();
void loadPoly_3(Polyhedron_3 * mesh);

void generateMesh(int id,string file, string texture);
void generateMeshFromObj(string file);
void addChildIntellignet(Node &root, Node n);
void deleteNode(Node* n);

ifstream roomFile;
int xSize;
int zSize;
int numItems;


Node* room=new Node();
vector<Node*> preorderNodeList;
////refers to index in the preorderNodeList
int currentNode;

vector<Mesh> meshList;

vec4 crossProduct(vec4 u,vec4 v);
float norm(vec4 v);
bool isConvex(vector<vec4> list);

int main(int argc, char** argv) {	


	string file=argv[1];
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(640, 480);
	glutCreateWindow("VFracturer OpenGL view");
	
	//Call GLEW only _after_ you get the window
	//I had to tell the author of your textbook that ;-)  -Cory
	glewInit();
	init();
	initRoom(file);
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keypress);
	glutIdleFunc(display);
	glutMainLoop();

	return 0;
}

void initRoom(string file)
{
	roomFile.open(file);
	/////read in room parameters
	vector<string> para=readLineSplitBySpace();
	xSize=atoi(para[0].c_str());
	zSize=atoi(para[1].c_str());
	numItems=atoi(para[2].c_str());
	loadItem();

	std::vector<Polyhedron_3 *> cells;
	genVoroCell(cells, 7);

	meshList.at(0).scale(1, 1, 1);

	Polyhedron_3 * P=new Polyhedron_3();
    VConverter::mesh_to_CGAL_Polyhedron_3(meshList.at(0),*P);
	VClipper clipper=VClipper(*P);
	//load each cell to the clipper
	for(int i = 0; i < cells.size(); ++i)
	{
		clipper.loadPatternCell(*cells[i]);
	}
	vector<Polyhedron_3 *> result;
	clipper.output(result);
	//load resulting poly_3 for display in OpenGL
	for(int i=0;i<result.size();++i)
	{
		loadPoly_3(result.at(i));
	}

	generatePreorderNodeList(room);
	cout<<"number of nodes: "<<preorderNodeList.size()<<endl;
	currentNode=0;
	preorderNodeList[currentNode]->select();
	roomFile.close();
}


void loadPoly_3(Polyhedron_3 * mesh)
{
	int x,y,z,mX,mZ,mat,id;x=y=z=0.0;
	float rot, r=0.0,g=0.0,b=1.0,sx,sy,sz;sx=sy=sz=1.0; rot=0.0;
	cout<<"loading Polyhedron_3..."<<endl;
	float rr,gg,bb;
	rr=gg=bb=0.5;
	id=meshList.size();

	Mesh m;
	VConverter::CGAL_Polyhedron_3_to_Mesh(*mesh,m);
	meshList.push_back(m);

	Node item=Node();
	item.transX=x; item.transZ=z;item.transY=y;
	item.rotationY=rot;
	item.scaleX=sx;item.scaleY=sy;item.scaleZ=sz;
	item.mesh=id;
	cout<<"mesh generated"<<"ID is  "<<id<<endl;
	addChildIntellignet(*room, item);
	
}
void loadItem()
{
	

	for (int i=0;i<numItems;++i)
	{
		int x,y,z,mX,mZ,mat;
		float rot, r=0.0,g=0.0,b=1.0,sx,sy,sz;
		int meshID=meshList.size();
		string fur;
		string line;
		vector<string> para;
		///////fur name
		getline(roomFile,line);
		getline(roomFile,line);
		para=splitStringOnSpace(line);
		fur=para[0].c_str(); 

		if(fur=="obj")
		{
			cout<<"loading obj"<<endl;
			float rr,gg,bb;
			int mat;
			rr=gg=bb=0.5;
	
			string file;
			string texture;

			getline(roomFile,line);
			file=line;
			getline(roomFile,line);
			mat=atoi(line.c_str());

			getline(roomFile,line);
			para=splitStringOnSpace(line);
			x=atoi(para[0].c_str())*xSize;y=atoi(para[1].c_str());z=atoi(para[2].c_str())*zSize;

			/////fur rot
			getline(roomFile,line);
			para=splitStringOnSpace(line);
			rot=atof(para[0].c_str());
			//////fur scale
			getline(roomFile,line);
			para=splitStringOnSpace(line);
			sx=atof(para[0].c_str());sy=atof(para[1].c_str());sz=atof(para[2].c_str());

			generateMeshFromObj(file);

			Node item=Node();
			item.transX=x; item.transZ=z;item.transY=y;
			item.rotationY=rot;
			item.scaleX=sx;item.scaleY=sy;item.scaleZ=sz;
			item.mesh=meshID;
			cout<<"mesh generated"<<"ID is  "<<meshID<<endl;
			addChildIntellignet(*room, item);
		}

	}

}
void generateMeshFromObj(string file)
{
	int faceCount=0;
	//TO-DO HEIGHT
	ifstream File;
	string line;
	string tag;
	vector<string> para;
	File.open(file);
	cout<<"obj file opened start parsing"<<endl;
	
	Mesh m=Mesh();
	float x,y,z;int v1,v2,v3;
	while(!File.eof())
	{
		getline(File,line);
		if(line.length()<1)
		{
			continue;
		}
		para=splitStringOnSpace(line);
		tag=para[0];

		if(tag=="v")
		{
			x=atof(para[1].c_str());y=atof(para[2].c_str());z=atof(para[3].c_str());
			m.addVertex(x,y,z);

		}
		else if(tag=="f")
		{
			m.init();
			v1=atoi(para[1].c_str())-1;v2=atoi(para[2].c_str())-1;v3=atoi(para[3].c_str())-1;
			vector<int> face;face.push_back(v1);face.push_back(v2);face.push_back(v3);
			m.addFace(face);
			faceCount++;
		}
	}
	
	int vNum=m.vertexList.size()/4;
	//color and indices *facelist
	for (int i=0;i<vNum;++i)
	{
		m.colorList.push_back(m.r);
		m.colorList.push_back(m.g);
		m.colorList.push_back(m.b);
	}
	cout<<"obj file loading completed"<<endl;
	meshList.push_back(m);
	File.close();
}
void addChildIntellignet(Node &root, Node n)
{
		root.addChild(n);
	    cout<<"child added succesfully"<<endl;

}
void deleteNode(Node* n)
{
	Node* p=n->parent;
	int index=0;;
	for (int i=0;i<p->children.size();++i)
	{
		Node* c=&p->children[i];
		if(c==n)
		{
			index=i;
		}
		
	}
	vector<Node> *childrenList=&p->children;
	childrenList->erase(childrenList->begin()+index);
	generatePreorderNodeList(room);
}
 vector<string> splitStringOnSpace(string line)
    {
        string buffer;
        stringstream ss(line); 
        vector<std::string> tokens; 
        while (ss >> buffer)
        {    tokens.push_back(buffer);
           }
        return tokens;
    }
 vector<string> readLineSplitBySpace()
{
    string line;
        getline(roomFile,line);
    string buffer;
    stringstream ss(line); 
    vector<std::string> tokens; 
    while (ss >> buffer)
    {    tokens.push_back(buffer);
        }
    return tokens;
}
void init() {
	
	//mesh buffer, all mesh share the same buffer
	glGenBuffers(1, &m_vbo);
	glGenBuffers(1, &m_cbo);
	glGenBuffers(1, &m_nbo);
	glGenBuffers(1, &m_ibo);
	
	//Everybody does this
	glClearColor(0, 0, 0, 1);
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0);
	glDepthFunc(GL_LEQUAL);

	
	//here is stuff for setting up our shaders
	const char* fragFile = "diffuseFrag.frag";
	const char* vertFile = "diffuseVert.vert";
	vertexShader = glCreateShader(GL_VERTEX_SHADER);
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	shaderProgram = glCreateProgram();
	
	//load up the source, compile and link the shader program
	const char* vertSource = textFileRead(vertFile);
	const char* fragSource = textFileRead(fragFile);
	glShaderSource(vertexShader, 1, &vertSource, 0);
	glShaderSource(fragmentShader, 1, &fragSource, 0);
	glCompileShader(vertexShader);
	glCompileShader(fragmentShader);

	//For your convenience, i decided to throw in some compiler/linker output helper functions
	//from CIS 565
	GLint compiled;
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(vertexShader);
	} 
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(fragmentShader);
	} 
	
	//set the attribute locations for our shaders
	glBindAttribLocation(shaderProgram, positionLocation, "vs_position");
	glBindAttribLocation(shaderProgram, normalLocation, "vs_normal");
	glBindAttribLocation(shaderProgram, colorLocation, "vs_color");

	//finish shader setup
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);

	//check for linking success
	GLint linked;
	glGetProgramiv(shaderProgram,GL_LINK_STATUS, &linked);
	if (!linked) 
	{
		printLinkInfoLog(shaderProgram);
	}

	//Get the uniform locations for our shaders, unfortunately they can not be set by us, we have
	//to ask OpenGL for them
	u_modelMatrixLocation = glGetUniformLocation(shaderProgram, "u_modelMatrix");
	u_projMatrixLocation = glGetUniformLocation(shaderProgram, "u_projMatrix");
	u_lightPos=glGetUniformLocation(shaderProgram, "u_lightPos");
	u_lightColor=glGetUniformLocation(shaderProgram, "u_lightColor");

	//Always remember that it doesn't do much good if you don't have OpenGL actually use the shaders
	glUseProgram(shaderProgram);

	resize(640,480);
	
	old = clock();
	std::cout<<"initialization succeeded.."<<std::endl;

}
void cleanup() {
	glDeleteBuffers(1, &m_vbo);
	glDeleteBuffers(1, &m_cbo);
	glDeleteBuffers(1, &m_nbo);
	glDeleteBuffers(1, &m_ibo);
	//delete all VBO for mesh

	//Tear down the shader program in reverse of building it
	glDetachShader(shaderProgram, vertexShader);
	glDetachShader(shaderProgram, fragmentShader);
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
	glDeleteProgram(shaderProgram);
}
void keypress(unsigned char key, int x, int y) {
	Node* n=preorderNodeList[currentNode];
	switch(key) {

	case 'q':
		cleanup();
		exit(0);
		break;
	case 'n':
		if(currentNode>=preorderNodeList.size())
			currentNode=0;
		preorderNodeList[currentNode]->deselect();

		currentNode++;
		if(currentNode>=preorderNodeList.size())
			currentNode=0;

		preorderNodeList[currentNode]->select();

		for (int i=0;i<preorderNodeList.size();++i)
		{
			cout<<preorderNodeList[i]->selected<<endl;
		}
		cout<<endl;
		break;
	case 'a':
	
		n->transX+=-0.5;
		break;

	case 'd':
	
		n->transX+=0.5;
		break;
	case 'w':
	
		n->transZ+=0.5;
		break;
	case 's':
	
		n->transZ+=-0.5;
		break;

	case '2':
		n->transY+=0.05;
		break;

	case '1':
		n->transY+=-0.05;
		break;

	case 'x':

		n->scaleX+=-0.5;
		break;

	case 'X':

		n->scaleX+=0.5;
		break;

	case 'y':

		n->scaleY+=-0.5;
		break;

	case 'Y':

		n->scaleY+=0.5;
		break;
	case 'z':

		n->scaleZ+=-0.5;
		break;

	case 'Z':

		n->scaleZ+=0.5;
		break;

	case 'r':

		n->rotationY+=-10.0;
		break;

	case 'R':

		n->rotationY+=10.0;
		break;

	case 'e':
		deleteNode(n);
		break;
	}

	glutPostRedisplay();
}

void display() {
	//Always and only do this at the start of a frame, it wipes the slate clean
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	clock_t newTime = clock();

	//part of the animation
	rotation += 150 * (static_cast<float>(newTime - old) / static_cast<float>(CLOCKS_PER_SEC));

	//create an identity matrix for the modelview matrix
	float rot=30;
	glm::mat4 modelView = glm::mat4(1.0);

	//adjusting camera view
	//modelView=glm::scale(modelView,vec3(0.3,0.3,0.3));
	modelView=glm::translate(modelView,vec3(0,0,-10.0));
	modelView=glm::rotate(modelView,rot,vec3(1,0,0));
	//modelView=glm::rotate(modelView,rotation/3,vec3(1,0,0));
	//modelView=glm::scale(modelView,vec3(0.5,0.5,0.5));
	//traverse the scenegraph to render the room
	traverseSceneGraph(*room,modelView);
	
	glutSwapBuffers();
	old = newTime;
}
void generatePreorderNodeList(Node* root)
{
	preorderNodeList.push_back(root);
	for(int i=0;i<root->children.size();++i)
	{
		Node* n=&root->children[i];
		generatePreorderNodeList(n);
	}

}
void traverseSceneGraph(Node& n, mat4 T )
{
	mat4 Ti=mat4(1.0);


	Ti=glm::translate(Ti,vec3(n.transX,n.transY,n.transZ));
	Ti=glm::rotate(Ti,n.rotationY,vec3(0,1,0));
	Ti=glm::scale(Ti,vec3(n.scaleX,n.scaleY,n.scaleZ));	


	T=T*Ti;

	if (n.mesh>=0)
	{
		if(n.selected)
			createMesh(T,n.mesh,0,1,0);
		else
			createMesh(T,n.mesh,n.red,n.green,n.blue);
	}

	for (int i=0;i<n.children.size();++i)
	{
		Node child=n.children[i];
		traverseSceneGraph(child,T);
		
	}
}

//id refers to mesh id
void createMesh(mat4 modelView, int id,float r,float g, float b)
{
	Mesh& m=meshList.at(id);
	int vSize=m.vertexList.size();
	int cSize=m.colorList.size();
	int iSize=m.faceList.size();

	//cout<<"v: "<<vSize<<"   c:   "<<cSize<<"   i:   "<<iSize<<endl;
	//activate our three kinds of information

	float * vertices=new float[vSize];
	for (int i=0;i<vSize;++i)
	{
		vertices[i]=m.vertexList.at(i);
		//cout<<vertices[i]<<endl;
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
	glBufferData(GL_ARRAY_BUFFER, vSize * sizeof(float), vertices, GL_STATIC_DRAW); 

	delete[] vertices;

	float * colors=new float[cSize];
	for (int i=0;i<cSize/3;++i)
	{
		colors[i*3]=r;colors[i*3+1]=g;colors[i*3+2]=b;
		//cout<<colors[i]<<endl;
	}

	
	glBindBuffer(GL_ARRAY_BUFFER, m_cbo);
	glBufferData(GL_ARRAY_BUFFER, cSize * sizeof(float), colors, GL_STATIC_DRAW); 

	delete[] colors;

	float * normals=new float[vSize];
	for (int i=0;i<vSize;++i)
	{
		normals[i]=m.normalList.at(i);
		//cout<<normals[i]<<endl;
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_nbo);
	glBufferData(GL_ARRAY_BUFFER, vSize * sizeof(float), normals, GL_STATIC_DRAW); 

	delete[] normals;

	unsigned short * indices=new unsigned short[iSize];
	//cout<<"loading ibo"<<endl;
	for (int i=0;i<iSize;++i)
	{
		indices[i]=m.faceList.at(i);
		//cout<<indices[i]<<endl;
	}
	//cout<<"Done loading ibo"<<endl;
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, iSize* sizeof(unsigned short), indices, GL_STATIC_DRAW);

	delete[] indices;

	glEnableVertexAttribArray(positionLocation);
	glEnableVertexAttribArray(colorLocation);
	glEnableVertexAttribArray(normalLocation);
	
	//we're using the vertex data first
	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
	//define how the vertex pointer should work, in our case we're accessing floats 4 at a time with no special pattern
	glVertexAttribPointer(positionLocation, 4, GL_FLOAT, 0, 0, static_cast<char*>(0));
	
	//now use color data, remember we're not using 4 at a time anymore
	glBindBuffer(GL_ARRAY_BUFFER, m_cbo);
	glVertexAttribPointer(colorLocation, 3, GL_FLOAT, 0, 0, static_cast<char*>(0));
	
	//one more time with the normals
	glBindBuffer(GL_ARRAY_BUFFER, m_nbo);
	glVertexAttribPointer(normalLocation, 4, GL_FLOAT, 0, 0, static_cast<char*>(0));
	
	//the last thing we need to do is setup our indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);

	//set the modelview uniform
	glUniformMatrix4fv(u_modelMatrixLocation,1, GL_FALSE, &modelView[0][0]);
	glUniform3fv(u_lightPos, 1, &lightPos[0]);
	//draw the elements
	
	glDrawElements(GL_TRIANGLES, iSize, GL_UNSIGNED_SHORT, 0);
	//glDrawArrays(GL_POINTS,0,iSize);
	
	//shut off the information since we're done drawing
	glDisableVertexAttribArray(positionLocation);
	glDisableVertexAttribArray(colorLocation);
	glDisableVertexAttribArray(normalLocation);
}
void resize(int width, int height) {
	//set the viewport, more boilerplate
	glViewport(0, 0, width, height);

	//
	glm::mat4 projection = glm::perspective(60.0f, static_cast<float>(width) / static_cast<float>(height), 0.1f, 30.0f);
	glm::mat4 camera = glm::lookAt(glm::vec3(0, 0, 10), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	projection = projection * camera;

	//set the projection matrix here, it only needs to be changed if the screen is resized otherwise it can stay the same
	glUniformMatrix4fv(u_projMatrixLocation, 1, GL_FALSE, &projection[0][0]);

	glutPostRedisplay();
}

//from swiftless.com
char* textFileRead(const char* fileName) {
    char* text;
    
    if (fileName != NULL) {
        FILE *file = fopen(fileName, "rt");
        
        if (file != NULL) {
            fseek(file, 0, SEEK_END);
            int count = ftell(file);
            rewind(file);
            
            if (count > 0) {
                text = (char*)malloc(sizeof(char) * (count + 1));
                count = fread(text, sizeof(char), count, file);
                text[count] = '\0';	//cap off the string with a terminal symbol, fixed by Cory
            }
            fclose(file);
        }
    }
    return text;
}

void printLinkInfoLog(int prog) 
{
	int infoLogLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;

	glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &infoLogLen);

	// should additionally check for OpenGL errors here

	if (infoLogLen > 0)
	{
		infoLog = new GLchar[infoLogLen];
		// error check for fail to allocate memory omitted
		glGetProgramInfoLog(prog,infoLogLen, &charsWritten, infoLog);
		std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
		delete [] infoLog;
	}
}

void printShaderInfoLog(int shader)
{
	int infoLogLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;

	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

	// should additionally check for OpenGL errors here

	if (infoLogLen > 0)
	{
		infoLog = new GLchar[infoLogLen];
		// error check for fail to allocate memory omitted
		glGetShaderInfoLog(shader,infoLogLen, &charsWritten, infoLog);
		std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
		delete [] infoLog;
	}

	// should additionally check for OpenGL errors here
}

vec4 crossProduct(vec4 u, vec4 v)
{
	vec4 i,j,k,ret;

	i[0]=1;i[1]=0;i[2]=0;i[3]=1;
	j[0]=0;j[1]=1;j[2]=0;j[3]=1;
	k[0]=0;k[1]=0;k[2]=1;k[3]=1;
	ret=(u[1]*v[2]-u[2]*v[1])*i+(u[2]*v[0]-u[0]*v[2])*j+(u[0]*v[1]-u[1]*v[0])*k;
	ret[3]=0;
	return ret;

}
float norm(vec4 v)
{
	float ret;
	ret=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	return ret;
}
bool isConvex(vector<vec4> list)
{
	float lastY;
	vec4 u1,v1,norm1;
	u1=list.at(0)-list.at(list.size()-2);
	v1=list.at(1)-list.at(0); 
	norm1=crossProduct(u1,v1);
	lastY=norm1[1];
	cout<<"y :"<<norm1[1]<<endl;

	for (int i=0;i<list.size()-2;++i)
	{
		vec4 u,v,norm;
		u=list.at(i+1)-list.at(i);
		v=list.at(i+2)-list.at(i+1); 
		norm=crossProduct(u,v);
		cout<<"y :"<<norm[1]<<endl;
		if(norm[1]*lastY<0)
		{
			return false;
		}

		//lastY=norm[1];
	}
	return true;
}
