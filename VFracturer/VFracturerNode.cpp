/**
**  by Zhen Gou, Zhenghan Mei 2014
**/
#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>
#include "VFracturerNode.h"
#include "VClipper.h"
#include "VOperations.h"

MTypeId VFracturerNode::id(0x00321);
MObject VFracturerNode::iX;
MObject VFracturerNode::iY;
MObject VFracturerNode::iZ;
MObject VFracturerNode::iR;
MObject VFracturerNode::iPer;
MObject VFracturerNode::dNum;
MObject VFracturerNode::pNum;
MObject VFracturerNode::pType;
MObject VFracturerNode::outputMesh;
MObject VFracturerNode::inputMesh;
MObject VFracturerNode::fractureMode;
MObject VFracturerNode::objPatternNum;
MObject VFracturerNode::objPatternFile;
const float RADIUS=0.25;

MString toString(int i);
MString toString(double f);
vector<string> splitStringOnSpace(string line);
void generateMeshFromObj(string file);

vector<Mesh> objPatternCellList;

VFracturerNode::VFracturerNode()
{
}
VFracturerNode::~VFracturerNode()
{
}

void *VFracturerNode::creator()
{
	return new VFracturerNode();
}

MStatus VFracturerNode::initialize()
{
	MFnUnitAttribute unitA;
	MFnNumericAttribute numA;
	MFnTypedAttribute typeA;
	
	iX=numA.create("impactX","iX",MFnNumericData::kFloat,0.0);
	iY=numA.create("impactY","iY",MFnNumericData::kFloat,0.0);
	iZ=numA.create("impactZ","iZ",MFnNumericData::kFloat,0.0);
	iR=numA.create("impactR","iR",MFnNumericData::kFloat,0.3);
	iPer=numA.create("impactPer","iP",MFnNumericData::kFloat,0.08);
	dNum=numA.create("decompNum","dNum",MFnNumericData::kInt,7);
	pNum=numA.create("patternNum","pNum",MFnNumericData::kInt,5);
	pType=numA.create("patternType","pType",MFnNumericData::kInt,1);
	inputMesh=typeA.create("inputMesh","im",MFnData::kMesh);
	outputMesh=typeA.create("outputMesh","ot",MFnData::kMesh);
	fractureMode=numA.create("fractureMode","fm",MFnNumericData::kInt,0);
	objPatternNum=numA.create("objPatternNum","objN",MFnNumericData::kInt,3);
	objPatternFile=typeA.create("objPatternFile","objF",MFnData::kString);

	addAttribute(VFracturerNode::iX);
	addAttribute(VFracturerNode::iY);
	addAttribute(VFracturerNode::iZ);
	addAttribute(VFracturerNode::iR);
	addAttribute(VFracturerNode::iPer);
	addAttribute(VFracturerNode::dNum);
	addAttribute(VFracturerNode::pNum);
	addAttribute(VFracturerNode::pType);
	addAttribute(VFracturerNode::outputMesh);
	addAttribute(VFracturerNode::inputMesh);
	addAttribute(VFracturerNode::fractureMode);
	addAttribute(VFracturerNode::objPatternNum);
	addAttribute(VFracturerNode::objPatternFile);


	//attributeAffects(VFracturerNode::iX,VFracturerNode::outputMesh);
	//attributeAffects(VFracturerNode::iY,VFracturerNode::outputMesh);
	//attributeAffects(VFracturerNode::iZ,VFracturerNode::outputMesh);
	attributeAffects(VFracturerNode::iR,VFracturerNode::outputMesh);
	//attributeAffects(VFracturerNode::iPer,VFracturerNode::outputMesh);
	attributeAffects(VFracturerNode::inputMesh,VFracturerNode::outputMesh);
	MGlobal::displayInfo("VFracturer node initialized");
	return MStatus::kSuccess;
}

MStatus VFracturerNode::compute(const MPlug & plug, MDataBlock& data)
{
	MObject pluginO=MFnPlugin::findPlugin("VFracturer");
	MFnPlugin plugin=MFnPlugin(pluginO);
	MString LP=plugin.loadPath();	
	MStatus status;

	if(plug==VFracturerNode::outputMesh)
	{
		//creating data handles
		MGlobal::displayInfo("VFractuer: input changed - Re-computing!");
		MDataHandle iXData=data.inputValue(iX); 
		MDataHandle iYData=data.inputValue(iY); 
		MDataHandle iZData=data.inputValue(iZ); 
		MDataHandle iRData=data.inputValue(iR); 
		MDataHandle iPerData=data.inputValue(iPer); 
		MDataHandle dNumData=data.inputValue(dNum); 
		MDataHandle pNumData=data.inputValue(pNum); 
		MDataHandle pTypeData=data.inputValue(pType); 
		MDataHandle outData=data.outputValue(outputMesh);
		MDataHandle inDataA=data.inputValue(inputMesh);

		MDataHandle modeData=data.inputValue(fractureMode); 
		MDataHandle objNData=data.inputValue(objPatternNum); 
		MDataHandle objFData=data.inputValue(objPatternFile); 
		//////////////////////////////////////////////////////////////////////////////////////////

		//extract data from handles
		float impX=iXData.asFloat();
		float impY=iYData.asFloat();
		float impZ=iZData.asFloat();
		float impR=iRData.asFloat();
		float impPer=iPerData.asFloat();
		int patternType = pTypeData.asInt();
		int num_pattern=pNumData.asInt();
		int num_decomp=dNumData.asInt();
		int mode=modeData.asInt();

		string objFileName("shard");
		MString temp=objFData.asString();
		if(temp.length()>0) objFileName=temp.asChar();
		int totalObjCellNum=objNData.asInt();
		//////////////////////////////////////////////////////////////////////////////////////////


		MFnMesh inMesh(inDataA.asMesh());
		Mesh mayaMesh;
		VConverter::MFnMesh_to_mesh(inMesh,mayaMesh);
		//debugging input parameters
		MGlobal::displayInfo("VFracturer: MFnMesh to Mesh converted, mesh has "+toString(mayaMesh.numVertices()));
		MGlobal::displayInfo("iX "+toString(impX)); MGlobal::displayInfo("iY "+toString(impY));MGlobal::displayInfo("iZ "+toString(impZ));MGlobal::displayInfo("iR "+toString(impR));MGlobal::displayInfo("iPer "+toString(impPer));
		MGlobal::displayInfo("dNum "+toString(num_decomp));MGlobal::displayInfo("pNum "+toString(num_pattern));

		//initialize empty mesh data for writing
		MFnMeshData meshCreator;
		MObject newOutput=meshCreator.create();
		MPointArray points;
		points.clear();
		MIntArray faceCounts;faceCounts.clear();
		MIntArray faceConnects;faceConnects.clear(); 
		//////////////////////////////////////////////////////////////////////////////////////////


		//initialize states and containers
		vec3 impactPos(impX,impY,impZ);
		float r=impR;
		float closenessPercent=impPer;
		std::vector<Polyhedron_3 *> cells; //TODO RELEASE MEMORY
		std::vector<Polyhedron_3 *> otherCells;//TODO RELEASE MEMORY
		std::vector<Polyhedron_3 *> patternCells; //TODO RELEASE MEMORY
		std::vector<Polyhedron_3 *> result;  //TODO RELEASE MEMORY
		Polyhedron_3 * impactCell;
		Polyhedron_3 * P=new Polyhedron_3();
		VConverter::mesh_to_CGAL_Polyhedron_3(mayaMesh,*P); float x,X,y,Y,z,Z;
		dimensionPolyhedron_3(*P,x,X,y,Y,z,Z);
		//////////////////////////////////////////////////////////////////////////////////////////

		//check impactpos valid
		if(!withinDimension(impactPos,vec3(X,Y,Z),vec3(x,y,z)) ) 
		{
			clipDimension(impactPos,vec3(X,Y,Z),vec3(x,y,z));
			MGlobal::displayInfo("VFracturer::impact position outside of object! position auto clipped to within dimension");
			//return MStatus(MStatus::kFailure);
		}

	
		VClipper clipper=VClipper(*P);

		//global decomp mode
		if(mode==0)
		{
			//genVoroCell(cells, otherCells, impactCell,num_decomp,vec3(X,Y,Z),vec3(x,y,z),impactPos,r,patternType);	
			genRegularVoroCell(cells,  num_decomp, vec3(X,Y,Z),vec3(x,y,z));
		
			//load each cell to the clipper
			for(int i = 0; i < cells.size(); ++i)
			{
				clipper.loadDecompositionCell(*cells[i]);
			}

//			clipper.loadDecompositionCell(*impactCell); 
			clipper.decomp(result);
		}

		//local impact mode
		else if(mode==1)
		{
			genVoroCell(cells, otherCells, impactCell,num_decomp,vec3(X,Y,Z),vec3(x,y,z),impactPos,r,patternType);	
		
			//load each cell to the clipper
			for(int i = 0; i < otherCells.size(); ++i)
			{
				clipper.loadDecompositionCell(*otherCells[i]);
			}

			clipper.loadImpactCell(*impactCell); dimensionPolyhedron_3(*impactCell,x,X,y,Y,z,Z);
			genVoroPatternCell(patternCells, num_pattern,vec3(X,Y,Z),vec3(x,y,z),impactPos,closenessPercent,patternType);	
			for(int i = 0; i < patternCells.size(); ++i)
			{
				clipper.loadPatternCell(*patternCells[i]);
			}
			clipper.output(result);
		}

		//obj pattern mode
		else if (mode == 2)
		{
			//generate our geometries
			for (int i=1;i<totalObjCellNum+1;++i)
			{
				string cellName=objFileName+std::to_string(static_cast<long long>(i))+".obj";
				generateMeshFromObj(cellName);///////////////////////TESTING
			}

			//find object center and scale for calibrating pattern cell
			float cx,cy,cz,sx,sy,sz;
			cx=(x+X)/2.0f;cy=(y+Y)/2.0f;cz=(z+Z)/2.0f;
			sx=abs(x-X);sy=abs(y-Y);sz=abs(z-Z);

			for (int i=0;i<objPatternCellList.size();++i)
			{
				Polyhedron_3 * cell=new Polyhedron_3();
				VConverter::mesh_to_CGAL_Polyhedron_3(objPatternCellList.at(i),*cell); 
				translate(*cell,cx,cy,cz);
				scale(*cell,sx,sy,sz);
				clipper.loadDecompositionCell(*cell);
				delete cell;
			}
			clipper.decomp(result);
		}
		//////////////////////////////////////////////////////////////////////////////////////////

		//load our mesh to MfnMesh for outputMesh of the node
		for(int j=0;j<result.size();++j)
		{
			Mesh mMesh;
			int offset=points.length();
			//MGlobal::displayInfo("offset is"+toString(offset));
			
			VConverter::CGAL_Polyhedron_3_to_Mesh(*(result.at(j)),mMesh);
				//load vertices
			for(int i=0;i<mMesh.numVertices();++i)
			{
				vec3 vertex=mMesh.getVertex(i);
				points.append(vertex[0],vertex[1],vertex[2]);
			}
			//load face
			 for(int i=0;i<mMesh.faceList.size();i+=3)
			 {
				 faceCounts.append(3);
				 faceConnects.append(mMesh.faceList.at(i)+offset);
				 faceConnects.append(mMesh.faceList.at(i+1)+offset);
				 faceConnects.append(mMesh.faceList.at(i+2)+offset);
			
			 }

		}

	
		MFnMesh outMesh;

		outMesh.create(points.length(),faceCounts.length(), points,faceCounts,faceConnects,newOutput,&status);
		outData.set(newOutput);
		data.setClean(plug);

		//release resource
		delete P;
		for(int i=0;i<result.size();++i)
			delete result.at(i);
		
	}


	else
		return MS::kUnknownParameter;
	return status;
}


MString toString(int i)
{
	char buffer[1024];
	sprintf_s(buffer,1024,"%d",i);
	return buffer;
}

MString toString(double f)
{
	char buffer[1024];
	sprintf_s(buffer,1024,"%f",f);
	return buffer;
}

void generateMeshFromObj(string file)
{
	int faceCount=0;
	ifstream File;
	string line;
	string tag;
	vector<string> para;
	File.open(file);
	MGlobal::displayInfo("obj file opened start parsing");
	
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
	MGlobal::displayInfo("obj file loading completed");
	objPatternCellList.push_back(m);
	File.close();
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