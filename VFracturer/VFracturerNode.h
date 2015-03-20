/**
**  by Zhen Gou, Zhenghan Mei 2014
**/
#ifndef VFRACTURERNODE_H
#define VFRACTURERNODE_H
#include "maya/MPxNode.h"
#include "maya/MFnNumericAttribute.h"
#include "maya/MFnTypedAttribute.h"
#include "maya/MFnUnitAttribute.h"
#include "maya/MTypeId.h"
#include "maya/MFnNumericData.h"
#include "maya/MFnStringData.h"
#include "maya/MFnMeshData.h"
#include "maya/MTime.h"
#include "maya/MGlobal.h"
#include "maya/MFnMesh.h"
#define MNoPluginEntry
#define MNoVersionString
#include "maya/MFnPlugin.h"
#include <maya/MPoint.h>
#include <maya/MPointArray.h>

#include <maya/MVector.h>
#include <maya/MVectorArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>

class VFracturerNode:public MPxNode
{
public:
	VFracturerNode();
	virtual ~VFracturerNode();
	virtual MStatus compute(const MPlug & plug, MDataBlock& data);
	static void *creator();
	static MStatus initialize();

	static MObject iX;
	static MObject iY;
	static MObject iZ;
	static MObject iR;
	static MObject iPer;
	static MObject dNum;
	static MObject pNum;
	static MObject pType;
	static MObject outputMesh;
	static MObject inputMesh;
	static MObject fractureMode; //0 for no local//1 for local// 2 for obj
	static MObject objPatternNum;
	static MObject objPatternFile;

	static MTypeId id;

};

#endif