/**
**  by Zhen Gou, Zhenghan Mei 2014
**/
#include <maya/MPxCommand.h>
#include <maya/MIOStream.h>
#include <maya/MString.h>
#include <maya/MArgList.h>
#include <maya/MSimple.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MDGModifier.h>
#include <maya/MPlugArray.h>
#include <maya/MVector.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MStringArray.h>
#include <maya/MFnPlugin.h>
#include <list>

#include "VFracturerNode.h"

MStatus initializePlugin( MObject obj )
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj, "VFracturerPlugin", "1.0", "Zhen Gou, Zhenghan Mei");
    // Register Command
	
  
	status=plugin.registerNode("VFracturerNode",VFracturerNode::id,VFracturerNode::creator,VFracturerNode::initialize);
	if(!status)
	{
		status.perror("register node failed");
	}
	MString path=MString("\"")+plugin.loadPath()+MString("/VFracturerUI.mel")+MString("\"");
	MGlobal::executeCommand("source "+path, true);
	//plugin.registerUI("LSystem","deleteMyUI");
    return status;
}

MStatus uninitializePlugin( MObject obj)
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj );

	status=plugin.deregisterNode(VFracturerNode::id);
	if (!status) {
	    status.perror("deregister node failed");
	    return status;
    }

    return status;
}


