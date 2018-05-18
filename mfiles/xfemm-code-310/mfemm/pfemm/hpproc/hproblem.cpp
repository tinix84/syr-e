/*
   This code is a modified version of an algorithm
   forming part of the software program Finite
   Element Method Magnetics (FEMM), authored by
   David Meeker. The original software code is
   subject to the Aladdin Free Public Licence
   version 8, November 18, 1999. For more information
   on FEMM see www.femm.info. This modified version
   is not endorsed in any way by the original
   authors of FEMM.

   This software has been modified to use the C++
   standard template libraries and remove all Microsoft (TM)
   MFC dependent code to allow easier reuse across
   multiple operating system platforms.

   Date Modified: 2014 - 03 - 21
   By:  Emoke Szelitzky
        Tibor Szelitzky
        Richard Crozier
   Contact:
        szelitzkye@gmail.com
        sztibi82@gmail.com
        richard.crozier@yahoo.co.uk
*/

#include <cstdlib>
#include <cmath>
#include "hproblem.h"
#include "fullmatrix.h"
#include "femmcomplex.h"

#ifndef muo
#define muo 1.2566370614359173e-6
#endif

#define ElementsPerSkinDepth 10

using namespace std;
/////////////////////////////////////////////////////////////////////////////
// CNode construction

CNode::CNode()
{
	x=0.;
	y=0.;
	IsSelected=false;
	BoundaryMarker=-1;
	InConductor=-1;
	InGroup=0;
}

CComplex CNode::CC()
{
	return (x+I*y);
}

double CNode::GetDistance(double xo, double yo)
{
	return sqrt((x-xo)*(x-xo) + (y-yo)*(y-yo));
}

void CNode::ToggleSelect()
{
	if (IsSelected==true) IsSelected=false;
	else IsSelected=true;
}

/////////////////////////////////////////////////////////////////////////////
// CMeshNode construction

CMeshNode::CMeshNode()
{
	x=y=T=0;
	Q=false;
	msk=0;
	IsSelected=false;
}

CComplex CMeshNode::CC()
{
	return (x+I*y);
}

double CMeshNode::GetDistance(double xo, double yo)
{
	return sqrt((x-xo)*(x-xo) + (y-yo)*(y-yo));
}

/////////////////////////////////////////////////////////////////////////////
// CSegment construction

CSegment::CSegment()
{
	n0=0;
	n1=0;
	IsSelected=false;
	BoundaryMarker=-1;
	InConductor=-1;
	InGroup=0;
}

void CSegment::ToggleSelect()
{
	if (IsSelected==true) IsSelected=false;
	else IsSelected=true;
}
/////////////////////////////////////////////////////////////////////////////
// CArcSegment construction

CArcSegment::CArcSegment()
{
	n0=0;
	n1=0;
	IsSelected=false;
	MaxSideLength=-1;
	ArcLength=90.;
	BoundaryMarker=-1;
	InConductor=-1;
	InGroup=0;
}

void CArcSegment::ToggleSelect()
{
	if (IsSelected==true) IsSelected=false;
	else IsSelected=true;
}

/////////////////////////////////////////////////////////////////////////////
// CNode construction

CBlockLabel::CBlockLabel()
{
	x=0.;
	y=0.;
	MaxArea=0.;
	IsSelected=false;
	InGroup=0;
	BlockType=-1;
	IsExternal=false;
	IsDefault=false;
}

void CBlockLabel::ToggleSelect()
{
	if (IsSelected==true) IsSelected=false;
	else IsSelected=true;
}

double CBlockLabel::GetDistance(double xo, double yo)
{
	return sqrt((x-xo)*(x-xo) + (y-yo)*(y-yo));
}

CMaterialProp::CMaterialProp()
{
		BlockName="New Material";
		Kx=Ky=1;
		Kt=0;
		qv=0;
}

CBoundaryProp::CBoundaryProp()
{
		BdryName="New Boundary";
		BdryFormat=0;
		Tset=Tinf=h=beta=qs=0;
		InConductor="<None>";
}

CComplex CMaterialProp::GetK(double t)
{
	int i,j;

	// Kx returned as real part;
	// Ky returned as imag part

	if (npts==0) return (Kx+I*Ky);
	if (npts==1) return (Im(Kn[0])*(1+I));
	if (t<=Re(Kn[0])) return (Im(Kn[0])*(1+I));
	if (t>=Re(Kn[npts-1])) return (Im(Kn[npts-1])*(1+I));

	for(i=0,j=1;j<npts;i++,j++)
	{
		if((t>=Re(Kn[i])) && (t<=Re(Kn[j])))
		{
			return (1+I)*(Im(Kn[i])+Im(Kn[j]-Kn[i])*Re(t-Kn[i])/Re(Kn[j]-Kn[i]));
		}
	}

	return (Kx+I*Ky);
}

CPointProp::CPointProp()
{
		PointName="New Point Property";
		InConductor="<None>";
		V=qp=0;
}

CCircuit::CCircuit()
{
		CircName="New Circuit";
		V=q=0;
		CircType=0;
}
