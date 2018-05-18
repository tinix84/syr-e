// fpproc.cpp : implementation of the FPProc class
//

// #include "stdafx.h"
// #include <afx.h>
// #include <afxtempl.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
// include the boost format lib to get nicer string handling capabilies
//#include "boost/format.hpp"
#include "problem.h"
#include "femmcomplex.h"
#include "fparse.h"
//#include "femm.h"
//#include "xyplot.h"
//#include "femmviewView.h"
#include "lua.h"
#include "lualib.h"
#include "fpproc.h"

#ifndef muo
#define muo 1.2566370614359173e-6
#endif

#ifndef Golden
#define Golden 0.3819660112501051517954131656
#endif

//#ifndef PLANAR
//#define PLANAR 0
//#endif
//
//#ifndef AXISYMMETRIC
//#define AXISYMMETRIC 1
//#endif

#ifndef _MSC_VER
#define _strnicmp strncasecmp
#endif

extern lua_State * lua;
extern void *pFemmviewdoc;
//extern CLuaConsoleDlg *LuaConsole;
extern bool bLinehook;

extern void lua_baselibopen (lua_State *L);
extern void lua_iolibopen (lua_State *L);
extern void lua_strlibopen (lua_State *L);
extern void lua_mathlibopen (lua_State *L);
extern void lua_dblibopen (lua_State *L);

// char *ParseDbl(char *t, double *f);
// char *ParseString(char *t, CString *s);
// char *ParseInt(char *t, int *f);

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//using namespace std;
using namespace femm;

// FPProc construction/destruction

double sqr(double x)
{
    return x*x;
}

// Replacement function for windows MFC function AfxMessageBox
//void FPProc::MsgBox(const char* message)
//{
//    printf(message);
//}

/**
 * Constuctor for the FPProc class.
 */
FPProc::FPProc()
{
    // set some default values for problem definition
    d_LineIntegralPoints = 400;
    d_ShiftH = true;
    Frequency = 0.;
    Depth = 1/0.0254;
    LengthUnits = 0;
    ProblemType = PLANAR;
    ProblemNote = "Add comments here.";
    FirstDraw = -1;
    A_High = 0.;
    A_Low = 0.;
    A_lb = 0.;
    A_ub = 0.;
    extRo = extRi = extZo = 0;
    Smooth = true;
    NumList = NULL;
    ConList = NULL;
    bHasMask = false;
    LengthConv = (double *)calloc(6,sizeof(double));
    LengthConv[0] = 0.0254;   //inches
    LengthConv[1] = 0.001;    //millimeters
    LengthConv[2] = 0.01;     //centimeters
    LengthConv[3] = 1.;       //meters
    LengthConv[4] = 2.54e-05; //mils
    LengthConv[5] = 1.e-06;   //micrometers
    Coords = CART;

    for(int i=0; i<9; i++)
    {
        d_PlotBounds[i][0] = d_PlotBounds[i][1] =
                                 PlotBounds[i][0] = PlotBounds[i][1] = 0;
    }

    // initialise the warning message function pointer to
    // point to the PrintWarningMsg function
    WarnMessage = &PrintWarningMsg;

    // determine path to bin directory
//     BinDir=((CFemmApp *)AfxGetApp())->GetExecutablePath();

//     ScanPreferences();

    // lua initialization stuff
    // initalise_lua();
}

/**
 * Destructor for the FPProc class.
 */
FPProc::~FPProc()
{
    unsigned int i;

    free(LengthConv);
    for(i=0; i<meshnode.size(); i++)
        if(ConList[i]!=NULL) free(ConList[i]);
    free(ConList);
    free(NumList);

    //if (pFemmviewdoc==this) pFemmviewdoc = NULL;
}

/**
 * Clear out all data associated with the last document to be loaded
 *
 */
void FPProc::ClearDocument()
{

    // clear out all current lines, nodes, and block labels
    if(NumList!=NULL)
    {
        for(unsigned int i=0; i<meshnode.size(); i++)
        {
            if(ConList[i]!=NULL)
            {
                free(ConList[i]);
            }
        }
        free(ConList);
        ConList = NULL;
        free(NumList);
        NumList = NULL;
    }
    nodelist.clear();
    linelist.clear();
    blocklist.clear();
    arclist.clear();
    nodeproplist.clear();
    lineproplist.clear();
    blockproplist.clear();
    circproplist.clear();
    meshnode.clear();
    meshelem.clear();
    contour.clear();

}

/**
 * Performs actions required when a new document is to be
 * loaded, including clearing out all existing data and
 * resetting various values to their defaults.
 */
bool FPProc::NewDocument()
{
    // TODO: add reinitialization code here
    // (SDI documents will reuse this document)

    // clear out all current lines, nodes, and block labels
    ClearDocument();

    // set problem attributes to generic ones;
    Frequency = 0.;
    LengthUnits = 0;
    Precision = 1e-8;
    ProblemType = PLANAR;
    ProblemNote = "Add comments here.";
    bHasMask = false;
    extRo = extRi = extZo = 0;
    Depth = -1;

    return true;
}

/////////////////////////////////////////////////////////////////////////////
// FPProc serialization

// void FPProc::Serialize(CArchive& ar)
// {
//     if (ar.IsStoring())
//     {
//         // TODO: add storing code here
//     }
//     else
//     {
//         // TODO: add loading code here
//     }
// }

/////////////////////////////////////////////////////////////////////////////
// FPProc diagnostics

#ifdef _DEBUG
void FPProc::AssertValid() const
{
    CDocument::AssertValid();
}

void FPProc::Dump(CDumpContext& dc) const
{
    CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// FPProc commands


bool FPProc::OpenDocument(string pathname)
{

    FILE *fp;
    int i,j,k,t;
    char s[1024],q[1024];
    char *v;
    double b,bi,br;
    double zr,zi;
    bool flag = false;
    CPointProp    PProp;
    CBoundaryProp BProp;
    CMaterialProp MProp;
    CCircuit      CProp;
    CNode         node;
    CSegment      segm;
    CArcSegment   asegm;
    CElement      elm;
    CBlockLabel   blk;
    CMeshNode     mnode;
    //CPoint        mline;

    // clear out all the document data and set defaults to standard values
    NewDocument();

    // attempt to open the file for reading
    if ((fp = fopen(pathname.c_str(),"rt")) == NULL)
    {
        WarnMessage("Couldn't read from specified .ans file");
        return false;
    }

    // parse the file
    while ((flag==false) && (fgets(s,1024,fp) != NULL))
    {
        sscanf(s,"%s",q);

        // Deal with flag for file format version
        if( _strnicmp(q,"[format]",8)==0 )
        {
            double vers;
            v=StripKey(s);
            sscanf(v,"%lf",&vers);
            vers = 10.*vers + 0.5;
            if( ((int) vers)!=40 )
            {
                WarnMessage("This file is from a different version of FEMM\nRe-analyze the problem using the current version.");
                fclose(fp);
                return false;
            }
            q[0] = '\0';
        }

        // Frequency of the problem
        if( _strnicmp(q,"[frequency]",11)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&Frequency);
            q[0] = '\0';
        }

        // Frequency of the problem
        if( _strnicmp(q,"[depth]",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&Depth);
            q[0] = '\0';
        }

        // Precision
        if( _strnicmp(q,"[precision]",11)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&Precision);
            q[0] = '\0';
        }

        // Units of length used by the problem
        if( _strnicmp(q,"[lengthunits]",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%s",q);
            if( _strnicmp(q,"inches",6)==0) LengthUnits=0;
            else if( _strnicmp(q,"millimeters",11)==0) LengthUnits=1;
            else if( _strnicmp(q,"centimeters",1)==0) LengthUnits=2;
            else if( _strnicmp(q,"mils",4)==0) LengthUnits=4;
            else if( _strnicmp(q,"microns",6)==0) LengthUnits=5;
            else if( _strnicmp(q,"meters",6)==0) LengthUnits=3;
            q[0] = '\0';
        }

        // Problem Type (planar or axisymmetric)
        if( _strnicmp(q,"[problemtype]",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%s",q);
            if( _strnicmp(q,"planar",6)==0) ProblemType=PLANAR;
            if( _strnicmp(q,"axisymmetric",3)==0) ProblemType=AXISYMMETRIC;
            q[0] = '\0';
        }

        // Coordinates (cartesian or polar)
        if( _strnicmp(q,"[coordinates]",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%s",q);
            if ( _strnicmp(q,"cartesian",4)==0) Coords=CART;
            if ( _strnicmp(q,"polar",5)==0) Coords=POLAR;
            q[0] = '\0';
        }

        // Comments
        if (_strnicmp(q,"[comment]",9)==0)
        {
            v=StripKey(s);
            // put in carriage returns;
            k=std::strlen(v);
            for(i=0; i<k; i++)
                if((v[i]=='\\') && (v[i+1]=='n'))
                {
                    v[i]=13;
                    v[i+1]=10;
                }

            for(i=0; i<k; i++)
                if(v[i]=='\"')
                {
                    v=v+i+1;
                    i=k;
                }
            k=std::strlen(v);
            if(k>0) for(i=k-1; i>=0; i--)
                {
                    if(v[i]=='\"')
                    {
                        v[i]=0;
                        i=-1;
                    }
                }
            ProblemNote=v;
            q[0] = '\0';
        }

        // properties for axisymmetric external region
        if( _strnicmp(q,"[extzo]",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&extZo);
            q[0] = '\0';
        }

        if( _strnicmp(q,"[extro]",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&extRo);
            q[0] = '\0';
        }

        if( _strnicmp(q,"[extri]",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&extRi);
            q[0] = '\0';
        }

        // Point Properties
        if( _strnicmp(q,"<beginpoint>",11)==0)
        {
            PProp.PointName="New Point Property";
            PProp.Jr=0.;
            PProp.Ji=0.;
            PProp.Ar=0.;
            PProp.Ai=0.;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<pointname>",11)==0)
        {
            v=StripKey(s);
            k=std::strlen(v);
            for(i=0; i<k; i++)
                if(v[i]=='\"')
                {
                    v=v+i+1;
                    i=k;
                }
            k=std::strlen(v);
            if(k>0) for(i=k-1; i>=0; i--)
                {
                    if(v[i]=='\"')
                    {
                        v[i]=0;
                        i=-1;
                    }
                }
            PProp.PointName=v;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<A_re>",6)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&PProp.Ar);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<A_im>",6)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&PProp.Ai);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<I_re>",6)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&PProp.Jr);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<I_im>",6)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&PProp.Ji);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<endpoint>",9)==0)
        {
            nodeproplist.push_back(PProp);
            q[0] = '\0';
        }

        // Boundary Properties;
        if( _strnicmp(q,"<beginbdry>",11)==0)
        {
            BProp.BdryName = "New Boundary";
            BProp.BdryFormat = 0;
            BProp.A0  = 0.;
            BProp.A1  = 0.;
            BProp.A2  = 0.;
            BProp.phi = 0.;
            BProp.Mu  = 0.;
            BProp.Sig = 0.;
            BProp.c0  = 0.;
            BProp.c1  = 0.;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<bdryname>",10)==0)
        {
            v=StripKey(s);
            k=std::strlen(v);
            for(i=0; i<k; i++)
                if(v[i]=='\"')
                {
                    v=v+i+1;
                    i=k;
                }
            k=std::strlen(v);
            if(k>0) for(i=k-1; i>=0; i--)
                {
                    if(v[i]=='\"')
                    {
                        v[i]=0;
                        i=-1;
                    }
                }
            BProp.BdryName = v;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<bdrytype>",10)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&BProp.BdryFormat);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<mu_ssd>",8)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.Mu);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<sigma_ssd>",11)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.Sig);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<A_0>",5)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.A0);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<A_1>",5)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.A1);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<A_2>",5)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.A2);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<phi>",5)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.phi);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<c0>",4)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.c0.re);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<c1>",4)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.c1.re);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<c0i>",5)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.c0.im);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<c1i>",5)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&BProp.c1.im);
            q[0] = '\0';
        }
        if( _strnicmp(q,"<endbdry>",9)==0)
        {
            lineproplist.push_back(BProp);
            q[0] = '\0';
        }


        // Block Properties;
        if( _strnicmp(q,"<beginblock>",12)==0)
        {
            MProp.BlockName="New Material";
            MProp.mu_x=1.;
            MProp.mu_y=1.;            // permeabilities, relative
            MProp.H_c=0.;                // magnetization, A/m
            MProp.Jr=0.;
            MProp.Ji=0.;                // applied current density, MA/m^2
            MProp.Cduct=0.;            // conductivity of the material, MS/m
            MProp.Lam_d=0.;            // lamination thickness, mm
            MProp.Theta_hn=0.;            // hysteresis angle, degrees
            MProp.Theta_hx=0.;            // hysteresis angle, degrees
            MProp.Theta_hy=0.;            // hysteresis angle, degrees
            MProp.NStrands=0;
            MProp.WireD=0;
            MProp.LamFill=1.;            // lamination fill factor;
            MProp.LamType=0;            // type of lamination;
            MProp.BHpoints=0;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<blockname>",10)==0)
        {
            v=StripKey(s);
            k=std::strlen(v);
            for(i=0; i<k; i++)
                if(v[i]=='\"')
                {
                    v=v+i+1;
                    i=k;
                }
            k=std::strlen(v);
            if(k>0) for(i=k-1; i>=0; i--)
                {
                    if(v[i]=='\"')
                    {
                        v[i]=0;
                        i=-1;
                    }
                }
            MProp.BlockName=v;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<mu_x>",6)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.mu_x);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<mu_y>",6)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.mu_y);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<H_c>",5)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.H_c);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<J_re>",6)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.Jr);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<J_im>",6)==0)
        {
            v=StripKey(s);
            if (Frequency!=0) sscanf(v,"%lf",&MProp.Ji);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<sigma>",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.Cduct);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<phi_h>",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.Theta_hn);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<phi_hx>",8)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.Theta_hx);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<phi_hy>",8)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.Theta_hy);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<d_lam>",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.Lam_d);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<LamFill>",8)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.LamFill);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<LamType>",9)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&MProp.LamType);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<NStrands>",10)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&MProp.NStrands);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<WireD>",7)==0)
        {
            v=StripKey(s);
            sscanf(v,"%lf",&MProp.WireD);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<BHPoints>",10)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&MProp.BHpoints);
            if (MProp.BHpoints>0)
            {
                //MProp.Hdata = (CComplex *)calloc(MProp.BHpoints,sizeof(CComplex));
                //MProp.Bdata =   (double *)calloc(MProp.BHpoints,sizeof(double));
                MProp.Hdata = (CComplex *)realloc(MProp.Hdata,MProp.BHpoints*sizeof(CComplex));
                MProp.Bdata = (double *)realloc(MProp.Bdata,MProp.BHpoints*sizeof(double));
                // reset the contents of the B and H data arrays to all zeros
                memset(MProp.Hdata, 0, MProp.BHpoints * sizeof(CComplex));
                memset(MProp.Bdata, 0, MProp.BHpoints * sizeof(double));
                for(j=0; j<MProp.BHpoints; j++)
                {
                    fgets(s,1024,fp);
                    MProp.Hdata[j]=0;
                    sscanf(s,"%lf\t%lf",&MProp.Bdata[j],&MProp.Hdata[j].re);
                }
            }
            q[0] = '\0';
        }

        if( _strnicmp(q,"<endblock>",9)==0)
        {
            if (MProp.BHpoints>0)
            {
                MProp.GetSlopes(Frequency*2.*PI);
            }

            blockproplist.push_back(MProp);

            // reinitialise the material property, and free allocated memory
            MProp.BHpoints=0;
            if (MProp.Bdata != NULL) free(MProp.Bdata);
            MProp.Bdata = NULL;
            if (MProp.Hdata != NULL) free(MProp.Hdata);
            MProp.Hdata = NULL;
            if (MProp.slope != NULL) free(MProp.slope);
            MProp.slope = NULL;
            q[0] = '\0';
        }

        // Circuit Properties
        if( _strnicmp(q,"<begincircuit>",14)==0)
        {
            CProp.CircName="New Circuit";
            CProp.CircType=0;
            CProp.Amps=0;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<circuitname>",13)==0)
        {
            v=StripKey(s);
            k=std::strlen(v);
            for(i=0; i<k; i++)
                if(v[i]=='\"')
                {
                    v=v+i+1;
                    i=k;
                }
            k=std::strlen(v);
            if(k>0) for(i=k-1; i>=0; i--)
                {
                    if(v[i]=='\"')
                    {
                        v[i]=0;
                        i=-1;
                    }
                }
            CProp.CircName=v;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<totalamps_re>",14)==0)
        {
            double inval;
            v=StripKey(s);
            sscanf(v,"%lf",&inval);
            CProp.Amps+=inval;
            q[0] = '\0';
        }

        if( _strnicmp(q,"<totalamps_im>",14)==0)
        {
            double inval;
            v=StripKey(s);
            sscanf(v,"%lf",&inval);
            if (Frequency!=0) CProp.Amps+=(I*inval);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<circuittype>",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&CProp.CircType);
            q[0] = '\0';
        }

        if( _strnicmp(q,"<endcircuit>",12)==0)
        {
            circproplist.push_back(CProp);
            q[0] = '\0';
        }

        // Points list;
        if(_strnicmp(q,"[numpoints]",11)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&k);
            for(i=0; i<k; i++)
            {
                fgets(s,1024,fp);
                sscanf(s,"%lf\t%lf\t%i\n",&node.x,&node.y,&t);
                node.BoundaryMarker=t-1;
                nodelist.push_back(node);
            }
            q[0] = '\0';
        }

        // read in segment list
        if(_strnicmp(q,"[numsegments]",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&k);
            for(i=0; i<k; i++)
            {
                int hidden = 0;
                fgets(s,1024,fp);
                sscanf(s,"%i\t%i\t%lf %i\t%i\n",&segm.n0,&segm.n1,&segm.MaxSideLength,&t,&hidden);
                segm.BoundaryMarker = t-1;
                if (hidden == 0)
                {
                    segm.Hidden = false;
                }
                else
                {
                    segm.Hidden = true;
                }
                linelist.push_back(segm);
            }
            q[0] = '\0';
        }

        // read in arc segment list
        if(_strnicmp(q,"[numarcsegments]",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&k);
            for(i=0; i<k; i++)
            {
                int hidden = 0;
                fgets(s,1024,fp);
                sscanf(s,"%i\t%i\t%lf\t%lf %i\t%i\n",&asegm.n0,&asegm.n1,
                       &asegm.ArcLength,&asegm.MaxSideLength,&t,&hidden);
                asegm.BoundaryMarker=t-1;
                if (hidden == 0)
                {
                    asegm.Hidden = false;
                }
                else
                {
                    asegm.Hidden = true;
                }
                arclist.push_back(asegm);
            }
            q[0] = '\0';
        }


        // read in list of holes;
        if(_strnicmp(q,"[numholes]",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&k);
            if(k>0)
            {
                blk.BlockType=-1;
                blk.MaxArea=0;
                for(i=0; i<k; i++)
                {
                    fgets(s,1024,fp);
                    sscanf(s,"%lf\t%lf\n",&blk.x,&blk.y);
                    //    blocklist.push_back(blk);
                    //  don't add holes to the list
                    //  of block labels because it messes up the
                    //  number of block labels.
                }
            }
            q[0] = '\0';
        }

        // read in regional attributes
        if(_strnicmp(q,"[numblocklabels]",13)==0)
        {
            v=StripKey(s);
            sscanf(v,"%i",&k);
            for(i=0; i<k; i++)
            {
                fgets(s,1024,fp);

                //some defaults
                blk.MaxArea=0.;
                blk.MagDir=0.;
                blk.MagDirFctn.clear();
                blk.Turns=1;
                blk.InCircuit=0;
                blk.InGroup=0;
                int isexternal = 0;
                blk.IsExternal=false;

                // scan in data
                v=ParseDbl(s,&blk.x);
                v=ParseDbl(v,&blk.y);
                v=ParseInt(v,&blk.BlockType);
                v=ParseDbl(v,&blk.MaxArea);
                v=ParseInt(v,&blk.InCircuit);
                v=ParseDbl(v,&blk.MagDir);
                v=ParseInt(v,&blk.InGroup);
                v=ParseInt(v,&blk.Turns);
                v=ParseInt(v,&isexternal);

                if (isexternal == 0)
                {
                    blk.IsExternal = false;
                }
                else
                {
                    blk.IsExternal = true;
                }

                v=ParseString(v,&blk.MagDirFctn);

                if (blk.MaxArea<0) blk.MaxArea=0;
                else blk.MaxArea=PI*blk.MaxArea*blk.MaxArea/4.;
                blk.BlockType-=1;
                blk.InCircuit-=1;
                blocklist.push_back(blk);
            }
            q[0] = '\0';
        }

        if(_strnicmp(q,"[solution]",10)==0)
        {
            flag = true;
            q[0] = '\0';
        }
    }

    // ensure memory is freed
    if (MProp.Bdata != NULL) free(MProp.Bdata);
    if (MProp.Hdata != NULL) free(MProp.Hdata);
    if (MProp.slope != NULL) free(MProp.slope);

    if (flag == false)
    {
        // The flag was never set to true during the while loop.
        // This means the "[solution]" string was never
        // encountered
        if(feof(fp))
        {
            // We read in the whole file but never found the start of
            // a solution section
            WarnMessage("No solution found in file.\n"); /* EOF */
        }
        else if(ferror(fp))
        {
            // There was some read error while trying to read the file
            WarnMessage("An error occured while reading file.\n"); /* Error */
        }
        fclose(fp);
        return false;
    }

    // read in meshnodes;
    fscanf(fp,"%i\n",&k);
    meshnode.resize(k);
    for(i=0; i<k; i++)
    {
        fgets(s,1024,fp);
        if (Frequency!=0)
            sscanf(s,"%lf\t%lf\t%lf\t%lf",&mnode.x,&mnode.y,
                   &mnode.A.re,&mnode.A.im);
        else
        {
            sscanf(s,"%lf\t%lf\t%lf",&mnode.x,&mnode.y,&mnode.A.re);
            mnode.A.im=0;
        }
        meshnode[i] = mnode;
    }

    // read in elements;
    fscanf(fp,"%i\n",&k);
    meshelem.resize(k);
    for(i=0; i<k; i++)
    {
        fgets(s,1024,fp);
        sscanf(s,"%i\t%i\t%i\t%i",&elm.p[0],&elm.p[1],&elm.p[2],&elm.lbl);
        elm.blk=blocklist[elm.lbl].BlockType;
        meshelem[i] = elm;
    }

    // read in circuit data;
    fscanf(fp,"%i\n",&k);
    for(i=0; i<k; i++)
    {
        fgets(s,1024,fp);
        if (Frequency==0)
        {
            sscanf(s,"%i\t%lf",&j,&zr);
            blocklist[i].Case=j;
            if (j==0) blocklist[i].dVolts=zr;
            else blocklist[i].J=zr;
        }
        else
        {
            sscanf(s,"%i\t%lf\t%lf",&j,&zr,&zi);
            blocklist[i].Case=j;
            if (j==0) blocklist[i].dVolts=zr + I*zi;
            else blocklist[i].J=zr + I*zi;
        }
    }

    fclose(fp);

    // scale depth to meters for internal computations;
    if(Depth==-1) Depth=1;
    else Depth*=LengthConv[LengthUnits];

    // element centroids and radii;
    for(i=0; i<(int)meshelem.size(); i++)
    {
        meshelem[i].ctr=Ctr(i);
        for(j=0,meshelem[i].rsqr=0; j<3; j++)
        {
            b=sqr(meshnode[meshelem[i].p[j]].x-meshelem[i].ctr.re)+
              sqr(meshnode[meshelem[i].p[j]].y-meshelem[i].ctr.im);
            if(b>meshelem[i].rsqr) meshelem[i].rsqr=b;
        }
    }

    // Compute magnetization direction in each element
    lua_State *LocalLua = lua_open(4096);
    lua_baselibopen(LocalLua);
    lua_strlibopen(LocalLua);
    lua_mathlibopen(LocalLua);
    char magbuff[4096];
    // create the formatter object in case of a lua defined mag direction
    //boost::format fmatter("x=%.17g\ny=%.17g\nr=x\nz=y\ntheta=%.17g\nR=%.17g\nreturn %s");
    for(i=0; i<(int)meshelem.size(); i++)
    {
        if(blocklist[meshelem[i].lbl].MagDirFctn.length()==0)
        {
            // The magnetisation direction in the associated blockis just
            // a number, so store it for this element
            meshelem[i].magdir = blocklist[meshelem[i].lbl].MagDir;
        }
        else
        {
            // The magnetization direction is defined by a lua calculation
            string str;
            CComplex X;

            // Get the element centroid
            X = meshelem[i].ctr;
            // generate the string using boost::format
//            fmatter % (X.re)
//            % (X.im)
//            % (arg(X)*180/PI)
//            % (abs(X))
//            % (blocklist[meshelem[i].lbl].MagDirFctn);
//            // get the created string
//            str = fmatter.str();
            // generate the string using snprintf
            std::snprintf(magbuff, 4096, "x=%.17g\ny=%.17g\nr=x\nz=y\ntheta=%.17g\nR=%.17g\nreturn %s",
                          (X.re) , (X.im) , (arg(X)*180/PI) , (abs(X)) , (blocklist[meshelem[i].lbl].MagDirFctn.c_str() ) );
            str = magbuff;
            // Have the lua interpreter evaluate the string
            lua_dostring(LocalLua,str.c_str());
            // Put the last number produced by lua into the element mag
            // direction
            meshelem[i].magdir = Re(lua_tonumber(LocalLua,-1));

            lua_pop(LocalLua, 1);
            // clear the string buffer for the next iteration
            memset(magbuff, 0, sizeof(magbuff));
        }
    }
    lua_close(LocalLua);

    // Find flux density in each element;
    for(i=0; i<(int)meshelem.size(); i++) GetElementB(meshelem[i]);

    // Find extreme values of A;
    A_Low = meshnode[0].A.re;
    A_High = meshnode[0].A.re;
    for(i=1; i<(int)meshnode.size(); i++)
    {
        if (meshnode[i].A.re>A_High) A_High=meshnode[i].A.re;
        if (meshnode[i].A.re<A_Low)  A_Low =meshnode[i].A.re;

        if(Frequency!=0)
        {
            if (meshnode[i].A.im<A_Low)  A_Low =meshnode[i].A.im;
            if (meshnode[i].A.im>A_High) A_High=meshnode[i].A.im;
        }
    }
    // save default values for extremes of A
    A_lb=A_Low;
    A_ub=A_High;

    if(Frequency!=0)  // compute frequency-dependent permeabilities for linear blocks;
    {

        CComplex deg45;
        deg45=1+I;
        CComplex K,halflag;
        double ds;
        double w=2.*PI*Frequency;

        for(k=0; k<(int)blockproplist.size(); k++)
        {
            if (blockproplist[k].LamType==0)
            {
                blockproplist[k].mu_fdx = blockproplist[k].mu_x*
                                          exp(-I*blockproplist[k].Theta_hx*PI/180.);
                blockproplist[k].mu_fdy = blockproplist[k].mu_y*
                                          exp(-I*blockproplist[k].Theta_hy*PI/180.);

                if(blockproplist[k].Lam_d!=0)
                {

                    halflag = exp(-I*blockproplist[k].Theta_hx*PI/360.);

                    ds = sqrt(2. / (0.4 * PI * w * blockproplist[k].Cduct * blockproplist[k].mu_x));

                    K = halflag*deg45*blockproplist[k].Lam_d*0.001/(2.*ds);

                    if (blockproplist[k].Cduct!=0)
                    {
                        blockproplist[k].mu_fdx=(blockproplist[k].mu_fdx*tanh(K)/K)*
                                                blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
                    }
                    else
                    {
                        blockproplist[k].mu_fdx=(blockproplist[k].mu_fdx)*
                                                blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
                    }

                    halflag=exp(-I*blockproplist[k].Theta_hy*PI/360.);
                    ds=sqrt(2./(0.4*PI*w*blockproplist[k].Cduct*blockproplist[k].mu_y));
                    K=halflag*deg45*blockproplist[k].Lam_d*0.001/(2.*ds);
                    if (blockproplist[k].Cduct!=0)
                    {
                        blockproplist[k].mu_fdy=(blockproplist[k].mu_fdy*tanh(K)/K)*
                                                blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
                    }
                    else
                    {
                        blockproplist[k].mu_fdy=(blockproplist[k].mu_fdy)*
                                                blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
                    }
                }
            }

        }
    }

    // compute fill factor associated with each block label
    for(k=0; k<(int)blocklist.size(); k++)
    {
        GetFillFactor(k);
    }

    // build list of elements connected to each node;
    // allocate connections list;
    NumList=(int *)calloc(meshnode.size(),sizeof(int));
    ConList=(int **)calloc(meshnode.size(),sizeof(int *));
    // find out number of connections to each node;
    for(i=0; i<(int)meshelem.size(); i++)
        for(j=0; j<3; j++)
            NumList[meshelem[i].p[j]]++;
    // allocate space for connections lists;
    for(i=0; i<(int)meshnode.size(); i++)
        ConList[i]=(int *)calloc(NumList[i],sizeof(int));
    // build list;
    for(i=0; i<(int)meshnode.size(); i++) NumList[i]=0;
    for(i=0; i<(int)meshelem.size(); i++)
        for(j=0; j<3; j++)
        {
            k=meshelem[i].p[j];
            ConList[k][NumList[k]]=i;
            NumList[k]++;
        }

    // find extreme values of J;
    {
        CComplex Jelm[3],Aelm[3];

        double J_Low, J_High;
        double Jr_Low, Jr_High;
        double Ji_Low, Ji_High;

        GetJA(0,Jelm,Aelm);
        Jr_Low=fabs(Jelm[0].re);
        Jr_High=Jr_Low;
        Ji_Low=fabs(Jelm[0].im);
        Ji_High=Ji_Low;
        J_Low=abs(Jelm[0]);
        J_High=J_Low;
        for(i=0; i<(int)meshelem.size(); i++)
        {
            GetJA(i,Jelm,Aelm);
            for(j=0; j<3; j++)
            {
                br=fabs(Jelm[j].re);
                bi=fabs(Jelm[j].im);
                b=abs(Jelm[j]);

                if(b>J_High) J_High=b;
                if(b<J_Low) J_Low=b;
                if(br>Jr_High) Jr_High=br;
                if(br<Jr_Low) Jr_Low=br;
                if(bi>Ji_High) Ji_High=bi;
                if(bi<Ji_Low) Ji_Low=bi;
            }
        }

        J_Low*=1.e-6;
        J_High*=1e-6;
        Jr_Low*=1.e-6;
        Jr_High*=1e-6;
        Ji_Low*=1.e-6;
        Ji_High*=1e-6;

        if (Frequency==0)
        {
            d_PlotBounds[2][0]=PlotBounds[2][0]=J_Low;
            d_PlotBounds[2][1]=PlotBounds[2][1]=J_High;
        }
        else
        {
            d_PlotBounds[6][0]=PlotBounds[6][0]=J_Low;
            d_PlotBounds[6][1]=PlotBounds[6][1]=J_High;
            d_PlotBounds[7][0]=PlotBounds[7][0]=Jr_Low;
            d_PlotBounds[7][1]=PlotBounds[7][1]=Jr_High;
            d_PlotBounds[8][0]=PlotBounds[8][0]=Ji_Low;
            d_PlotBounds[8][1]=PlotBounds[8][1]=Ji_High;
        }
    }

    // Find extreme values of B and H;
    {
        double Br_Low, Br_High;
        double Bi_Low, Bi_High;
        double H_Low;
        double Hr_Low, Hr_High;
        double Hi_Low, Hi_High;
        CComplex h1,h2;

        Br_Low  = sqrt(sqr(meshelem[0].B1.re) + sqr(meshelem[0].B2.re));
        Br_High = Br_Low;
        Bi_Low  = sqrt(sqr(meshelem[0].B1.im)+ sqr(meshelem[0].B2.im));
        Bi_High = Bi_Low;
        B_Low   = sqrt(Br_Low*Br_Low + Bi_Low*Bi_Low);
        B_High  = B_Low;

        if (Frequency!=0)
            GetH(meshelem[0].B1,meshelem[0].B2,h1,h2,0);
        else
        {
            h1 = 0;
            h2 = 0;
            GetH(meshelem[0].B1.re,meshelem[0].B2.re,h1.re,h2.re,0);
        }

        Hr_Low  = sqrt(sqr(h1.re) + sqr(h2.re));
        Hr_High = Hr_Low;
        Hi_Low  = sqrt(sqr(h1.im)+ sqr(h2.im));
        Hi_High = Hi_Low;
        H_Low   = sqrt(Hr_Low*Hr_Low + Hi_Low*Hi_Low);
        H_High  = H_Low;

        for(i=0; i<(int)meshelem.size(); i++)
        {
            GetNodalB(meshelem[i].b1,meshelem[i].b2,meshelem[i]);
            for(j=0; j<3; j++)
            {
                br=sqrt(sqr(meshelem[i].b1[j].re) +
                        sqr(meshelem[i].b2[j].re));
                bi=sqrt(sqr(meshelem[i].b1[j].im) +
                        sqr(meshelem[i].b2[j].im));
                b=sqrt(br*br+bi*bi);

                if(b>B_High)   B_High=b;
                if(b<B_Low)   B_Low=b;
                if(br>Br_High) Br_High=br;
                if(br<Br_Low) Br_Low=br;
                if(bi>Bi_High) Bi_High=bi;
                if(bi<Bi_Low) Bi_Low=bi;
            }

            // getting lazy--just consider element averages for H
            if (Frequency!=0)
                GetH(meshelem[i].B1,meshelem[i].B2,h1,h2,i);
            else
                GetH(meshelem[i].B1.re,meshelem[i].B2.re,h1.re,h2.re,i);

            br=sqrt(sqr(h1.re) + sqr(h2.re));
            bi=sqrt(sqr(h1.im) + sqr(h2.im));
            b=sqrt(br*br+bi*bi);
            if(b>H_High)   H_High=b;
            if(b<H_Low)   H_Low=b;
            if(br>Hr_High) Hr_High=br;
            if(br<Hr_Low) Hr_Low=br;
            if(bi>Hi_High) Hi_High=bi;
            if(bi<Hi_Low) Hi_Low=bi;

        }

        if (Frequency==0)
        {
            d_PlotBounds[0][0]=PlotBounds[0][0]=B_Low;
            d_PlotBounds[0][1]=PlotBounds[0][1]=B_High;
            d_PlotBounds[1][0]=PlotBounds[1][0]=H_Low;
            d_PlotBounds[1][1]=PlotBounds[1][1]=H_High;
        }
        else
        {
            d_PlotBounds[0][0]=PlotBounds[0][0]=B_Low;
            d_PlotBounds[0][1]=PlotBounds[0][1]=B_High;
            d_PlotBounds[1][0]=PlotBounds[1][0]=Br_Low;
            d_PlotBounds[1][1]=PlotBounds[1][1]=Br_High;
            d_PlotBounds[2][0]=PlotBounds[2][0]=Bi_Low;
            d_PlotBounds[2][1]=PlotBounds[2][1]=Bi_High;
            d_PlotBounds[3][0]=PlotBounds[3][0]=H_Low;
            d_PlotBounds[3][1]=PlotBounds[3][1]=H_High;
            d_PlotBounds[4][0]=PlotBounds[4][0]=Hr_Low;
            d_PlotBounds[4][1]=PlotBounds[4][1]=Hr_High;
            d_PlotBounds[5][0]=PlotBounds[5][0]=Hi_Low;
            d_PlotBounds[5][1]=PlotBounds[5][1]=Hi_High;
        }
    }

//    // Choose bounds based on the type of contour plot
//    // currently in play
//    POSITION pos = GetFirstViewPosition();
//    CFemmviewView *theView=(CFemmviewView *)GetNextView(pos);
//
//    if(Frequency==0)
//    {
//        if (theView->DensityPlot==2) theView->DensityPlot=1;
//        if (theView->DensityPlot>1)  theView->DensityPlot=0;
//    }

    // compute total resulting current for circuits with an a priori defined
    // voltage gradient;  Need this to display circuit results & impedance.
    for(i=0; i<(int)circproplist.size(); i++)
    {
        CComplex Jelm[3],Aelm[3];
        double a;

        if(circproplist[i].CircType>1)
            for(j=0,circproplist[i].Amps=0.; j<(int)meshelem.size(); j++)
            {
                if(blocklist[meshelem[j].lbl].InCircuit==i)
                {

                    GetJA(j,Jelm,Aelm);
                    // Convert area units to metres
                    a = ElmArea(j) * sqr(LengthConv[LengthUnits]);
                    // Add the current in the element (J * Elemnet Area) to the total
                    for(k=0; k<3; k++) circproplist[i].Amps += a * Jelm[k]/3;
                }
            }
    }

    // Build adjacency information for each element.
    FindBoundaryEdges();

    // Check to see if any regions are multiply defined
    // (i.e. tagged by more than one block label). If so,
    // display an error message and mark the problem blocks.
    for(k=0,bMultiplyDefinedLabels=false; k<(int)blocklist.size(); k++)
    {
        // test if the label is inside the meshed region, by attempting to find
        // which triangle it is in, if it's outside the problem region it will
        // be ignored anyway
        if( (i = InTriangle(blocklist[k].x,blocklist[k].y)) >= 0 )
        {
            // the label is in the problem domain, test if the label assigned
            // to the element which the label is in has the same value as the
            // label number
            if(meshelem[i].lbl != k)
            {
                // if the label number assigned to the element is not the same as
                // the block label numer, there must be multiply defined labels for
                // the region

                // select the offending region
                blocklist[meshelem[i].lbl].IsSelected=TRUE;

                // if it the first multiply defined label we have found, issue a warning
                // and set the appropriate flag to true
                if (!bMultiplyDefinedLabels)
                {
                    string msg = "Some regions in the problem have been defined\n";
                    msg += "by more than one block label.  These potentially\n";
                    msg += "problematic regions will appear as selected in\n";
                    msg += "the initial view.";
                    WarnMessage(msg.c_str());
                    bMultiplyDefinedLabels = true;
                }
            }
        }
    }


    // Get some information needed to compute energy stored in
    // permanent magnets with a nonlinear demagnetization curve
    if (Frequency==0)
    {
        for(k=0; k<(int)blockproplist.size(); k++)
        {
            if ((blockproplist[k].H_c>0) && (blockproplist[k].BHpoints>0))
            {
                blockproplist[k].Nrg = blockproplist[k].GetCoEnergy(blockproplist[k].GetB(blockproplist[k].H_c));
            }
        }
    }

    FirstDraw=true;
    return true;
}


int FPProc::InTriangle(double x, double y)
{
    static int k;
    int j,hi,lo,sz;
    double z;

    sz=meshelem.size();
    if((k<0) || (k>=sz)) k=0;

    // In most applications, the triangle we're looking
    // for is nearby the last one we found.  Since the elements
    // are ordered in a banded structure, we want to check the
    // elements nearby the last one selected first.
    if (InTriangleTest(x,y,k)) return k;

    hi=k;
    lo=k;

    for(j=0; j<sz; j+=2)
    {
        hi++;
        if(hi>=sz) hi=0;
        lo--;
        if(lo<0)   lo=sz-1;

        z=(meshelem[hi].ctr.re-x)*(meshelem[hi].ctr.re-x) +
          (meshelem[hi].ctr.im-y)*(meshelem[hi].ctr.im-y);
        if(z<=meshelem[hi].rsqr)
        {
            if(InTriangleTest(x,y,hi))
            {
                k=hi;
                return k;
            }
        }

        z=(meshelem[lo].ctr.re-x)*(meshelem[lo].ctr.re-x) +
          (meshelem[lo].ctr.im-y)*(meshelem[lo].ctr.im-y);
        if(z<=meshelem[lo].rsqr)
        {
            if(InTriangleTest(x,y,lo))
            {
                k=lo;
                return k;
            }
        }

    }

    return (-1);
}

bool FPProc::GetPointValues(double x, double y, CPointVals &u)
{
    int k;

    // find the mesh triangle in which x,y resides, if any
    k = InTriangle(x,y);

    if (k<0)
    {
        // we are not in a triangle so return false
        return false;
    }

    // we are in a triangle 'k', we can interpolate the point
    // values in the triangle
    GetPointValues(x,y,k,u);

    return true;
}

bool FPProc::GetPointValues(double x, double y, int k, CPointVals &u)
{
    int i,j,n[3],lbl;
    double a[3],b[3],c[3],da,ravg;

    for(i=0; i<3; i++)
    {
        // get the nodes of the mesh element 'k'
        n[i] = meshelem[k].p[i];
    }

    a[0] = meshnode[n[1]].x * meshnode[n[2]].y - meshnode[n[2]].x * meshnode[n[1]].y;
    a[1] = meshnode[n[2]].x * meshnode[n[0]].y - meshnode[n[0]].x * meshnode[n[2]].y;
    a[2] = meshnode[n[0]].x * meshnode[n[1]].y - meshnode[n[1]].x * meshnode[n[0]].y;
    b[0] = meshnode[n[1]].y - meshnode[n[2]].y;
    b[1] = meshnode[n[2]].y - meshnode[n[0]].y;
    b[2] = meshnode[n[0]].y - meshnode[n[1]].y;
    c[0] = meshnode[n[2]].x - meshnode[n[1]].x;
    c[1] = meshnode[n[0]].x - meshnode[n[2]].x;
    c[2] = meshnode[n[1]].x - meshnode[n[0]].x;

    da = ( b[0]*c[1] - b[1]*c[0] );

    ravg = LengthConv[LengthUnits]*
           (meshnode[n[0]].x + meshnode[n[1]].x + meshnode[n[2]].x)/3.;

    // interpolate the flux density B at the given point in the element
    GetPointB(x,y,u.B1,u.B2,meshelem[k]);

    u.Hc=0;
//    if(blockproplist[meshelem[k].blk].LamType>2)
    u.ff=blocklist[meshelem[k].lbl].FillFactor;
//    else u.ff=-1;

    if (Frequency==0)
    {
        u.A = 0;
        if(ProblemType==PLANAR)
        {
            for(i=0; i<3; i++)
                u.A.re += meshnode[n[i]].A.re * (a[i] + b[i] * x + c[i] * y) / (da);
        }
        else
        {
            /* Old way that I interpolated potential in axi case:
                // interpolation from A from nodal points.
                // note that the potential that's actually stored
                // for axisymmetric problems is 2*Pi*r*A, so divide
                // by nodal r to get 2*Pi*A at the nodes.  Linearly
                // interpolate this, then multiply by r at the point
                // of interest to get back to 2*Pi*r*A.
                for(i=0,rp=0;i<3;i++){
                    r=meshnode[n[i]].x;
                    rp+=meshnode[n[i]].x*(a[i]+b[i]*x+c[i]*y)/da;
                    if (r>1.e-6) u.A.re+=meshnode[n[i]].A.re*
                        (a[i]+b[i]*x+c[i]*y)/(r*da);
                }
                u.A.re*=rp;
            */


            // a ``smarter'' interpolation.  One based on A can't
            // represent constant flux density very well.
            // This works, but I should re-write it in a more
            // efficient form--it's doing a lot of the work
            // twice, because the a-b-c stuff that's already
            // been computed is ignored.
            double v[6];
            double R[3];
//            double Z[3];
            double p,q;

            for(i=0; i<3; i++)
            {
                R[i]=meshnode[n[i]].x;
//                Z[i]=meshnode[n[i]].y;
            }

            // corner nodes
            v[0]=meshnode[n[0]].A.re;
            v[2]=meshnode[n[1]].A.re;
            v[4]=meshnode[n[2]].A.re;

            // construct values for mid-side nodes;
            if ((R[0]<1.e-06) && (R[1]<1.e-06))
                v[1]=(v[0]+v[2])/2.;
            else
                v[1]=(R[1]*(3.*v[0] + v[2]) + R[0]*(v[0] + 3.*v[2]))/
                     (4.*(R[0] + R[1]));

            if ((R[1]<1.e-06) && (R[2]<1.e-06))
                v[3]=(v[2]+v[4])/2.;
            else
                v[3]=(R[2]*(3.*v[2] + v[4]) + R[1]*(v[2] + 3.*v[4]))/
                     (4.*(R[1] + R[2]));

            if ((R[2]<1.e-06) && (R[0]<1.e-06))
                v[5]=(v[4]+v[0])/2.;
            else
                v[5]=(R[0]*(3.*v[4] + v[0]) + R[2]*(v[4] + 3.*v[0]))/
                     (4.*(R[2] + R[0]));

            // compute location in element transformed onto
            // a unit triangle;
            p=(b[1]*x+c[1]*y + a[1])/da;
            q=(b[2]*x+c[2]*y + a[2])/da;

            // now, interpolate to get potential...
            u.A.re = v[0] - p*(3.*v[0] - 4.*v[1] + v[2]) +
                     2.*p*p*(v[0] - 2.*v[1] + v[2]) -
                     q*(3.*v[0] + v[4] - 4.*v[5]) +
                     2.*q*q*(v[0] + v[4] - 2.*v[5]) +
                     4.*p*q*(v[0] - v[1] + v[3] - v[5]);

            /*        // "simple" way to do it...
                    // problem is that this mucks up things
                    // near the centerline, where things ought
                    // to look pretty quadratic.
                    for(i=0;i<3;i++)
                        u.A.re+=meshnode[n[i]].A.re*(a[i]+b[i]*x+c[i]*y)/(da);
            */
        }
        u.mu1.im = 0;
        u.mu2.im = 0;
        GetMu(u.B1.re,u.B2.re,u.mu1.re,u.mu2.re,k);
        u.H1 = u.B1/(Re(u.mu1)*muo);
        u.H2 = u.B2/(Re(u.mu2)*muo);
        u.Je = 0;
        u.Js = blockproplist[meshelem[k].blk].Jr;
        lbl  = meshelem[k].lbl;
        j = blocklist[lbl].InCircuit;
        if(j>=0)
        {
            if(blocklist[lbl].Case==0)
            {
                if (ProblemType==PLANAR)
                    u.Js-=Re(blocklist[meshelem[k].lbl].o)*
                          blocklist[lbl].dVolts;
                else
                {

                    int tn;
                    double R[3];
                    for(tn=0; tn<3; tn++)
                    {
                        R[tn] = meshnode[n[tn]].x;
                        if (R[tn]<1.e-6) R[tn] = ravg;
                        else R[tn] *= LengthConv[LengthUnits];
                    }
                    for(ravg=0.,tn=0; tn<3; tn++)
                        ravg+=(1./R[tn])*(a[tn]+b[tn]*x+c[tn]*y)/(da);
                    u.Js-=Re(blocklist[meshelem[k].lbl].o)*
                          blocklist[lbl].dVolts*ravg;
                }
            }
            else u.Js+=blocklist[lbl].J;
        }
        u.c = Re(blocklist[meshelem[k].lbl].o);
        u.E = blockproplist[meshelem[k].blk].DoEnergy(u.B1.re,u.B2.re);

        // correct H and energy stored in magnet for second-quadrant
        // representation of a PM.
        if (blockproplist[meshelem[k].blk].H_c!=0)
        {
            int bk=meshelem[k].blk;

            u.Hc = blockproplist[bk].H_c*exp(I*PI*meshelem[k].magdir/180.);
            u.H1 = u.H1-Re(u.Hc);
            u.H2 = u.H2-Im(u.Hc);

            // in the linear case:
            if (blockproplist[bk].BHpoints==0)
                u.E = 0.5*muo*(u.mu1.re*u.H1.re*u.H1.re + u.mu2.re*u.H2.re*u.H2.re);
            else
            {
                u.E = u.E + blockproplist[bk].Nrg
                      - blockproplist[bk].H_c*Re((u.B1.re+I*u.B2.re)/exp(I*PI*meshelem[k].magdir/180.));
            }

            // If considering the magnet as an equivalent coil, add Hc to the demagnetizing field
            if (!d_ShiftH)
            {
                u.H1 = u.H1 + Re(u.Hc);
                u.H2 = u.H2 + Im(u.Hc);
                u.Hc = 0;
            }
        }

        // add in "local" stored energy for wound that would be subject to
        // prox and skin effect for nonzero frequency cases.
        if (blockproplist[meshelem[k].blk].LamType>2)
        {
            CComplex J;
            J=u.Js*1.e6;

            u.E+=Re(J*J)*Im(blocklist[meshelem[i].lbl].o)/2.;
        }

        u.Ph=0;
        u.Pe=0;
        return true;
    }

    if(Frequency!=0)
    {
        u.A=0;
        if(ProblemType==PLANAR)
        {
            for(i=0; i<3; i++)
                u.A+=meshnode[n[i]].A*(a[i]+b[i]*x+c[i]*y)/(da);
        }
        else
        {
            CComplex v[6];
            double R[3];
//            double Z[3];
            double p,q;

            for(i=0; i<3; i++)
            {
                R[i]=meshnode[n[i]].x;
//                Z[i]=meshnode[n[i]].y;
            }

            // corner nodes
            v[0]=meshnode[n[0]].A;
            v[2]=meshnode[n[1]].A;
            v[4]=meshnode[n[2]].A;

            // construct values for mid-side nodes;
            if ((R[0]<1.e-06) && (R[1]<1.e-06))
                v[1]=(v[0]+v[2])/2.;
            else
                v[1]=(R[1]*(3.*v[0] + v[2]) + R[0]*(v[0] + 3.*v[2]))/
                     (4.*(R[0] + R[1]));

            if ((R[1]<1.e-06) && (R[2]<1.e-06))
                v[3]=(v[2]+v[4])/2.;
            else
                v[3]=(R[2]*(3.*v[2] + v[4]) + R[1]*(v[2] + 3.*v[4]))/
                     (4.*(R[1] + R[2]));

            if ((R[2]<1.e-06) && (R[0]<1.e-06))
                v[5]=(v[4]+v[0])/2.;
            else
                v[5]=(R[0]*(3.*v[4] + v[0]) + R[2]*(v[4] + 3.*v[0]))/
                     (4.*(R[2] + R[0]));

            // compute location in element transformed onto
            // a unit triangle;
            p=(b[1]*x+c[1]*y + a[1])/da;
            q=(b[2]*x+c[2]*y + a[2])/da;

            // now, interpolate to get potential...
            u.A = v[0] - p*(3.*v[0] - 4.*v[1] + v[2]) +
                  2.*p*p*(v[0] - 2.*v[1] + v[2]) -
                  q*(3.*v[0] + v[4] - 4.*v[5]) +
                  2.*q*q*(v[0] + v[4] - 2.*v[5]) +
                  4.*p*q*(v[0] - v[1] + v[3] - v[5]);
        }

        GetMu(u.B1,u.B2,u.mu1,u.mu2,k);
        u.H1 = u.B1/(u.mu1*muo);
        u.H2 = u.B2/(u.mu2*muo);

        u.Js=blockproplist[meshelem[k].blk].Jr +
             I*blockproplist[meshelem[k].blk].Ji;
        lbl=meshelem[k].lbl;
        j=blocklist[lbl].InCircuit;
        if(j>=0)
        {
            if(blocklist[lbl].Case==0)
            {
                if (ProblemType==PLANAR)
                    u.Js-=blocklist[meshelem[k].lbl].o*blocklist[lbl].dVolts;
                else
                {

                    int tn;
                    double R[3];
                    for(tn=0; tn<3; tn++)
                    {
                        R[tn]=meshnode[n[tn]].x;
                        if (R[tn]<1.e-6) R[tn]=ravg;
                        else R[tn]*=LengthConv[LengthUnits];
                    }
                    for(ravg=0.,tn=0; tn<3; tn++)
                        ravg+=(1./R[tn])*(a[tn]+b[tn]*x+c[tn]*y)/(da);
                    u.Js-=blocklist[meshelem[k].lbl].o*
                          blocklist[lbl].dVolts*ravg;
                }
            }
            else u.Js+=blocklist[lbl].J;
        }

        // report just loss-related part of conductivity.
        if (blockproplist[meshelem[k].blk].Cduct!=0)
            u.c=1./Re(1./(blocklist[meshelem[k].lbl].o));
        else u.c=0;

        if (blockproplist[meshelem[k].blk].Lam_d!=0) u.c=0;

        // only add in eddy currents if the region is solid
        if (blocklist[meshelem[k].lbl].FillFactor<0)
            u.Je=-I*Frequency*2.*PI*u.c*u.A;

        if(ProblemType!=0)
        {
            if(x!=0)
                u.Je/=(2.*PI*x*LengthConv[LengthUnits]);
            else u.Je=0;
        }

        CComplex z;
        z=(u.H1*u.B1.Conj()) + (u.H2*u.B2.Conj());
        u.E=0.25*z.re;

        // add in "local" stored energy for wound that would be subject to
        // prox and skin effect for nonzero frequency cases.
        if (blockproplist[meshelem[k].blk].LamType>2)
        {
            CComplex J;
            J=u.Js*1.e6;

            u.E += Re(J*conj(J))*(Im(1./blocklist[meshelem[k].lbl].o)/(2.e6*PI*Frequency))/4.;
        }




        u.Ph=Frequency*PI*z.im;
        u.Pe=0;
        if (u.c!=0)
        {
            z=u.Js + u.Je;
            u.Pe=1.e06*(z.re*z.re + z.im*z.im)/(u.c*2.);
        }

        return true;
    }

    return false;
}

void FPProc::GetPointB(double x, double y, CComplex &B1, CComplex &B2,
                       CElement &elm)
{
    // elm is a reference to the element that contains the point of interest.
    int i,n[3];
    double da,a[3],b[3],c[3];

    if(Smooth==false)
    {
        B1=elm.B1;
        B2=elm.B2;
        return;
    }

    for(i=0; i<3; i++) n[i]=elm.p[i];
    a[0]=meshnode[n[1]].x * meshnode[n[2]].y - meshnode[n[2]].x * meshnode[n[1]].y;
    a[1]=meshnode[n[2]].x * meshnode[n[0]].y - meshnode[n[0]].x * meshnode[n[2]].y;
    a[2]=meshnode[n[0]].x * meshnode[n[1]].y - meshnode[n[1]].x * meshnode[n[0]].y;
    b[0]=meshnode[n[1]].y - meshnode[n[2]].y;
    b[1]=meshnode[n[2]].y - meshnode[n[0]].y;
    b[2]=meshnode[n[0]].y - meshnode[n[1]].y;
    c[0]=meshnode[n[2]].x - meshnode[n[1]].x;
    c[1]=meshnode[n[0]].x - meshnode[n[2]].x;
    c[2]=meshnode[n[1]].x - meshnode[n[0]].x;
    da=(b[0]*c[1]-b[1]*c[0]);

    B1.Set(0,0);
    B2.Set(0,0);
    for(i=0; i<3; i++)
    {
        B1+=(elm.b1[i]*(a[i]+b[i]*x+c[i]*y)/da);
        B2+=(elm.b2[i]*(a[i]+b[i]*x+c[i]*y)/da);
    }
}

void FPProc::GetNodalB(CComplex *b1, CComplex *b2,CElement &elm)
{
    // elm is a reference to the element that contains the point of interest.
    CComplex p;
    CComplex tn,bn,bt,btu,btv,u1,u2,v1,v2;
    int i,j,k,l,q,m,pt,nxt;
    i=j=k=l=q=m=pt=nxt = 0;
    double r,R,z;
    r=R=z = 0;
    CElement *e;
    int flag;

    // find nodal values of flux density via a patch method.
    for(i=0; i<3; i++)
    {

        k=elm.p[i];
        p.Set(meshnode[k].x,meshnode[k].y);
        b1[i].Set(0,0);
        b2[i].Set(0,0);
        for(j=0,m=0; j<NumList[k]; j++)
            if(elm.lbl==meshelem[ConList[k][j]].lbl) m++;
            else
            {
                if(Frequency==0)
                {
                    if ((blockproplist[elm.blk].mu_x==
                            blockproplist[meshelem[ConList[k][j]].blk].mu_x) &&
                            (blockproplist[elm.blk].mu_y==
                             blockproplist[meshelem[ConList[k][j]].blk].mu_y) &&
                            (blockproplist[elm.blk].H_c==
                             blockproplist[meshelem[ConList[k][j]].blk].H_c) &&
                            (elm.magdir==meshelem[ConList[k][j]].magdir)) m++;
                    else if ((elm.blk==meshelem[ConList[k][j]].blk) &&
                             (elm.magdir==meshelem[ConList[k][j]].magdir)) m++;
                }
                else if ((blockproplist[elm.blk].mu_fdx==
                          blockproplist[meshelem[ConList[k][j]].blk].mu_fdx) &&
                         (blockproplist[elm.blk].mu_fdy==
                          blockproplist[meshelem[ConList[k][j]].blk].mu_fdy)) m++;
            }

        if(m==NumList[k]) // normal smoothing method for points
        {
            // away from any boundaries
            for(j=0,R=0; j<NumList[k]; j++)
            {
                m=ConList[k][j];
                z=1./abs(p-Ctr(m));
                R+=z;
                b1[i]+=(z*meshelem[m].B1);
                b2[i]+=(z*meshelem[m].B2);
            }
            b1[i]/=R;
            b2[i]/=R;
        }

        else
        {
            R=0;
            v1=0;
            v2=0;

            //scan ccw for an interface...
            e=&elm;
            for(q=0; q<NumList[k]; q++)
            {
                //find ccw side of the element;
                for(j=0; j<3; j++) if(e->p[j]==k) pt=j;
                pt--;
                if(pt<0) pt=2;
                pt=e->p[pt];

                //scan to find element adjacent to this side;
                for(j=0,nxt=-1; j<NumList[k]; j++)
                {
                    if(&meshelem[ConList[k][j]]!=e)
                    {
                        for(l=0; l<3; l++)
                            if(meshelem[ConList[k][j]].p[l]==pt)
                                nxt=ConList[k][j];
                    }
                }

                if(nxt==-1)
                {
                    // a special-case punt
                    q=NumList[k];
                    b1[i]=(e->B1);
                    b2[i]=(e->B2);
                    v1=1;
                    v2=1;
                }
                else if(elm.lbl!=meshelem[nxt].lbl)
                {
                    // we have found two elements on either side of the interface
                    // now, we take contribution from B at the center of the
                    // interface side
                    tn.Set(meshnode[pt].x-meshnode[k].x,
                           meshnode[pt].y-meshnode[k].y);
                    r=(meshnode[pt].x+meshnode[k].x)*LengthConv[LengthUnits]/2.;
                    bn=(meshnode[pt].A-meshnode[k].A)/
                       (abs(tn)*LengthConv[LengthUnits]);
                    if(ProblemType==AXISYMMETRIC)
                    {
                        bn/=(-2.*PI*r);
                    }
                    z=0.5/abs(tn);
                    tn/=abs(tn);

                    // for the moment, kludge with bt...
                    bt=e->B1*tn.re + e->B2*tn.im;

                    R+=z;
                    b1[i]+=(z*tn.re*bt);
                    b2[i]+=(z*tn.im*bt);
                    b1[i]+=(z*tn.im*bn);
                    b2[i]+=(-z*tn.re*bn);
                    v1=tn;
                    q=NumList[k];
                }
                else e=&meshelem[nxt];
            }

            //scan cw for an interface...
            if(v2==0) // catches the "special-case punt" where we have
            {
                // already set nodal B values....
                e=&elm;
                for(q=0; q<NumList[k]; q++)
                {
                    //find cw side of the element;
                    for(j=0; j<3; j++) if(e->p[j]==k) pt=j;
                    pt++;
                    if(pt>2) pt=0;
                    pt=e->p[pt];

                    //scan to find element adjacent to this side;
                    for(j=0,nxt=-1; j<NumList[k]; j++)
                    {
                        if(&meshelem[ConList[k][j]]!=e)
                        {
                            for(l=0; l<3; l++)
                                if(meshelem[ConList[k][j]].p[l]==pt)
                                    nxt=ConList[k][j];
                        }
                    }
                    if (nxt==-1)
                    {
                        // a special-case punt
                        q=NumList[k];
                        b1[i]=(e->B1);
                        b2[i]=(e->B2);
                        v1=1;
                        v2=1;
                    }
                    else if(elm.lbl!=meshelem[nxt].lbl)
                    {
                        // we have found two elements on either side of the interface
                        // now, we take contribution from B at the center of the
                        // interface side
                        tn.Set(meshnode[pt].x-meshnode[k].x,
                               meshnode[pt].y-meshnode[k].y);
                        r=(meshnode[pt].x+meshnode[k].x)*LengthConv[LengthUnits]/2.;
                        bn=(meshnode[pt].A-meshnode[k].A)/
                           (abs(tn)*LengthConv[LengthUnits]);
                        if(ProblemType==AXISYMMETRIC)
                        {
                            bn/=(-2.*PI*r);
                        }
                        z=0.5/abs(tn);
                        tn/=abs(tn);

                        // for the moment, kludge with bt...
                        bt=e->B1*tn.re + e->B2*tn.im;

                        R+=z;
                        b1[i]+=(z*tn.re*bt);
                        b2[i]+=(z*tn.im*bt);
                        b1[i]+=(z*tn.im*bn);
                        b2[i]+=(-z*tn.re*bn);
                        v2=tn;
                        q=NumList[k];
                    }
                    else e=&meshelem[nxt];
                }
                b1[i]/=R;
                b2[i]/=R;
            }

            // check to see if angle of corner is too sharp to apply
            // this rule; really only does right if the interface is flat;
            flag=false;
            // if there is only one edge, approx is ok;
            if ((abs(v1)<0.9) || (abs(v2)<0.9)) flag=true;
            // if the interfaces make less than a 10 degree angle, things are ok;
            if ( (-v1.re*v2.re-v1.im*v2.im) > 0.985) flag=true;

            // Otherwise, punt...
            if(flag==false)
            {
                bn=0;
                k=elm.p[i];
                for(j=0; j<NumList[k]; j++)
                {
                    if(elm.lbl==meshelem[ConList[k][j]].lbl)
                    {
                        m=ConList[k][j];
                        bt.re=sqrt(meshelem[m].B1.re*meshelem[m].B1.re +
                                   meshelem[m].B2.re*meshelem[m].B2.re);
                        bt.im=sqrt(meshelem[m].B1.im*meshelem[m].B1.im +
                                   meshelem[m].B2.im*meshelem[m].B2.im);
                        if(bt.re>bn.re) bn.re=bt.re;
                        if(bt.im>bn.im) bn.im=bt.im;
                    }
                }

                R=sqrt(elm.B1.re*elm.B1.re + elm.B2.re*elm.B2.re);
                if(R!=0)
                {
                    b1[i].re=bn.re/R * elm.B1.re;
                    b2[i].re=bn.re/R * elm.B2.re;
                }
                else
                {
                    b1[i].re=0;
                    b2[i].re=0;
                }

                R=sqrt(elm.B1.im*elm.B1.im + elm.B2.im*elm.B2.im);
                if(R!=0)
                {
                    b1[i].im=bn.im/R * elm.B1.im;
                    b2[i].im=bn.im/R * elm.B2.im;
                }
                else
                {
                    b1[i].im=0;
                    b2[i].im=0;
                }
            }
        }

        // check to see if the point has a point current; if so, just
        // use element average values;
        if (nodeproplist.size()!=0)
            for(j=0; j<(int)nodelist.size(); j++)
            {
                if (abs(p-(nodelist[j].x+nodelist[j].y*I))<1.e-08)
                    if(nodelist[j].BoundaryMarker>=0)
                    {
                        if ((nodeproplist[nodelist[j].BoundaryMarker].Jr!=0) ||
                                (nodeproplist[nodelist[j].BoundaryMarker].Ji!=0))
                        {
                            b1[i]=elm.B1;
                            b2[i]=elm.B2;
                        }
                    }
            }


        //check for special case of node on r=0 axisymmetric; set Br=0;
        if ((fabs(p.re)<1.e-06) && (ProblemType==AXISYMMETRIC)) b1[i].Set(0.,0);
    }
}

void FPProc::GetElementB(CElement &elm)
{
    int i,n[3];
    double b[3],c[3],da;

    for(i=0; i<3; i++) n[i]=elm.p[i];

    b[0]=meshnode[n[1]].y - meshnode[n[2]].y;
    b[1]=meshnode[n[2]].y - meshnode[n[0]].y;
    b[2]=meshnode[n[0]].y - meshnode[n[1]].y;
    c[0]=meshnode[n[2]].x - meshnode[n[1]].x;
    c[1]=meshnode[n[0]].x - meshnode[n[2]].x;
    c[2]=meshnode[n[1]].x - meshnode[n[0]].x;
    da=(b[0]*c[1]-b[1]*c[0]);

    if (ProblemType==PLANAR)
    {
        elm.B1=0;
        elm.B2=0;
        for(i=0; i<3; i++)
        {
            elm.B1+=meshnode[n[i]].A*c[i]/(da*LengthConv[LengthUnits]);
            elm.B2-=meshnode[n[i]].A*b[i]/(da*LengthConv[LengthUnits]);
        }
        return;
    }
    else
    {
        CComplex v[6],dp,dq;
        double R[3],r; //Z[3]

        for(i=0,r=0; i<3; i++)
        {
            R[i]=meshnode[n[i]].x;
//            Z[i]=meshnode[n[i]].y;
            r+=R[i]/3.;
        }

        // corner nodes
        v[0]=meshnode[n[0]].A;
        v[2]=meshnode[n[1]].A;
        v[4]=meshnode[n[2]].A;

        // construct values for mid-side nodes;
        if ((R[0]<1.e-06) && (R[1]<1.e-06)) v[1]=(v[0]+v[2])/2.;
        else v[1]=(R[1]*(3.*v[0] + v[2]) + R[0]*(v[0] + 3.*v[2]))/
                      (4.*(R[0] + R[1]));

        if ((R[1]<1.e-06) && (R[2]<1.e-06)) v[3]=(v[2]+v[4])/2.;
        else v[3]=(R[2]*(3.*v[2] + v[4]) + R[1]*(v[2] + 3.*v[4]))/
                      (4.*(R[1] + R[2]));

        if ((R[2]<1.e-06) && (R[0]<1.e-06)) v[5]=(v[4]+v[0])/2.;
        else v[5]=(R[0]*(3.*v[4] + v[0]) + R[2]*(v[4] + 3.*v[0]))/
                      (4.*(R[2] + R[0]));

        // derivatives w.r.t. p and q:
        dp=(-v[0] + v[2] + 4.*v[3] - 4.*v[5])/3.;
        dq=(-v[0] - 4.*v[1] + 4.*v[3] + v[4])/3.;

        // now, compute flux.
        da*=2.*PI*r*LengthConv[LengthUnits]*LengthConv[LengthUnits];
        elm.B1=-(c[1]*dp+c[2]*dq)/da;
        elm.B2= (b[1]*dp+b[2]*dq)/da;

        return;
    }
}


// void FPProc::OnReload()
// {
//     // TODO: Add your command handler code here
//     CString pname = GetPathName();
//     if(pname.GetLength()>0)
//     {
//         NewDocument();
//         SetPathName(pname,false);
//         OpenDocument(pname);
//     }
// }

int FPProc::ClosestNode(double x, double y)
{
    int i,j;
    double d0,d1;

    if(nodelist.size()==0) return -1;

    j=0;
    d0=nodelist[0].GetDistance(x,y);
    for(i=0; i<(int)nodelist.size(); i++)
    {
        d1=nodelist[i].GetDistance(x,y);
        if(d1<d0)
        {
            d0=d1;
            j=i;
        }
    }

    return j;
}

//void FPProc::GetLineValues(CXYPlot &p,int PlotType,int NumPlotPoints)
//{
//    double *q,z,u,dz;
//    CComplex pt,n,t;
//    int i,j,k,m,elm;
//    CPointVals v;
//    bool flag;
//
//    q=(double *)calloc(contour.size(),sizeof(double));
//    for(i=1,z=0.; i<contour.size(); i++)
//    {
//        z+=abs(contour[i]-contour[i-1]);
//        q[i]=z;
//    }
//    dz=z/(NumPlotPoints-1);
//
//    /*
//        m_XYPlotType.AddString("Potential");
//        m_XYPlotType.AddString("|B|        (Magnitude of flux density)");
//        m_XYPlotType.AddString("B . n      (Normal flux density)");
//        m_XYPlotType.AddString("B . t      (Tangential flux density)");
//        m_XYPlotType.AddString("|H|        (Magnitude of field intensity)");
//        m_XYPlotType.AddString("H . n      (Normal field intensity)");
//        m_XYPlotType.AddString("H . t      (Tangential field intensity)");
//        m_XYPlotType.AddString("J_eddy
//    */
//
//    if(Frequency==0)
//    {
//        switch (PlotType)
//        {
//        case 0:
//            p.Create(NumPlotPoints,2);
//            if (ProblemType==PLANAR) strcpy(p.lbls[1],"Potential, Wb/m");
//            else strcpy(p.lbls[1],"Flux, Wb");
//            break;
//        case 1:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"|B|, Tesla");
//            break;
//        case 2:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"B.n, Tesla");
//            break;
//        case 3:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"B.t, Tesla");
//            break;
//        case 4:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"|H|, Amp/m");
//            break;
//        case 5:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"H.n, Amp/m");
//            break;
//        case 6:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"H.t, Amp/m");
//            break;
//        default:
//            p.Create(NumPlotPoints,2);
//            break;
//        }
//    }
//    else
//    {
//        switch (PlotType)
//        {
//        case 0:
//            p.Create(NumPlotPoints,4);
//            if(ProblemType==PLANAR)
//            {
//                strcpy(p.lbls[1],"|A|, Wb/m");
//                strcpy(p.lbls[2],"Re[A], Wb/m");
//                strcpy(p.lbls[3],"Im[A], Wb/m");
//            }
//            else
//            {
//                strcpy(p.lbls[1],"|Flux|, Wb");
//                strcpy(p.lbls[2],"Re[Flux], Wb");
//                strcpy(p.lbls[3],"Im[Flux], Wb");
//            }
//            break;
//        case 1:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"|B|, Tesla");
//            break;
//        case 2:
//            p.Create(NumPlotPoints,4);
//            strcpy(p.lbls[1],"|B.n|, Tesla");
//            strcpy(p.lbls[2],"Re[B.n], Tesla");
//            strcpy(p.lbls[3],"Im[B.n], Tesla");
//            break;
//        case 3:
//            p.Create(NumPlotPoints,4);
//            strcpy(p.lbls[1],"|B.t|, Tesla");
//            strcpy(p.lbls[2],"Re[B.t], Tesla");
//            strcpy(p.lbls[3],"Im[B.t], Tesla");
//            break;
//        case 4:
//            p.Create(NumPlotPoints,2);
//            strcpy(p.lbls[1],"|H|, Amp/m");
//            break;
//        case 5:
//            p.Create(NumPlotPoints,4);
//            if(ProblemType==PLANAR)
//            {
//                strcpy(p.lbls[1],"|H.n|, Amp/m");
//                strcpy(p.lbls[2],"Re[H.n], Amp/m");
//                strcpy(p.lbls[3],"Im[H.n], Amp/m");
//            }
//            break;
//        case 6:
//            p.Create(NumPlotPoints,4);
//            strcpy(p.lbls[1],"|H.t|, Amp/m");
//            strcpy(p.lbls[2],"Re[H.t], Amp/m");
//            strcpy(p.lbls[3],"Im[H.t], Amp/m");
//            break;
//        case 7:
//            p.Create(NumPlotPoints,4);
//            strcpy(p.lbls[1],"|Je|, MA/m^2");
//            strcpy(p.lbls[2],"Re[Je], MA/m^2");
//            strcpy(p.lbls[3],"Im[Je], MA/m^2");
//            break;
//        case 8:
//            p.Create(NumPlotPoints,4);
//            strcpy(p.lbls[1],"|Js+Je|, MA/m^2");
//            strcpy(p.lbls[2],"Re[Js+Je], MA/m^2");
//            strcpy(p.lbls[3],"Im[Js+Je], MA/m^2");
//            break;
//        default:
//            p.Create(NumPlotPoints,2);
//            break;
//        }
//    }
//
//    switch(LengthUnits)
//    {
//    case 1:
//        strcpy(p.lbls[0],"Length, mm");
//        break;
//    case 2:
//        strcpy(p.lbls[0],"Length, cm");
//        break;
//    case 3:
//        strcpy(p.lbls[0],"Length, m");
//        break;
//    case 4:
//        strcpy(p.lbls[0],"Length, mils");
//        break;
//    case 5:
//        strcpy(p.lbls[0],"Length, um");
//        break;
//    default:
//        strcpy(p.lbls[0],"Length, inches");
//        break;
//    }
//
//    for(i=0,k=1,z=0,elm=-1; i<NumPlotPoints; i++,z+=dz)
//    {
//        while((z>q[k]) && (k<(contour.size()-1))) k++;
//        u=(z-q[k-1])/(q[k]-q[k-1]);
//        pt=contour[k-1]+u*(contour[k]-contour[k-1]);
//        t=contour[k]-contour[k-1];
//        t/=abs(t);
//        n = I*t;
//        pt+=(n*1.e-06);
//
//        if (elm<0) elm=InTriangle(pt.re,pt.im);
//        else if (InTriangleTest(pt.re,pt.im,elm)==false)
//        {
//            flag=false;
//            for(j=0; j<3; j++)
//                for(m=0; m<NumList[meshelem[elm].p[j]]; m++)
//                {
//                    elm=ConList[meshelem[elm].p[j]][m];
//                    if (InTriangleTest(pt.re,pt.im,elm)==true)
//                    {
//                        flag=true;
//                        m=100;
//                        j=3;
//                    }
//                }
//            if (flag==false) elm=InTriangle(pt.re,pt.im);
//        }
//        if(elm>=0)
//            flag=GetPointValues(pt.re,pt.im,elm,v);
//        else flag=false;
//
//        p.M[i][0]=z;
//        if ((Frequency==0) && (flag!=false))
//        {
//            switch (PlotType)
//            {
//            case 0:
//                p.M[i][1]=v.A.re;
//                break;
//            case 1:
//                p.M[i][1]=sqrt(v.B1.Abs()*v.B1.Abs() + v.B2.Abs()*v.B2.Abs());
//                break;
//            case 2:
//                p.M[i][1]= n.re*v.B1.re + n.im*v.B2.re;
//                break;
//            case 3:
//                p.M[i][1]= t.re*v.B1.re + t.im*v.B2.re;
//                break;
//            case 4:
//                p.M[i][1]=sqrt(v.H1.Abs()*v.H1.Abs() + v.H2.Abs()*v.H2.Abs());
//                break;
//            case 5:
//                p.M[i][1]= n.re*v.H1.re + n.im*v.H2.re;
//                break;
//            case 6:
//                p.M[i][1]= t.re*v.H1.re + t.im*v.H2.re;
//                break;
//            default:
//                p.M[i][1]=0;
//                break;
//            }
//        }
//        else if (flag!=false)
//        {
//            switch (PlotType)
//            {
//            case 0:
//                p.M[i][1]=v.A.Abs();
//                p.M[i][2]=v.A.re;
//                p.M[i][3]=v.A.im;
//                break;
//            case 1:
//                p.M[i][1]=sqrt(v.B1.Abs()*v.B1.Abs() + v.B2.Abs()*v.B2.Abs());
//                break;
//            case 2:
//                p.M[i][2]= n.re*v.B1.re + n.im*v.B2.re;
//                p.M[i][3]= n.re*v.B1.im + n.im*v.B2.im;
//                p.M[i][1]=sqrt(p.M[i][2]*p.M[i][2] +
//                               p.M[i][3]*p.M[i][3]);
//                break;
//            case 3:
//                p.M[i][2]= t.re*v.B1.re + t.im*v.B2.re;
//                p.M[i][3]= t.re*v.B1.im + t.im*v.B2.im;
//                p.M[i][1]=sqrt(p.M[i][2]*p.M[i][2] +
//                               p.M[i][3]*p.M[i][3]);
//                break;
//            case 4:
//                p.M[i][1]=sqrt(v.H1.Abs()*v.H1.Abs() + v.H2.Abs()*v.H2.Abs());
//                break;
//            case 5:
//                p.M[i][2]= n.re*v.H1.re + n.im*v.H2.re;
//                p.M[i][3]= n.re*v.H1.im + n.im*v.H2.im;
//                p.M[i][1]=sqrt(p.M[i][2]*p.M[i][2] +
//                               p.M[i][3]*p.M[i][3]);
//                break;
//            case 6:
//                p.M[i][2]= t.re*v.H1.re + t.im*v.H2.re;
//                p.M[i][3]= t.re*v.H1.im + t.im*v.H2.im;
//                p.M[i][1]=sqrt(p.M[i][2]*p.M[i][2] +
//                               p.M[i][3]*p.M[i][3]);
//                break;
//            case 7:
//                p.M[i][2]= v.Je.re;
//                p.M[i][3]= v.Je.im;
//                p.M[i][1]= abs(v.Je);
//                break;
//            case 8:
//                p.M[i][2]= v.Je.re+v.Js.re;
//                p.M[i][3]= v.Je.im+v.Js.im;
//                p.M[i][1]= abs(v.Je+v.Js);
//                break;
//            default:
//                p.M[i][1]=0;
//                break;
//            }
//        }
//    }
//
//    free(q);
//}

bool FPProc::InTriangleTest(double x, double y, int i)
{
    int j,k;
    double z;
    bool InFlag;

    if(i<0) return false;

    for(j=0,InFlag=true; ((j<3) && (InFlag==true)); j++)
    {
        k=j+1;
        if(k==3) k=0;
        z=(meshnode[meshelem[i].p[k]].x-meshnode[meshelem[i].p[j]].x)*
          (y-meshnode[meshelem[i].p[j]].y) -
          (meshnode[meshelem[i].p[k]].y-meshnode[meshelem[i].p[j]].y)*
          (x-meshnode[meshelem[i].p[j]].x);
        if(z<0) InFlag=false;
    }

    return InFlag;
}

CPointVals::CPointVals()
{
    A=0;            // vector potential
    B1=0;
    B2=0;        // flux density
    mu1=1;
    mu2=1;    // permeability
    H1=0;
    H2=0;        // field intensity
    Je=0;
    Js=0;        // eddy current and source current densities
    c=0;            // conductivity
    E=0;            // energy stored in the magnetic field
    Ph=0;            // power dissipated by hysteresis
    Pe=0;
    ff=1;
}

CComplex FPProc::Ctr(int i)
{
    CComplex p,c;
    int j;

    for(j=0,c=0; j<3; j++)
    {
        p.Set(meshnode[ meshelem[i].p[j] ].x/3., meshnode[ meshelem[i].p[j] ].y/3.);
        c+=p;
    }

    return c;
}

double FPProc::ElmArea(int i)
{
    int j,n[3];
    double b0,b1,c0,c1;

    for(j=0; j<3; j++) n[j]=meshelem[i].p[j];

    b0=meshnode[n[1]].y - meshnode[n[2]].y;
    b1=meshnode[n[2]].y - meshnode[n[0]].y;
    c0=meshnode[n[2]].x - meshnode[n[1]].x;
    c1=meshnode[n[0]].x - meshnode[n[2]].x;
    return (b0*c1-b1*c0)/2.;

}

double FPProc::ElmArea(CElement *elm)
{
    int j,n[3];
    double b0,b1,c0,c1;

    for(j=0; j<3; j++) n[j]=elm->p[j];

    b0=meshnode[n[1]].y - meshnode[n[2]].y;
    b1=meshnode[n[2]].y - meshnode[n[0]].y;
    c0=meshnode[n[2]].x - meshnode[n[1]].x;
    c1=meshnode[n[0]].x - meshnode[n[2]].x;
    return (b0*c1-b1*c0)/2.;

}

CComplex FPProc::GetJA(int k,CComplex *J,CComplex *A)
{
    // returns current density with contribution from all sources in
    // units of MA/m^2

    int i,blk,lbl,crc;
    double r,c,rn;
    r=c=rn = 0;
    CComplex Javg;

    blk=meshelem[k].blk;
    lbl=meshelem[k].lbl;
    crc=blocklist[lbl].InCircuit;

    // first, get A
    for(i=0; i<3; i++)
    {
        if(ProblemType==PLANAR) A[i]=(meshnode[meshelem[k].p[i]].A);
        else
        {
            rn=meshnode[meshelem[k].p[i]].x*LengthConv[LengthUnits];
            if(fabs(rn/LengthConv[LengthUnits])<1.e-06) A[i]=0;
            else A[i]=(meshnode[meshelem[k].p[i]].A)/(2.*PI*rn);
        }
    }

    if(ProblemType==AXISYMMETRIC) r = Re(Ctr(k))*LengthConv[LengthUnits];

    // contribution from explicitly specified J
    for(i=0; i<3; i++) J[i]=blockproplist[blk].Jr+I*blockproplist[blk].Ji;
    Javg=blockproplist[blk].Jr+I*blockproplist[blk].Ji;

    c=blockproplist[blk].Cduct;
    if ((blockproplist[blk].Lam_d!=0) && (blockproplist[blk].LamType==0)) c=0;
    if (blocklist[lbl].FillFactor>0) c=0;


    // contribution from eddy currents;
    if(Frequency!=0)
        for(i=0; i<3; i++)
        {
            J[i]-=I*Frequency*2.*PI*c*A[i];
            Javg-=I*Frequency*2.*PI*c*A[i]/3.;
        }


    // contribution from circuit currents //
    if(crc>=0)
    {
        if(blocklist[lbl].Case==0)  // specified voltage
        {
            if(ProblemType==PLANAR)
            {
                for(i=0; i<3; i++)
                    J[i]-=c*blocklist[lbl].dVolts;
                Javg-=c*blocklist[lbl].dVolts;
            }
            else
            {
                for(i=0; i<3; i++)
                {
                    rn=meshnode[meshelem[k].p[i]].x;
                    if(fabs(rn/LengthConv[LengthUnits])<1.e-06)
                        J[i]-=c*blocklist[lbl].dVolts/r;
                    else
                        J[i]-=c*blocklist[lbl].dVolts/(rn*LengthConv[LengthUnits]);

                }
                Javg-=c*blocklist[lbl].dVolts/r;
            }
        }
        else
        {
            for(i=0; i<3; i++) J[i]+=blocklist[lbl].J; // specified current
            Javg+=blocklist[lbl].J;
        }

    }

    // convert results to A/m^2
    for(i=0; i<3; i++) J[i]*=1.e06;

    return (Javg*1.e06);
}

CComplex FPProc::PlnInt(double a, CComplex *u, CComplex *v)
{
    int i;
    CComplex z[3],x;

    z[0]=2.*u[0]+u[1]+u[2];
    z[1]=u[0]+2.*u[1]+u[2];
    z[2]=u[0]+u[1]+2.*u[2];

    for(i=0,x=0; i<3; i++) x+=v[i]*z[i];
    return a*x/12.;
}

CComplex FPProc::AxiInt(double a, CComplex *u, CComplex *v,double *r)
{
    int i;
    static CComplex M[3][3];
    CComplex x, z[3];

    M[0][0]=6.*r[0]+2.*r[1]+2.*r[2];
    M[0][1]=2.*r[0]+2.*r[1]+1.*r[2];
    M[0][2]=2.*r[0]+1.*r[1]+2.*r[2];
    M[1][1]=2.*r[0]+6.*r[1]+2.*r[2];
    M[1][2]=1.*r[0]+2.*r[1]+2.*r[2];
    M[2][2]=2.*r[0]+2.*r[1]+6.*r[2];
    M[1][0]=M[0][1];
    M[2][0]=M[0][2];
    M[2][1]=M[1][2];

    for(i=0; i<3; i++) z[i]=M[i][0]*u[0]+M[i][1]*u[1]+M[i][2]*u[2];
    for(i=0,x=0; i<3; i++) x+=v[i]*z[i];
    return PI*a*x/30.;
}

CComplex FPProc::HenrotteVector(int k)
{
    int i,n[3];
    double b[3],c[3],da;
    CComplex v;

    for(i=0; i<3; i++)
    {
        n[i] = meshelem[k].p[i];
    }

    b[0]=meshnode[n[1]].y - meshnode[n[2]].y;
    b[1]=meshnode[n[2]].y - meshnode[n[0]].y;
    b[2]=meshnode[n[0]].y - meshnode[n[1]].y;
    c[0]=meshnode[n[2]].x - meshnode[n[1]].x;
    c[1]=meshnode[n[0]].x - meshnode[n[2]].x;
    c[2]=meshnode[n[1]].x - meshnode[n[0]].x;

    da = (b[0] * c[1] - b[1] * c[0]);

    for(i=0,v=0; i<3; i++)
    {
        v -= meshnode[n[i]].msk * (b[i] + I * c[i]) / (da * LengthConv[LengthUnits]);  // grad
    }

    return v;
}

CComplex FPProc::BlockIntegral(int inttype)
{
    int i,k;
    CComplex c,y,z,J,mu1,mu2,B1,B2,H1,H2,F1,F2;
    CComplex A[3],Jn[3],U[3],V[3];
    double a,sig,R = 0;
    double r[3] = {0, 0, 0};

    z=0;
    for(i=0; i<3; i++) U[i]=1.;

    if(inttype==6) z= BlockIntegral(3) + BlockIntegral(4); //total losses
    else for(i=0; i<(int)meshelem.size(); i++)
        {
            if(blocklist[meshelem[i].lbl].IsSelected==true)
            {

                // compute some useful quantities employed by most integrals...
                J=GetJA(i,Jn,A);
                a=ElmArea(i)*pow(LengthConv[LengthUnits],2.);
                if(ProblemType==AXISYMMETRIC)
                {
                    for(k=0; k<3; k++)
                        r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
                    R=(r[0]+r[1]+r[2])/3.;
                }

                // now, compute the desired integral;
                switch(inttype)
                {
                case 0: //  A.J
                    for(k=0; k<3; k++) V[k]=Jn[k].Conj();
                    if(ProblemType==PLANAR)
                        y=PlnInt(a,A,V)*Depth;
                    else
                        y=AxiInt(a,A,V,r);
                    z+=y;

                    break;

                case 11: // x (or r) direction Lorentz force, SS part.
                    B2=meshelem[i].B2;
                    y= -(B2.re*J.re + B2.im*J.im);
                    if (ProblemType==AXISYMMETRIC) y=0;
                    else y*=Depth;
                    if(Frequency!=0) y*=0.5;
                    z+=(a*y);
                    break;

                case 12: // y (or z) direction Lorentz force, SS part.
                    for(k=0; k<3; k++) V[k]=Re(meshelem[i].B1*Jn[k].Conj());
                    if(ProblemType==PLANAR)
                        y=PlnInt(a,U,V)*Depth;
                    else
                        y=AxiInt(-a,U,V,r);
                    if(Frequency!=0) y*=0.5;
                    z+=y;

                    break;

                case 13: // x (or r) direction Lorentz force, 2x part.
                    if((Frequency!=0) && (ProblemType==PLANAR))
                    {
                        B2=meshelem[i].B2;
                        y= -(B2.re*J.re - B2.im*J.im) - I*(B2.re*J.im+B2.im*J.re);
                        z+=0.5*(a*y*Depth);
                    }
                    break;

                case 14: // y (or z) direction Lorentz force, 2x part.
                    if (Frequency!=0)
                    {
                        B1=meshelem[i].B1;
                        B2=meshelem[i].B2;
                        y= (B1.re*J.re - B1.im*J.im) + I*(B1.re*J.im+B1.im*J.re);
                        if(ProblemType==AXISYMMETRIC) y=(-y*2.*PI*R);
                        else y*=Depth;
                        z+=(a*y)/2.;
                    }
                    break;

                case 16: // Lorentz Torque, 2x
                    if ((Frequency!=0) && (ProblemType==PLANAR))
                    {
                        B1=meshelem[i].B1;
                        B2=meshelem[i].B2;
                        c=Ctr(i)*LengthConv[LengthUnits];
                        y= c.re*((B1.re*J.re - B1.im*J.im) + I*(B1.re*J.im+B1.im*J.re))
                           +c.im*((B2.re*J.re - B2.im*J.im) + I*(B2.re*J.im+B2.im*J.re));
                        z+=0.5*(a*y*Depth);
                    }
                    break;

                case 15: // Lorentz Torque, SS part.
                    if(ProblemType==PLANAR)
                    {
                        B1=meshelem[i].B1;
                        B2=meshelem[i].B2;
                        c=Ctr(i)*LengthConv[LengthUnits];
                        y= c.im*(B2.re*J.re + B2.im*J.im) + c.re*(B1.re*J.re + B1.im*J.im);
                        if(Frequency!=0) y*=0.5;
                        z+=(a*y*Depth);
                    }
                    break;

                case 1: // integrate A over the element;
                    if(ProblemType==AXISYMMETRIC)
                        y=AxiInt(a,U,A,r);
                    else
                        for(k=0,y=0; k<3; k++) y+=a*Depth*A[k]/3.;

                    z+=y;
                    break;

                case 2: // stored energy
                    if(ProblemType==AXISYMMETRIC) a*=(2.*PI*R);
                    else a*=Depth;
                    B1=meshelem[i].B1;
                    B2=meshelem[i].B2;
                    if(Frequency!=0)
                    {
                        // have to compute the energy stored in a special way for
                        // wound regions subject to prox and skin effects
                        if (blockproplist[meshelem[i].blk].LamType>2)
                        {
                            CComplex mu;
                            mu=muo*blocklist[meshelem[i].lbl].mu;
                            double u=Im(1./blocklist[meshelem[i].lbl].o)/(2.e6*PI*Frequency);
                            y=a*Re(B1*conj(B1)+B2*conj(B2))*Re(1./mu)/4.;
                            y+=a*Re(J*conj(J))*u/4.;
                        }
                        else y=a*blockproplist[meshelem[i].blk].DoEnergy(B1,B2);
                    }
                    else
                    {
                        // correct H and energy stored in magnet for second-quadrant
                        // representation of a PM.
                        if (blockproplist[meshelem[i].blk].H_c!=0)
                        {
                            int bk=meshelem[i].blk;

                            // in the linear case:
                            if (blockproplist[bk].BHpoints==0)
                            {
                                CComplex Hc;
                                mu1=blockproplist[bk].mu_x;
                                mu2=blockproplist[bk].mu_y;
                                H1=B1/(mu1*muo);
                                H2=B2/(mu2*muo);
                                Hc = blockproplist[bk].H_c*exp(I*PI*meshelem[i].magdir/180.);
                                H1=H1-Re(Hc);
                                H2=H2-Im(Hc);
                                y = a*0.5*muo*(mu1.re*H1.re*H1.re + mu2.re*H2.re*H2.re);
                            }
                            else  // the material is nonlinear
                            {
                                y=blockproplist[bk].DoEnergy(B1.re,B2.re);
                                y = y + blockproplist[bk].Nrg
                                    - blockproplist[bk].H_c*Re((B1.re+I*B2.re)/exp(I*PI*meshelem[i].magdir/180.));
                                y*=a;
                            }
                        }
                        else y=a*blockproplist[meshelem[i].blk].DoEnergy(B1.re,B2.re);

                        // add in "local" stored energy for wound that would be subject to
                        // prox and skin effect for nonzero frequency cases.
                        if (blockproplist[meshelem[i].blk].LamType>2)
                        {
                            double u=Im(blocklist[meshelem[i].lbl].o);
                            y+=a*Re(J*J)*u/2.;
                        }
                    }
                    y*=AECF(i); // correction for axisymmetric external region;

                    z+=y;
                    break;

                case 3:  // Hysteresis & Laminated eddy current losses
                    if(Frequency!=0)
                    {
                        if(ProblemType==AXISYMMETRIC) a*=(2.*PI*R);
                        else a*=Depth;
                        B1=meshelem[i].B1;
                        B2=meshelem[i].B2;
                        GetMu(B1,B2,mu1,mu2,i);
                        H1=B1/(mu1*muo);
                        H2=B2/(mu2*muo);

                        y=a*PI*Frequency*Im(H1*B1.Conj() + H2*B2.Conj());
                        z+=y;
                    }
                    break;

                case 4: // Resistive Losses
                    sig=1.e06/Re(1./blocklist[meshelem[i].lbl].o);
                    if((blockproplist[meshelem[i].blk].Lam_d!=0) &&
                            (blockproplist[meshelem[i].blk].LamType==0)) sig=0;
                    if(sig!=0)
                    {

                        if (ProblemType==PLANAR)
                        {
                            for(k=0; k<3; k++) V[k]=Jn[k].Conj()/sig;
                            y=PlnInt(a,Jn,V)*Depth;
                        }

                        if(ProblemType==AXISYMMETRIC)
                            y=2.*PI*R*a*J*conj(J)/sig;

                        if(Frequency!=0) y/=2.;
                        z+=y;
                    }
                    break;

                case 5: // cross-section area
                    z+=a;
                    break;

                case 10: // volume
                    if(ProblemType==AXISYMMETRIC) a*=(2.*PI*R);
                    else a*=Depth;
                    z+=a;
                    break;

                case 7: // total current in block;
                    z+=a*J;

                    break;

                case 8: // integrate x or r part of b over the block
                    if(ProblemType==AXISYMMETRIC) a*=(2.*PI*R);
                    else a*=Depth;
                    z+=(a*meshelem[i].B1);
                    break;

                case 9: // integrate y or z part of b over the block
                    if(ProblemType==AXISYMMETRIC) a*=(2.*PI*R);
                    else a*=Depth;
                    z+=(a*meshelem[i].B2);
                    break;

                case 17: // Coenergy
                    if(ProblemType==AXISYMMETRIC) a*=(2.*PI*R);
                    else a*=Depth;
                    B1=meshelem[i].B1;
                    B2=meshelem[i].B2;
                    if(Frequency!=0)
                    {
                        // have to compute the energy stored in a special way for
                        // wound regions subject to prox and skin effects
                        if (blockproplist[meshelem[i].blk].LamType>2)
                        {
                            CComplex mu;
                            mu=muo*blocklist[meshelem[i].lbl].mu;
                            double u=Im(1./blocklist[meshelem[i].lbl].o)/(2.e6*PI*Frequency);
                            y=a*Re(B1*conj(B1)+B2*conj(B2))*Re(1./mu)/4.;
                            y+=a*Re(J*conj(J))*u/4.;
                        }
                        else y=a*blockproplist[meshelem[i].blk].DoCoEnergy(B1,B2);
                    }
                    else
                    {
                        y=a*blockproplist[meshelem[i].blk].DoCoEnergy(B1.re,B2.re);

                        // add in "local" stored energy for wound that would be subject to
                        // prox and skin effect for nonzero frequency cases.
                        if (blockproplist[meshelem[i].blk].LamType>2)
                        {
                            double u=Im(blocklist[meshelem[i].lbl].o);
                            y+=a*Re(J*J)*u/2.;
                        }
                    }
                    y*=AECF(i); // correction for axisymmetric external region;

                    z+=y;
                    break;

                case 24: // Moment of Inertia-like integral

                    // For axisymmetric problems, compute the moment
                    // of inertia about the r=0 axis.
                    if(ProblemType==AXISYMMETRIC)
                    {
                        for(k=0; k<3; k++) V[k]=r[k];
                        y=AxiInt(a,V,V,r);
                    }

                    // For planar problems, compute the moment of
                    // inertia about the z=axis.
                    else
                    {
                        for(k=0; k<3; k++)
                        {
                            U[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
                            V[k]=meshnode[meshelem[i].p[k]].y*LengthConv[LengthUnits];
                        }
                        y =U[0]*U[0] + U[1]*U[1] + U[2]*U[2];
                        y+=U[0]*U[1] + U[0]*U[2] + U[1]*U[2];
                        y+=V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
                        y+=V[0]*V[1] + V[0]*V[2] + V[1]*V[2];
                        y*=(a*Depth/6.);
                    }

                    z+=y;
                    break;

                default:
                    break;
                }
            }

CComplex temp;
            // integrals that need to be evaluated over all elements,
            // regardless of which elements are actually selected.
            if((inttype>=18) || (inttype<=23))
            {
                a=ElmArea(i)*pow(LengthConv[LengthUnits],2.);
                if(ProblemType==AXISYMMETRIC)
                {
                    for(k=0; k<3; k++)
                        r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
                    R=(r[0]+r[1]+r[2])/3.;
                    a*=(2.*PI*R);
                }
                else a*=Depth;

                switch(inttype)
                {

                case 18: // x (or r) direction Henrotte force, SS part.
                    if(ProblemType!=0) break;

                    B1 = meshelem[i].B1;

                    B2 = meshelem[i].B2;

                    c = HenrotteVector(i);

                    y = (((B1*conj(B1)) - (B2*conj(B2)))*Re(c) + 2.*Re(B1*conj(B2))*Im(c))/(2.*muo);

                    if(Frequency!=0)
                    {
                        y/=2.;
                    }

                    y*=AECF(i); // correction for axisymmetric external region;

                    z+=(a*y);
                    break;

                case 19: // y (or z) direction Henrotte force, SS part.

                    B1=meshelem[i].B1;
                    B2=meshelem[i].B2;
                    c=HenrotteVector(i);

                    y=(((B2*conj(B2)) - (B1*conj(B1)))*Im(c) + 2.*Re(B1*conj(B2))*Re(c))/(2.*muo);

                    y*=AECF(i); // correction for axisymmetric external region;

                    if(Frequency!=0) y/=2.;
                    z+=(a*y);

                    break;

                case 20: // x (or r) direction Henrotte force, 2x part.

                    if(ProblemType!=0) break;
                    B1=meshelem[i].B1;
                    B2=meshelem[i].B2;
                    c=HenrotteVector(i);
                    z+=a*((((B1*B1) - (B2*B2))*Re(c) + 2.*B1*B2*Im(c))/(4.*muo)) * AECF(i);

                    break;

                case 21: // y (or z) direction Henrotte force, 2x part.

                    B1=meshelem[i].B1;
                    B2=meshelem[i].B2;
                    c=HenrotteVector(i);
                    z+= a*((((B2*B2) - (B1*B1))*Im(c) + 2.*B1*B2*Re(c))/(4.*muo)) * AECF(i);

                    break;

                case 22: // Henrotte torque, SS part.
                    if(ProblemType!=0) break;
                    B1=meshelem[i].B1;
                    B2=meshelem[i].B2;
                    c=HenrotteVector(i);

                    F1 = (((B1*conj(B1)) - (B2*conj(B2)))*Re(c) +
                          2.*Re(B1*conj(B2))*Im(c))/(2.*muo);
                    F2 = (((B2*conj(B2)) - (B1*conj(B1)))*Im(c) +
                          2.*Re(B1*conj(B2))*Re(c))/(2.*muo);

                    for(c=0,k=0; k<3; k++)
                        c+=meshnode[meshelem[i].p[k]].CC()*LengthConv[LengthUnits]/3.;

                    y=Re(c)*F2 -Im(c)*F1;
                    if(Frequency!=0) y/=2.;
                    y*=AECF(i);
                    z+=(a*y);

                    break;

                case 23: // Henrotte torque, 2x part.

                    if(ProblemType!=0) break;
                    B1=meshelem[i].B1;
                    B2=meshelem[i].B2;
                    c=HenrotteVector(i);
                    F1 = (((B1*B1) - (B2*B2))*Re(c) + 2.*B1*B2*Im(c))/(4.*muo);
                    F2 = (((B2*B2) - (B1*B1))*Im(c) + 2.*B1*B2*Re(c))/(4.*muo);

                    for(c=0,k=0; k<3; k++)
                        c+=meshnode[meshelem[i].p[k]].CC()*LengthConv[LengthUnits]/3;

                    z+=a*(Re(c)*F2 -Im(c)*F1)*AECF(i);

                    break;

                default:
                    break;
                }
            }
        }

    return z;
}

void FPProc::LineIntegral(int inttype, CComplex *z)
{
// inttype    Integral
//        0    B.n
//        1    H.t
//        2    Contour length
//        3    Stress Tensor Force
//        4    Stress Tensor Torque
//        5    (B.n)^2

    // inttype==0 => B.n
    if(inttype==0)
    {
        CComplex a0,a1;
        CPointVals u;
        double l;
        int i,k;

        k=contour.size();
        GetPointValues(contour[0].re,contour[0].im, u);
        a0=u.A;
        GetPointValues(contour[k-1].re,contour[k-1].im,u);
        a1=u.A;
        if(ProblemType==PLANAR)
        {
            for(i=0,l=0; i<k-1; i++)
                l+=abs(contour[i+1]-contour[i]);
            l*=LengthConv[LengthUnits];
            z[0] = (a0-a1)*Depth;
            if(l!=0) z[1]= z[0]/(l*Depth);
        }
        else
        {
            for(i=0,l=0; i<k-1; i++)
                l+=(PI*(contour[i].re+contour[i+1].re)*
                    abs(contour[i+1]-contour[i]));
            l*=pow(LengthConv[LengthUnits],2.);
            z[0]= a1-a0;
            if(l!=0) z[1]= z[0]/l;
        }
    }

    // inttype==1 => H.t
    if(inttype==1)
    {
        CComplex n,t,pt,Ht;
        CPointVals v;
        double dz,u,l;
        int i,j,k,m,elm;
        int NumPlotPoints=d_LineIntegralPoints;
        bool flag;

        z[0]=0;
        for(k=1; k<(int)contour.size(); k++)
        {
            dz=abs(contour[k]-contour[k-1])/((double) NumPlotPoints);
            for(i=0,elm=-1; i<NumPlotPoints; i++)
            {
                u=(((double) i)+0.5)/((double) NumPlotPoints);
                pt=contour[k-1] + u*(contour[k] - contour[k-1]);
                t=contour[k]-contour[k-1];
                t/=abs(t);
                n=I*t;
                pt+=n*1.e-06;

                if (elm<0) elm=InTriangle(pt.re,pt.im);
                else if (InTriangleTest(pt.re,pt.im,elm)==false)
                {
                    flag=false;
                    for(j=0; j<3; j++)
                        for(m=0; m<NumList[meshelem[elm].p[j]]; m++)
                        {
                            elm=ConList[meshelem[elm].p[j]][m];
                            if (InTriangleTest(pt.re,pt.im,elm)==true)
                            {
                                flag=true;
                                m=100;
                                j=3;
                            }
                        }
                    if (flag==false) elm=InTriangle(pt.re,pt.im);
                }
                if(elm>=0)
                    flag=GetPointValues(pt.re,pt.im,elm,v);
                else flag=false;

                if(flag==true)
                {
                    Ht = t.re*v.H1 + t.im*v.H2;
                    z[0]+=(Ht*dz*LengthConv[LengthUnits]);
                }
            }


            for(i=0,l=0; i<(int)contour.size()-1; i++)
                l+=abs(contour[i+1]-contour[i]);
            l*=LengthConv[LengthUnits];
            if(l!=0) z[1]=z[0]/l;
        }
    }

    // inttype==2 => Contour Length
    if(inttype==2)
    {
        int i,k;
        k=contour.size();
        for(i=0,z[0].re=0; i<k-1; i++)
            z[0].re+=abs(contour[i+1]-contour[i]);
        z[0].re*=LengthConv[LengthUnits];

        if(ProblemType==AXISYMMETRIC)
        {
            for(i=0,z[0].im=0; i<k-1; i++)
                z[0].im+=(PI*(contour[i].re+contour[i+1].re)*
                          abs(contour[i+1]-contour[i]));
            z[0].im*=pow(LengthConv[LengthUnits],2.);
        }
        else
        {
            z[0].im=z[0].re*Depth;
        }
    }

    // inttype==3 => Stress Tensor Force
    if(inttype==3)
    {
        CComplex n,t,pt,Hn,Bn,BH,dF1,dF2;
        CPointVals v;
        double dz,dza,u;
        int i,j,k,m,elm;
        int NumPlotPoints=d_LineIntegralPoints;
        bool flag;

        for(i=0; i<4; i++) z[i]=0;

        for(k=1; k<(int)contour.size(); k++)
        {
            dz=abs(contour[k]-contour[k-1])/((double) NumPlotPoints);
            for(i=0,elm=-1; i<NumPlotPoints; i++)
            {
                u=(((double) i)+0.5)/((double) NumPlotPoints);
                pt=contour[k-1] + u*(contour[k] - contour[k-1]);
                t=contour[k]-contour[k-1];
                t/=abs(t);
                n=I*t;
                pt+=n*1.e-06;

                if (elm<0) elm=InTriangle(pt.re,pt.im);
                else if (InTriangleTest(pt.re,pt.im,elm)==false)
                {
                    flag=false;
                    for(j=0; j<3; j++)
                        for(m=0; m<NumList[meshelem[elm].p[j]]; m++)
                        {
                            elm=ConList[meshelem[elm].p[j]][m];
                            if (InTriangleTest(pt.re,pt.im,elm)==true)
                            {
                                flag=true;
                                m=100;
                                j=3;
                            }
                        }
                    if (flag==false) elm=InTriangle(pt.re,pt.im);
                }
                if(elm>=0)
                    flag=GetPointValues(pt.re,pt.im,elm,v);
                else flag=false;

                if(flag==true)
                {
                    if(Frequency==0)
                    {
                        Hn= n.re*v.H1 + n.im*v.H2;
                        Bn= n.re*v.B1 + n.im*v.B2;
                        BH= v.B1*v.H1 + v.B2*v.H2;
                        dF1=v.H1*Bn + v.B1*Hn - n.re*BH;
                        dF2=v.H2*Bn + v.B2*Hn - n.im*BH;

                        dza=dz*LengthConv[LengthUnits];
                        if(ProblemType==AXISYMMETRIC)
                        {
                            dza*=2.*PI*pt.re*LengthConv[LengthUnits];
                            dF1=0;
                        }
                        else dza*=Depth;

                        z[0]+=(dF1*dza/2.);
                        z[1]+=(dF2*dza/2.);
                    }
                    else
                    {
                        Hn=n.re*v.H1 + n.im*v.H2;
                        Bn=n.re*v.B1 + n.im*v.B2;
                        BH = v.B1*v.H1 + v.B2*v.H2;
                        dF1 = v.H1*Bn + v.B1*Hn - n.re*BH;
                        dF2 = v.H2*Bn + v.B2*Hn - n.im*BH;

                        dza=dz*LengthConv[LengthUnits];
                        if(ProblemType==AXISYMMETRIC)
                        {
                            dza*=2.*PI*pt.re*LengthConv[LengthUnits];
                            dF1=0;
                        }
                        else dza*=Depth;

                        z[0]+=(dF1*dza/4.);
                        z[1]+=(dF2*dza/4.);

                        BH  = v.B1*v.H1.Conj() +v.B2*v.H2.Conj();

                        if (ProblemType!=AXISYMMETRIC)
                            dF1 = v.H1*Bn.Conj() + v.B1*Hn.Conj() - n.re*BH;
                        dF2=  v.H2*Bn.Conj() + v.B2*Hn.Conj() - n.im*BH;


                        z[2]+=(dF1*dza/4.);
                        z[3]+=(dF2*dza/4.);
                    }
                }
            }

        }
    }

    // inttype==4 => Stress Tensor Torque
    if(inttype==4)
    {
        CComplex n,t,pt,Hn,Bn,BH,dF1,dF2,dT;
        CPointVals v;
        double dz,dza,u;
        int i,j,k,m,elm;
        int NumPlotPoints=d_LineIntegralPoints;
        bool flag;

        for(i=0; i<2; i++) z[i].Set(0,0);

        for(k=1; k<(int)contour.size(); k++)
        {
            dz=abs(contour[k]-contour[k-1])/((double) NumPlotPoints);
            for(i=0,elm=-1; i<NumPlotPoints; i++)
            {
                u=(((double) i)+0.5)/((double) NumPlotPoints);
                pt=contour[k-1]+ u*(contour[k] - contour[k-1]);
                t=contour[k]-contour[k-1];
                t/=abs(t);
                n=I*t;
                pt+=n*1.e-6;

                if (elm<0) elm=InTriangle(pt.re,pt.im);
                else if (InTriangleTest(pt.re,pt.im,elm)==false)
                {
                    flag=false;
                    for(j=0; j<3; j++)
                        for(m=0; m<NumList[meshelem[elm].p[j]]; m++)
                        {
                            elm=ConList[meshelem[elm].p[j]][m];
                            if (InTriangleTest(pt.re,pt.im,elm)==true)
                            {
                                flag=true;
                                m=100;
                                j=3;
                            }
                        }
                    if (flag==false) elm=InTriangle(pt.re,pt.im);
                }

                if(elm>=0)
                {
                    flag=GetPointValues(pt.re,pt.im,elm,v);
                }
                else
                {
                    flag=false;
                }

                if(flag==true)
                {
                    if(Frequency==0)
                    {
                        Hn= n.re*v.H1 + n.im*v.H2;
                        Bn= n.re*v.B1 + n.im*v.B2;
                        BH= v.B1*v.H1 + v.B2*v.H2;
                        dF1=v.H1*Bn + v.B1*Hn - n.re*BH;
                        dF2=v.H2*Bn + v.B2*Hn - n.im*BH;
                        dT= pt.re*dF2 - dF1*pt.im;
                        dza=dz*LengthConv[LengthUnits]*LengthConv[LengthUnits];

                        z[0]+=(dT*dza*Depth/2.);
                    }
                    else
                    {
                        Hn=n.re*v.H1 + n.im*v.H2;
                        Bn=n.re*v.B1 + n.im*v.B2;
                        BH = v.B1*v.H1 + v.B2*v.H2;
                        dF1 = v.H1*Bn + v.B1*Hn - n.re*BH;
                        dF2 = v.H2*Bn + v.B2*Hn - n.im*BH;
                        dT=pt.re*dF2 - dF1*pt.im;
                        dza=dz*LengthConv[LengthUnits]*LengthConv[LengthUnits];

                        z[0]+=(dT*dza*Depth/4.);

                        BH  = v.B1*v.H1.Conj() +v.B2*v.H2.Conj();
                        dF1 = v.H1*Bn.Conj() + v.B1*Hn.Conj() - n.re*BH;
                        dF2=  v.H2*Bn.Conj() + v.B2*Hn.Conj() - n.im*BH;
                        dT= pt.re*dF2 - dF1*pt.im ;

                        z[1]+=(dT*dza*Depth/4.);

                    }
                }
            }
        }

    }

    // inttype==5 => (B.n)^2
    if(inttype==5)
    {
        CComplex n,t,pt,Ht;
        CPointVals pvals;
        double dz,u,l;
        int i,j,k,m,elm;
        int NumPlotPoints = d_LineIntegralPoints;
        bool flag;

        z[0] = 0;
        // loop through each segment in the contour intgrating over each in turn
        for(k=1; k<(int)contour.size(); k++)
        {
            // break the segment up into differential subsegments for the integral
            dz = abs(contour[k]-contour[k-1]) / ((double) NumPlotPoints);

            // loop through each subsegment performing the integral
            for(i=0,elm=-1; i<NumPlotPoints; i++)
            {
                // get a point location in the middle of the subsegment
                u = (((double) i) + 0.5) / ((double) NumPlotPoints);
                pt = contour[k-1] + u*(contour[k] - contour[k-1]);
                // get a unit vector tangential to the segment at the sample point
                t = contour[k]-contour[k-1];
                t /= abs(t);
                // get a unit vector normal to the segment at the sample point
                n = I * t;
                // shift the sample point a little in the normal direction
                pt += n * 1.e-06;

                if (elm < 0)
                {
                    // This is the first run and we must find which mesh element we are in
                    elm = InTriangle(pt.re,pt.im);
                }
                else if (InTriangleTest(pt.re,pt.im,elm) == false)
                {
                    // This is not the first run, but we are no longer in the same element
                    // and must rediscover what element we are in
                    flag = false;
                    // first check neighbouring elements to the last element as it will
                    // save time over checking the whole mesh
                    for(j=0; j<3; j++)
                    {
                        for(m=0; m < NumList[meshelem[elm].p[j]]; m++)
                        {
                            elm = ConList[meshelem[elm].p[j]][m];

                            if (InTriangleTest(pt.re,pt.im,elm) == true)
                            {
                                // the current element was a neighbour of the previous element
                                flag = true;
                                m = 100;
                                j = 3;
                            }
                        }
                    }

                    if (flag == false)
                    {
                        // new element was not a neighbour of the old element, so we must
                        // search the whole mesh
                        elm = InTriangle(pt.re,pt.im);
                    }
                } // if (elm < 0)

                if (elm >= 0)
                {
                    // Get the point values at the sample location
                    flag = GetPointValues(pt.re,pt.im,elm,pvals);
                }
                else
                {
                    flag = false;
                }

                if(flag == true)
                {
                    // get the dot product of the normal and the B field at the sample point
                    Ht = n.re * pvals.B1 + n.im * pvals.B2;
                    // add the square of the field times the
                    z[0] += (Ht * Ht.Conj() * dz * LengthConv[LengthUnits]);
                }

            } // for(i=0,elm=-1; i<NumPlotPoints; i++)

            // now we will also calculate the average over the contour
            for(i=0,l=0; i<(int)contour.size()-1; i++)
            {
                // add up the contour segment lengths
                l += abs(contour[i+1]-contour[i]);
            }
            // convert the length units
            l *= LengthConv[LengthUnits];

            if(l!=0)
            {
                // divide the integral by the contour length
                z[1] = z[0] / l;
            }

        }
    }

    return;
}


int FPProc::ClosestArcSegment(double x, double y)
{
    double d0,d1;
    int i,j;

    if(arclist.size()==0) return -1;

    j=0;
    d0=ShortestDistanceFromArc(CComplex(x,y),arclist[0]);
    for(i=0; i<(int)arclist.size(); i++)
    {
        d1=ShortestDistanceFromArc(CComplex(x,y),arclist[i]);
        if(d1<d0)
        {
            d0=d1;
            j=i;
        }
    }

    return j;
}

void FPProc::GetCircle(CArcSegment &arc,CComplex &c, double &R)
{
    CComplex a0,a1,t;
    double d,tta;

    a0.Set(nodelist[arc.n0].x,nodelist[arc.n0].y);
    a1.Set(nodelist[arc.n1].x,nodelist[arc.n1].y);
    d=abs(a1-a0);            // distance between arc endpoints

    // figure out what the radius of the circle is...
    t=(a1-a0)/d;
    tta=arc.ArcLength*PI/180.;
    R=d/(2.*sin(tta/2.));
    c=a0 + (d/2. + I*sqrt(R*R-d*d/4.))*t; // center of the arc segment's circle...
}

double FPProc::ShortestDistanceFromArc(CComplex p, CArcSegment &arc)
{
    double R,d,l,z;
    CComplex a0,a1,c,t;

    a0.Set(nodelist[arc.n0].x,nodelist[arc.n0].y);
    a1.Set(nodelist[arc.n1].x,nodelist[arc.n1].y);
    GetCircle(arc,c,R);
    d=abs(p-c);

    if(d==0) return R;

    t=(p-c)/d;
    l=abs(p-c-R*t);
    z=arg(t/(a0-c))*180/PI;
    if ((z>0) && (z<arc.ArcLength)) return l;

    z=abs(p-a0);
    l=abs(p-a1);
    if(z<l) return z;
    return l;
}

double FPProc::ShortestDistanceFromSegment(double p, double q, int segm)
{
    double t,x[3],y[3];

    x[0]=nodelist[linelist[segm].n0].x;
    y[0]=nodelist[linelist[segm].n0].y;
    x[1]=nodelist[linelist[segm].n1].x;
    y[1]=nodelist[linelist[segm].n1].y;

    t=((p-x[0])*(x[1]-x[0]) + (q-y[0])*(y[1]-y[0]))/
      ((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]));

    if (t>1.) t=1.;
    if (t<0.) t=0.;

    x[2]=x[0]+t*(x[1]-x[0]);
    y[2]=y[0]+t*(y[1]-y[0]);

    return sqrt((p-x[2])*(p-x[2]) + (q-y[2])*(q-y[2]));
}

// bool FPProc::ScanPreferences()
// {
//     FILE *fp;
//     string fname;
//
//     fname= (string)BinDir + "femmview.cfg";
//
//     d_LineIntegralPoints = 100;
//     d_ShiftH =
//
//     fp=fopen(fname,"rt");
//     if (fp!=NULL)
//     {
//         bool flag=false;
//         char s[1024];
//         char q[1024];
//         char *v;
//
//         // parse the file
//         while (fgets(s,1024,fp)!=NULL)
//         {
//             sscanf(s,"%s",q);
//
//             if( _strnicmp(q,"<LineIntegralPoints>",20)==0)
//             {
//               v=StripKey(s);
//               sscanf(v,"%i",&d_LineIntegralPoints);
//               q[0] = '\0';
//             }
//
//             if( _strnicmp(q,"<ShiftH>",8)==0)
//             {
//               v=StripKey(s);
//               sscanf(v,"%i",&d_ShiftH);
//               q[0] = '\0';
//             }
//         }
//         fclose(fp);
//     }
//     else return false;
//
//     return true;
// }

void FPProc::BendContour(double angle, double anglestep)
{
    if (angle==0) return;
    if (anglestep==0) anglestep=1;

    int k,n;
    double d,tta,dtta,R;
    CComplex c,a0,a1;

    // check to see if there are at least enough
    // points to have made one line;
    k = contour.size()-1;
    if (k<1) return;

    // restrict the angle of the contour to 180 degrees;
    if ((angle<-180.) || (angle>180.))
    {
        return;
    }
    n = (int) ceil(fabs(angle/anglestep));
    tta = angle*PI/180.;
    dtta = tta/((double) n);

    // pop last point off of the contour;
    a1 = contour[k];
    contour.erase(contour.begin()+k);
    a0 = contour[k-1];

    // compute location of arc center;
    // and radius of the circle that the
    // arc lives on.
    d = abs(a1-a0);
    R = d / ( 2. * sin(fabs(tta/2.)) );

    if(tta>0)
    {
        c = a0 + (R/d) * (a1-a0) * exp(I*(PI-tta)/2.);
    }
    else
    {
        c = a0 + (R/d) * (a1-a0) * exp(-I*(PI+tta)/2.);
    }

    // add the points on the contour
    for(k=1; k<=n; k++)
    {
        contour.push_back( c + (a0 - c) * exp(k * I * dtta) );
    }
}


//bool FPProc::OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
//{
//    // TODO: Add your specialized code here and/or call the base class
//    if (bLinehook!=false) return true;
//    return CDocument::OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
//}

CComplex FPProc::GetStrandedVoltageDrop(int lbl)
{
    // Derive the voltage drop associated with a stranded and
    // current-carrying region.

    int i,k;
    CComplex dVolts,rho;
    CComplex A[3],J[3],U[3],V[3];
    double a,atot;
    double r[3];

    U[0]=1;
    U[1]=1;
    U[2]=1;

    for(i=0,dVolts=0,atot=0; i<(int)meshelem.size(); i++)
    {
        if(meshelem[i].lbl==lbl)
        {
            rho=blocklist[meshelem[i].lbl].o*1.e6;
            if(Frequency==0) rho=Re(rho);
            if (rho!=0) rho=(1./rho);

            GetJA(i,J,A);
            a=ElmArea(i)*LengthConv[LengthUnits]*LengthConv[LengthUnits];
            atot+=a;

            if(ProblemType==AXISYMMETRIC)
            {
                for(k=0; k<3; k++)
                    r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
            }

            for(k=0; k<3; k++) V[k]=(2.*PI*I*Frequency*A[k] + rho*J[k]);
            if(ProblemType==PLANAR) dVolts+=PlnInt(a,V,U)*Depth;
            else dVolts+=AxiInt(a,V,U,r);
        }
    }
    dVolts*=( ((double) blocklist[lbl].Turns) / atot);

    return dVolts;
}


void FPProc::GetFillFactor(int lbl)
{
    // Get the fill factor associated with a stranded and
    // current-carrying region.  For AC problems, also compute
    // the apparent conductivity and permeability for use in
    // post-processing the voltage.

    CMaterialProp* bp= &blockproplist[blocklist[lbl].BlockType];
    CBlockLabel* bl= &blocklist[lbl];
    double lc=LengthConv[LengthUnits]*LengthConv[LengthUnits];
    double atot,awire,w,d,o,fill,dd,W,R,c1,c2,c3,c4;
    atot=awire=w=d=o=fill=dd=W=R=c1=c2=c3=c4=0;
    int i,wiretype;
    CComplex ufd,ueff,ofd;

    // default values
    if (abs(bl->Turns)>1)
        bl->FillFactor=1;
    else
        bl->FillFactor=-1;
    bl->o=bp->Cduct;
    bl->mu=0.;

    if (blockproplist[blocklist[lbl].BlockType].LamType<3) return;

    // compute total area of associated block
    for(i=0,atot=0; i<(int)meshelem.size(); i++)
        if(meshelem[i].lbl==lbl) atot+=ElmArea(i)*lc;
    if (atot==0) return;

    wiretype=bp->LamType-3;
    // wiretype = 0 for magnet wire
    // wiretype = 1 for stranded but non-litz wire
    // wiretype = 2 for litz wire
    // wiretype = 3 for rectangular wire

    if(wiretype==3) // rectangular wire
    {
        W=2.*PI*Frequency;
        d=bp->WireD*0.001;
        bl->FillFactor=fabs(d*d*((double) bl->Turns)/atot);
        dd=d/sqrt(bl->FillFactor);    // foil pitch
        fill=d/dd;                    // fill for purposes of equivalent foil analysis
        o=bp->Cduct*(d/dd)*1.e6;    // effective foil conductivity in S/m

        if(Frequency==0)
        {
            bl->o  = bp->Cduct*bl->FillFactor + I*(dd-d)*dd*muo/6.;
            // for frequency=0 problems, imaginary part of conductivity
            // is used to store a factor used for local stored energy computation.
            // The factor is equal to Im[Normal[Series[1/sigma,{W,0,1}]]]/W
            bl->mu = 1;
            return;
        }

        // effective permeability for the equivalent foil.  Note that this is
        // the same equation as effective permeability of a lamination...
        if (o!=0)
        {
            ufd=muo*tanh(sqrt(I*W*o*muo)*d/2.)/(sqrt(I*W*o*muo)*d/2.);
            ueff=(fill*ufd+(1.-fill)*muo);
            bl->o= 1./(muo/(fill*o*ufd) + I*dd*dd*(1.-fill)*muo*W/4. - I*dd*dd*ueff*W/12.);
            bl->o*=1.e-6; // represent conductivity in units of MS/m for consistency with other parts of code.
            bl->mu=ueff/muo;
        }
        else
        {
            // This is a non-physical case because there is current but conductivity is zero;
            // Treat by idealizing as an imaginary conductivity that takes account of the
            // locally stored energy in the winding due to the disribution of turns.
            bl->mu=1;
            bl->o=6./(I*W*(dd-d)*dd*muo);
        }
        return;
    }


    // procedure for round wires;

    switch (wiretype)
    {
        // wiretype = 0 for magnet wire
    case 0:
        R=bp->WireD*0.0005;
        awire=PI*R*R*((double) bp->NStrands)*((double) bl->Turns);
        break;

        // wiretype = 1 for stranded but non-litz wire
    case 1:
        R=bp->WireD*0.0005*sqrt((double) bp->NStrands);
        awire=PI*R*R*((double) bl->Turns);
        break;

        // wiretype = 2 for litz wire
    case 2:
        R=bp->WireD*0.0005;
        awire=PI*R*R*((double) bp->NStrands)*((double) bl->Turns);
        break;
    }
    bl->FillFactor=fabs(awire/atot);
    fill=bl->FillFactor;

    // preliminary definitions
    w=2.*PI*Frequency;                        // frequency in rad/s
    o=bp->Cduct*1.e6;                        // conductivity in S/m
    W=w*o*muo*R*R/2.;                        // non-dimensionalized frequency
    dd=(1.6494541661869013*R)/sqrt(fill);    // foil pitch in equivalent foil geometry


    if(Frequency==0)
    {
        bl->o  = bp->Cduct*fill
                 + ((I/2.)*muo*R*R*log(1.5299240194394943/sqrt(fill)))/fill
                 - ((I/12.)*muo*dd*dd);
        // for frequency=0 problems, the imaginary part of conductivity
        // is used to store a factor used for local stored energy computation.
        bl->mu = 1.; // relative permeability of the block.
        return;
    }

    if (bp->Cduct==0)
    {
        // This is a non-physical case because there is current but conductivity is zero;
        // Treat by idealizing as an imaginary conductivity that takes account of the
        // locally stored energy in the winding due to the disribution of turns.
        bl->o    = 1./(((I/2.)*w*muo*R*R*log(1.5299240194394943/sqrt(fill)))/fill - (I/12.)*muo*dd*dd);
        bl->mu    = 1;
        return;
    }

    // fit for frequency-dependent permeability...
    c1=0.7756067409818643 + fill*(0.6873854335408803 + fill*(0.06841584481674128 -0.07143732702512284*fill));
    c2=1.5*fill/c1;
    ufd=c2*(tanh(sqrt(c1*I*W))/sqrt(c1*I*W))+(1.-c2); // relative frequency-dependent permeability
    bl->mu=ufd;

    // fit for frequency-dependent conductivity....
    c3=0.8824642871525136+fill*(-0.008605512994838827+fill*(0.7223208744682307-0.2157183942377177*fill));
    c4=log(1.5299240194394943/sqrt(fill))-c3/3.;
    ofd=o*fill/(I*c4*W+sqrt(I*c3*W)*(1./tanh(sqrt(I*c3*W))));    // fit to curves in Skin4.nb
    ofd=1./(1./ofd-I*w*ufd*muo*dd*dd/12.);                        // don't double-book local stored energy;
    bl->o=ofd*1.e-6;                                            // return frequency-dependent conductivity in MS/m
}

CComplex FPProc::GetStrandedLinkage(int lbl)
{
    // This is a routine for the special case of determining
    // the flux linkage of a stranded conductor at zero frequency
    // when the conductor is carrying zero current.
    int i,k;
    CComplex FluxLinkage;
    CComplex A[3],J[3],U[3];//,V[3];
    double a,atot;
    double r[3];

    U[0]=1;
    U[1]=1;
    U[2]=1;

    for(i=0,FluxLinkage=0,atot=0; i<(int)meshelem.size(); i++)
    {
        if(meshelem[i].lbl==lbl)
        {
            GetJA(i,J,A);
            a=ElmArea(i)*LengthConv[LengthUnits]*LengthConv[LengthUnits];
            atot+=a;

            if(ProblemType==AXISYMMETRIC)
            {
                for(k=0; k<3; k++)
                    r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
            }

            if(ProblemType==PLANAR) FluxLinkage+=PlnInt(a,A,U)*Depth;
            else FluxLinkage+=AxiInt(a,A,U,r);
        }
    }
    FluxLinkage*=( ((double) blocklist[lbl].Turns) / atot);

    return FluxLinkage;
}

CComplex FPProc::GetSolidAxisymmetricLinkage(int lbl)
{
    // This is a routine for the special case of determining
    // the flux linkage of a solid and axisymmetric conductor
    // t zero frequency when the conductor is carrying zero
    // current.  The trick here is to take account of the distribution
    // of the current that would be there, if there was any current.
    // In solid axisymmetric regions, the inner radius of the
    // region tends to carry a higher current density that the outer
    // edges, because the length of conductor that the current has
    // to traverse is smaller on the inner edge.

    int i,k;
    CComplex FluxLinkage;
    CComplex Aa,A[3],J[3],U[3];//,V[3];
    double a,atot,R;
    double r[3];

    U[0]=1;
    U[1]=1;
    U[2]=1;

    for(i=0,FluxLinkage=0,atot=0; i<(int)meshelem.size(); i++)
    {
        if(meshelem[i].lbl==lbl)
        {
            GetJA(i,J,A);
            Aa=(A[0]+A[1]+A[2])/3.;
            a=ElmArea(i)*LengthConv[LengthUnits]*LengthConv[LengthUnits];

            for(k=0; k<3; k++)
                r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
            R=(r[0]+r[1]+r[2])/3.;

            atot+=a/R;
            FluxLinkage+=2.*PI*R*a*(Aa/R);
        }
    }
    FluxLinkage*=( ((double) blocklist[lbl].Turns) / atot);

    return FluxLinkage;
}

CComplex FPProc::GetParallelLinkage(int numcirc)
{
    // routine for deducing the flux linkage of a "parallel-connected"
    // "circuit" in the annoying special case in which the
    // current carried in the circuit is zero and the frequency is zero.
    // This routine takes care of the case in which the current is divvied
    // up based on the conductivity and size of the various regions

    int i,k;
    CComplex FluxLinkage;
    CComplex Aa,A[3],J[3],U[3];//,V[3];
    double a,atot,R,c;
    a=atot=R=c = 0;
    double r[3];

    U[0]=1;
    U[1]=1;
    U[2]=1;

    for(i=0,FluxLinkage=0,atot=0; i<(int)meshelem.size(); i++)
    {
        if(blocklist[meshelem[i].lbl].InCircuit==numcirc)
        {
            c=blockproplist[meshelem[i].blk].Cduct;
            GetJA(i,J,A);
            a=ElmArea(i)*LengthConv[LengthUnits]*LengthConv[LengthUnits];

            if(ProblemType==AXISYMMETRIC)
            {
                for(k=0; k<3; k++)
                    r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
                R=(r[0]+r[1]+r[2])/3.;
                Aa=(A[0]+A[1]+A[2])/3.;
            }

            if(ProblemType==PLANAR)
            {
                FluxLinkage+=PlnInt(a,A,U)*Depth*c;
                atot+=(a*c);
            }
            else
            {
                FluxLinkage+=2.*PI*R*c*(Aa/R);
                atot+=(a*c/R);
            }
        }
    }
    FluxLinkage/=atot;

    return FluxLinkage;
}

CComplex FPProc::GetParallelLinkageAlt(int numcirc)
{
    // routine for deducing the flux linkage of a "parallel-connected"
    // "circuit" in the annoying special case in which the
    // current carried in the circuit is zero and the frequency is zero.
    // This routine takes care of the "punt" case in which all regions
    // in the "circuit" have been assigned a zero conductivity.
    // In this case, an even current density is applied to all regions
    // that are marked with the circuit (for both axi and planar cases).

    int i,k;
    CComplex FluxLinkage;
    CComplex Aa,A[3],J[3],U[3];//,V[3];
    double a,atot; //c,R;
    double r[3];

    U[0]=1;
    U[1]=1;
    U[2]=1;

    for(i=0,FluxLinkage=0,atot=0; i<(int)meshelem.size(); i++)
    {
        if(blocklist[meshelem[i].lbl].InCircuit==numcirc)
        {
//            c=blockproplist[meshelem[i].blk].Cduct;
            GetJA(i,J,A);
            a=ElmArea(i)*LengthConv[LengthUnits]*LengthConv[LengthUnits];
            atot+=a;

            if(ProblemType==AXISYMMETRIC)
            {
                for(k=0; k<3; k++)
                    r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
//                R=(r[0]+r[1]+r[2])/3.;
                Aa=(A[0]+A[1]+A[2])/3.;
            }

            if(ProblemType==PLANAR)    FluxLinkage+=PlnInt(a,A,U)*Depth;
            else FluxLinkage+=AxiInt(a,A,U,r);
        }
    }
    FluxLinkage/=atot;

    return FluxLinkage;
}

CComplex FPProc::GetVoltageDrop(int circnum)
{
    int i;
    CComplex Volts;

    Volts=0;

    // if the circuit is a "series" circuit...
    if(circproplist[circnum].CircType==1)
    {
        for(i=0; i<(int)blocklist.size(); i++)
        {
            if (blocklist[i].InCircuit==circnum)
            {
                // if this region is solid...
                if (blocklist[i].Case==0)
                {
                    // still need to multiply by turns, since it indicates
                    // the direction of the current;
                    if(ProblemType==AXISYMMETRIC)
                        Volts -= (2.*PI*blocklist[i].dVolts*blocklist[i].Turns);
                    else
                        Volts -= (Depth*blocklist[i].dVolts*blocklist[i].Turns);
                }
                // or if this region is stranded
                else Volts+=GetStrandedVoltageDrop(i);
            }
        }
    }
    else if(circproplist[circnum].CircType==0)
    {
        // root through the block labels until
        // we find one that gives us the voltage drop.
        bool flag=false;
        for(i=0; i<(int)blocklist.size(); i++)
        {
            if ((blocklist[i].InCircuit==circnum) && (blocklist[i].Case==0))
            {
                if(ProblemType==AXISYMMETRIC)
                    Volts -= (2.*PI*blocklist[i].dVolts);
                else
                    Volts -= (Depth*blocklist[i].dVolts);
                flag=true;
                i=blocklist.size();
            }
        }

        // but perhaps no voltage was found.  This could be a parallel circuit
        // in which the conductivity in every block in the circuit is equal
        // to zero.  We can still have flux linkage and voltage drop here.
        // have to root through things in a brute force way to get the voltage drop.
        if(flag==false)
        {
            int k;
            CComplex FluxLinkage;
            CComplex A[3],J[3],U[3];
            double a,atot;
            double r[3];

            U[0]=1;
            U[1]=1;
            U[2]=1;

            for(i=0,FluxLinkage=0,atot=0; i<(int)meshelem.size(); i++)
            {
                if(blocklist[meshelem[i].lbl].InCircuit==circnum)
                {
                    GetJA(i,J,A);
                    a=ElmArea(i)*LengthConv[LengthUnits]*LengthConv[LengthUnits];
                    atot+=a;

                    if(ProblemType==AXISYMMETRIC)
                    {
                        for(k=0; k<3; k++)
                            r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
                    }
                    if(ProblemType==PLANAR)
                        FluxLinkage+=PlnInt(a,A,U)*Depth;
                    else FluxLinkage+=AxiInt(a,A,U,r);
                }
            }
            Volts=(2.*PI*Frequency/atot)*FluxLinkage;
        }
    }

    return Volts;
}

CComplex FPProc::GetFluxLinkage(int circnum)
{
    int i,k;
    CComplex FluxLinkage;
    CComplex A[3],J[3];
    double a,r[3];

    // in the "normal" case, we can just use Integral of A.J
    // and divide through by i.conj to get the flux linkage.
    if((circproplist[circnum].Amps.re!=0) || (circproplist[circnum].Amps.im!=0))
    {
        for(i=0,FluxLinkage=0; i<(int)meshelem.size(); i++)
        {
            if(blocklist[meshelem[i].lbl].InCircuit==circnum)
            {
                GetJA(i,J,A);
                a=ElmArea(i)*LengthConv[LengthUnits]*LengthConv[LengthUnits];
                if(ProblemType==AXISYMMETRIC)
                {
                    for(k=0; k<3; k++)
                        r[k]=meshnode[meshelem[i].p[k]].x*LengthConv[LengthUnits];
                }

                // for a multiturn region, there can be some "local" flux linkage due to the complex-valued
                // part of the conductivity.
                if(Im(blocklist[meshelem[i].lbl].o)!=0)
                {
                    double u;
                    if(Frequency==0) u=Im(blocklist[meshelem[i].lbl].o);
                    else u=Im(1.e-6/blocklist[meshelem[i].lbl].o)/(2.*PI*Frequency);
                    for(k=0; k<3; k++) A[k]+=u*J[k];
                }

                for(k=0; k<3; k++) J[k]=J[k].Conj();
                if(ProblemType==PLANAR) FluxLinkage+=PlnInt(a,A,J)*Depth;
                else FluxLinkage+=AxiInt(a,A,J,r);


            }
        }

        FluxLinkage/=conj(circproplist[circnum].Amps);
    }
    else
    {
        // Rats!  The circuit of interest is not carrying any current.
        // However, the circuit can still have a non-zero flux linkage.
        // due to mutual inductance. Now, we have to go through
        // some annoying manipulations to back out the flux linkage.

        // For a non-zero frequency, things aren't too bad.  Any
        // voltage that we have must be solely due to flux linkage.
        // To get the flux linkage, just divide the voltage by the frequency.
        if (Frequency!=0)
            FluxLinkage=GetVoltageDrop(circnum)/(2.*PI*Frequency);

        // The zero frequency case is the interesting one, because
        // there are lots of annoying special cases.
        else
        {
            // First, deal with "series" circuits...
            if(circproplist[circnum].CircType==1)
            {
                for(i=0,FluxLinkage=0; i<(int)blocklist.size(); i++)
                {
                    if(blocklist[i].InCircuit==circnum)
                    {
                        if((blocklist[i].Case==1) || (ProblemType==PLANAR))
                            FluxLinkage+=GetStrandedLinkage(i);
                        // in solid, axisymmetric case, the current distribution
                        // isn't even.  We have to do some extra work to get
                        // the flux linkage in this case.
                        else FluxLinkage+=GetSolidAxisymmetricLinkage(i);
                    }
                }
            }
            else
            {
                bool flag;

                // first, see whether some of the blocks have nonzero conductivity;
                // flag=true means that at least one of the blocks has a
                // nonzero conductivity.
                for(i=0,flag=false; i<(int)blocklist.size(); i++)
                    if((blocklist[i].Case==0) &&
                            (blocklist[i].InCircuit==circnum)) flag=true;

                // if there is at least one nonzero conductivity block, we can use
                // the GetParallelLinkage routine, which is more or less driving
                // all the blocks with a ficticious voltage gradient.
                if (flag) FluxLinkage=GetParallelLinkage(i);
                // otherwise, treat the "punt" case, where every part of the
                // parallel "circuit" is just assumed to have the same applied
                // current density;
                else FluxLinkage=GetParallelLinkageAlt(i);
            }
        }
    }

    return FluxLinkage;
}

void FPProc::GetMagnetization(int n, CComplex &M1, CComplex &M2)
{
    // Puts the piece-wise constant magnetization for an element into
    // M1 and M2.  The magnetization could be useful for some kinds of
    // postprocessing, e.g. computation of field gradients by integrating
    // gradient contributions from each non-air element;

    CComplex mu1,mu2,Hc;
    CComplex b1,b2;

    b1=meshelem[n].B1;
    b2=meshelem[n].B2;
    Hc=0;
    mu1=0;
    mu2=0;

    if(Frequency==0)
    {
        GetMu(Re(b1),Re(b2),mu1.re,mu2.re,n);
        Hc=blockproplist[meshelem[n].blk].H_c*exp(I*meshelem[n].magdir*PI/180.);
    }
    else GetMu(b1,b2,mu1,mu2,n);

    M1 = b1*(mu1-1)/(mu1*muo) + Re(Hc);
    M2 = b2*(mu2-1)/(mu2*muo) + Im(Hc);
}

double FPProc::AECF(int k)
{
    // Computes the permeability correction factor for axisymmetric
    // external regions.  This is sort of a kludge, but it's the best
    // way I could fit it in.  The structure of the code wasn't really
    // designed to have a permeability that varies with position in a
    // continuous way.

    if (ProblemType!=AXISYMMETRIC) {
        return 1.; // no correction for planar problems
    }

    if (!blocklist[meshelem[k].lbl].IsExternal) {
        return 1; // only need to correct for external regions
    }

    double r=abs(meshelem[k].ctr-I*extZo);
    return (r*r*extRi)/(extRo*extRo*extRo); // permeability gets divided by this factor;
}

// versions of GetMu that sort out whether or not the AECF should be applied,
// as well as the corrections required for wound regions.
void FPProc::GetMu(CComplex b1, CComplex b2,CComplex &mu1, CComplex &mu2, int i)
{
    if(blockproplist[meshelem[i].blk].LamType>2) // is a region subject to prox effects
    {
        mu1=blocklist[meshelem[i].lbl].mu;
        mu2=mu1;
    }
    else blockproplist[meshelem[i].blk].GetMu(b1,b2,mu1,mu2);

    double aecf=AECF(i);
    mu1/=aecf;
    mu2/=aecf;
}

void FPProc::GetMu(double b1, double b2, double &mu1, double &mu2, int i)
{
    blockproplist[meshelem[i].blk].GetMu(b1,b2,mu1,mu2);
    double aecf=AECF(i);
    mu1/=aecf;
    mu2/=aecf;
}

void FPProc::GetH(double b1, double b2, double &h1, double &h2, int k)
{
    double mu1,mu2;
    CComplex Hc;

    GetMu(b1,b2,mu1,mu2,k);
    h1 = b1/(mu1*muo);
    h2 = b2/(mu2*muo);
    if ((d_ShiftH) && (blockproplist[meshelem[k].blk].H_c!=0))
    {
        Hc = blockproplist[meshelem[k].blk].H_c*
             exp(I*PI*meshelem[k].magdir/180.);
        h1=h1-Re(Hc);
        h2=h2-Im(Hc);
    }
}

void FPProc::GetH(CComplex b1, CComplex b2, CComplex &h1, CComplex &h2, int k)
{
    CComplex mu1,mu2;

    GetMu(b1,b2,mu1,mu2,k);
    h1 = b1/(mu1*muo);
    h2 = b2/(mu2*muo);
}

void FPProc::FindBoundaryEdges()
{
    int i, j;
    static int plus1mod3[3] = {1, 2, 0};
    static int minus1mod3[3] = {2, 0, 1};

    // Init all elements' neigh to be unfinished.
    for(i = 0; i < (int)meshelem.size(); i ++)
    {
        for(j = 0; j < 3; j ++)
            meshelem[i].n[j] = 0;
    }

    int orgi, desti;
    int ei, ni;
    bool done;

    // Loop all elements, to find and set there neighs.
    for(i = 0; i < (int)meshelem.size(); i ++)
    {
        for(j = 0; j < 3; j ++)
        {
            if(meshelem[i].n[j] == 0)
            {
                // Get this edge's org and dest node index,
                orgi = meshelem[i].p[plus1mod3[j]];
                desti = meshelem[i].p[minus1mod3[j]];
                done = false;
                // Find this edge's neigh from the org node's list
                for(ni = 0; ni < NumList[orgi]; ni ++)
                {
                    // Find a Element around org node contained dest node of this edge.
                    ei = ConList[orgi][ni];
                    if(ei == i) continue; // Skip myself.
                    // Check this Element's 3 vert to see if there exist dest node.
                    if(meshelem[ei].p[0] == desti)
                    {
                        done = true;
                        break;
                    }
                    else if(meshelem[ei].p[1] == desti)
                    {
                        done = true;
                        break;
                    }
                    else if(meshelem[ei].p[2] == desti)
                    {
                        done = true;
                        break;
                    }
                }
                if(!done)
                {
                    // This edge must be a Boundary Edge.
                    meshelem[i].n[j] = 1;
                }
            } // Finish One Edge
        } // End of One Element Loop
    } // End of Main Loop

}

