/*************************************************************************
Geometry Functions for Structural Bioinformatics
This is a reimplementation of Dr. Yang Cao's GeometryTools.hpp
*************************************************************************/

#ifndef Geometryhpp
#define Geometryhpp

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib> 

using namespace std;

const double PI=3.14159265358979323846;
const double Extra=1.0e-4;
const double UpMax=1.0e+10;

template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
}

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
}

//Functions using vectors
inline void crossproduct(const double *c1, const double *c2, double *cc)
{
    cc[0] = c1[1] * c2[2] - c1[2] * c2[1];
    cc[1] = c1[2] * c2[0] - c1[0] * c2[2];
    cc[2] = c1[0] * c2[1] - c2[0] * c1[1];
}

inline double innerproduct(const double *c1, const double *c2)
{
    return c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
}

inline void subtract(const double *c1, const double *c2, double *cc)
{
    cc[0] = c1[0] - c2[0];
    cc[1] = c1[1] - c2[1];
    cc[2] = c1[2] - c2[2];
}

inline double deg2rad(double deg)
{
    return deg*PI/180;
}

inline double rad2deg(double rad)
{
    return rad*180/PI;
}

inline void vectorsum(const double *c1, const double *c2, double *cc)
{
    cc[0] = c1[0] + c2[0];
    cc[1] = c1[1] + c2[1];
    cc[2] = c1[2] + c2[2];
}

inline bool norm(const double *c, double *cc)
{
    double len = c[0]*c[0] + c[1]*c[1] +c[2]*c[2];
    if(len<Extra) 
    {
        cc[0] = 1; cc[1] = 0; cc[2] = 0;
        cerr<<"Error in norm(). length~=0\n";
        return false;
    }
    len = 1/sqrt(len);
    cc[0] = c[0]*len;
    cc[1] = c[1]*len;
    cc[2] = c[2]*len;
    return true;
}

/************************Rotation Matrix****************************
* Give the Rotation Axis and the Rotation Angle
* Generate the Rotation Matrix
*************************End***************************************/
inline bool RotationMatrixB(const double *axis, double angle, double rotmtx[3][3])
{
    rotmtx[0][0]=1; rotmtx[0][1]=0; rotmtx[0][2]=0;
    rotmtx[1][0]=0; rotmtx[1][1]=1; rotmtx[1][2]=0;
    rotmtx[2][0]=0; rotmtx[2][1]=0; rotmtx[2][2]=1;

    double ouc[3];
    if(!norm(axis, ouc)) 
    {
        cerr<<"Error 2 in RotationMatrixB()"<<endl;
        return false;
    }
    double c=cos(angle);
    double s=sin(angle);
    double t=1-c;
    rotmtx[0][0]=t*ouc[0]*ouc[0]+c;
    rotmtx[0][1]=t*ouc[0]*ouc[1]+s*ouc[2];
    rotmtx[0][2]=t*ouc[0]*ouc[2]-s*ouc[1];
   
    rotmtx[1][0]=t*ouc[1]*ouc[0]-s*ouc[2];
    rotmtx[1][1]=t*ouc[1]*ouc[1]+c;
    rotmtx[1][2]=t*ouc[1]*ouc[2]+s*ouc[0];
   
    rotmtx[2][0]=t*ouc[2]*ouc[0]+s*ouc[1];
    rotmtx[2][1]=t*ouc[2]*ouc[1]-s*ouc[0];
    rotmtx[2][2]=t*ouc[2]*ouc[2]+c;
    return true;
}

inline bool MatrixTimesTransVector(double mtx[3][3],
    const double *tvt, double *vct)
{
    int i=0, j=0;
    vct[0]=vct[1]=vct[2]=0;
    for(i=0; i<3; ++i) for(j=0; j<3; ++j) vct[i]+=mtx[i][j]*tvt[j];
    return true;
}

/***************************CoordinateRotation********************
* For a PointA, Given the Coordinates of 2 points(axisA, axisB) 
* in line of Rotation Axis and the Rotation Angle 
* Generate the coordinate of PointB after operation
*****************************END**********************************/
bool CoordinateRotation(const double *pointA, const double *axisA,
    const double *axisB, const double angle, double *pointB)
{
    pointB[0]=pointA[0];
    pointB[1]=pointA[1];
    pointB[2]=pointA[2];

    double axis[3];
    axis[0]=axisB[0]-axisA[0]; 
    axis[1]=axisB[1]-axisA[1];
    axis[2]=axisB[2]-axisA[2]; 
   
    double rotmtx[3][3];
    if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
    double point_A[3];
    point_A[0]=pointA[0]-axisA[0];
    point_A[1]=pointA[1]-axisA[1];
    point_A[2]=pointA[2]-axisA[2];
   
    double point_B[3];
    MatrixTimesTransVector(rotmtx, point_A, point_B);
    pointB[0]=point_B[0]+axisA[0];
    pointB[1]=point_B[1]+axisA[1];
    pointB[2]=point_B[2]+axisA[2];
    return true;
}

/************************Points to Dihedral************************/
//If use inline here there will be an error when compiling.
//Return the angle of c1-c2-c3-c4. Unit: radian
//points c1-c2-c3-c4 should be in two different planes.
//or there may be faults.
//by Yang Cao
//
//Change by chengxin: if they are on the same plane, return -2*PI
inline double Points2Dihedral(const double*c1, const double*c2, const double*c3, const double*c4)
{
    double vector1[3];
    double vector2[3];
    double vector3[3];
    subtract(c1, c2, vector1);
    subtract(c2, c3, vector2);
    subtract(c3, c4, vector3);
   
    double v1[3];
    double v2[3];
    crossproduct(vector2, vector1, v1);
    crossproduct(vector3, vector2, v2);
   
    double v3[3];
    double v4[3];
    if(!norm (v1, v3)|| !norm (v2, v4))
    {
        cerr<<"Error in Points2Dihedral 2"<<endl;
        return -2*PI;
    }
   
    double dihedral = innerproduct(v3, v4);
   
    if(dihedral>1+Extra || dihedral<-1-Extra) cerr<<"Error, double Points2Dihedral()"<<endl;
    if (dihedral>1) dihedral=1;
    else if(dihedral<-1) dihedral=-1;
   
    double v5[3];
    crossproduct(v4, v3, v5);
    double direction = innerproduct(v5, vector2);
   
    if (direction>0) return  acos(dihedral);
    else             return -acos(dihedral);
}

inline double Points2Distance(const double *c1, const double *c2)
{
    double a=(c1[0]-c2[0]);
    double b=(c1[1]-c2[1]);
    double c=(c1[2]-c2[2]);
    return sqrt(a*a+b*b+c*c);
}

//Return the angle of c1-c2-c3. Unit: radian
//Angle <c1c2c3 
inline double Points2Angle(const double *c1, const double *c2, const double *c3)
{
    double a, b, c, tmp1, tmp2, tmp3, alpha;
   
    tmp1=c1[0]-c2[0];
    tmp2=c1[1]-c2[1];
    tmp3=c1[2]-c2[2];
    a=tmp1*tmp1+tmp2*tmp2+ tmp3*tmp3;
   
    tmp1=c2[0]-c3[0];
    tmp2=c2[1]-c3[1];
    tmp3=c2[2]-c3[2];
    b=tmp1*tmp1+tmp2*tmp2+ tmp3*tmp3;
   
    tmp1=c3[0]-c1[0];
    tmp2=c3[1]-c1[1];
    tmp3=c3[2]-c1[2];
    c=tmp1*tmp1+tmp2*tmp2+ tmp3*tmp3;
   
    alpha=(a+b-c)/(2*sqrt(a*b));
   
    if (alpha>1+Extra || alpha<-1-Extra) cerr<<"Error, double Points2Angle()\n";
    if      (alpha> 1) alpha=1;
    else if (alpha<-1) alpha=-1;
    return acos(alpha);
}

//Return the angle of <c1-c2,c3-c4>. Unit: radian
inline double Points4Angle(const double *c1, const double *c2,
    const double *c3, const double *c4)
{
    double a, b, a1, a2, a3, b1, b2, b3, alpha;
   
    a1=c2[0]-c1[0];
    a2=c2[1]-c1[1];
    a3=c2[2]-c1[2];
    a = a1*a1 + a2*a2 + a3*a3;
   
    b1=c4[0]-c3[0];
    b2=c4[1]-c3[1];
    b3=c4[2]-c3[2];
    b = b1*b1 + b2*b2 + b3*b3;
   
    alpha=(a1*b1+a2*b2+a3*b3)/sqrt(a*b);
   
    if(alpha>1+Extra || alpha<-1-Extra)
        cerr<<"Error, double Points2Angle()\n";
   
    if (alpha>1)      alpha=1;
    else if(alpha<-1) alpha=-1;
    return acos(alpha);
}

/***************************CoordinateRotation********************
* For a PointA, Given the Coordinates of 2 points(axisA, axisB) 
* in line of Rotation Axis and the Rotation Angle 
* Generate the coordinate of PointB after operation
*****************************END**********************************/
bool GroupRotation(const double *axisA, const double *axisB, const double angle,
    double **pointB, const size_t rstart, const size_t rend)
{
    double axis[3];
    axis[0]=axisB[0]-axisA[0]; 
    axis[1]=axisB[1]-axisA[1];
    axis[2]=axisB[2]-axisA[2]; 
   
    double rotmtx[3][3];
    if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
    double point_A[3];
    for (size_t i=rstart; i<rend; i++)
    {
        point_A[0]=pointB[i][0]-axisA[0];
        point_A[1]=pointB[i][1]-axisA[1];
        point_A[2]=pointB[i][2]-axisA[2];
      
        pointB[i][0]=axisA[0]+rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];
        pointB[i][1]=axisA[1]+rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];
        pointB[i][2]=axisA[2]+rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];
    }
   
    return true;
}

/******************************************************************************
Coordinate Translation for pointB[rstart:rend,:] alone trans
******************************************************************************/
void GroupTranslation(const double *trans, double ** pointB, 
    const size_t rstart, const size_t rend)
{
   for(size_t i=rstart; i<rend; ++i)
   {
        pointB[i][0]+=trans[0];
        pointB[i][1]+=trans[1];
        pointB[i][2]+=trans[2];
   }
}
/* return 0 if successful; return 1 if otherwise.
 * both ang and di are in the unit of degree, not rad */
int calculateCoordinates(double *D, 
    const double *refA, const double *refB, const double *refC, 
    const double l, const double ang, const double di)
{
    double CA[3];
    double CB[3];
    double CC[3];
    double tmp_array[3];
    subtract(refA,refB,CA);
    subtract(refB,refC,CB);
    subtract(refA,refC,CC);
    if (innerproduct(CA,CA)<Extra || innerproduct(CB,CB)<Extra || innerproduct(CC,CC)<Extra)
    {
        D[0]=refC[0];
        D[1]=refC[0];
        D[2]=refC[0];
        return 1;
    }

    /* Plane Parameters */
    double A=(CA[1]*CB[2]) - (CA[2]*CB[1]);
    double B=(CA[2]*CB[0]) - (CA[0]*CB[2]);
    double G=(CA[0]*CB[1]) - (CA[1]*CB[0]);

    /* Dot Product Constant */
    double F=sqrt(CB[0]*CB[0] + CB[1]*CB[1] + CB[2]*CB[2]
        ) * l * cos(deg2rad(ang));

    double cons=B*CB[2]-CB[1]*G;
    cons=sqrt(cons*cons*(-(F*F)*(A*A+B*B+G*G)+(B*B*(CB[0]*CB[0]+CB[2]*CB[2]
        ) + A*A*(CB[1]*CB[1]+CB[2]*CB[2])- (2*A*CB[0]*CB[2]*G) + (
        CB[0]*CB[0]+ CB[1]*CB[1])*G*G- (2*B*CB[1])*(A*CB[0]+CB[2]*G))*l*l));
    double denom=B*B*(CB[0]*CB[0]+CB[2]*CB[2]) + A*A*(CB[1]*CB[1]+CB[2]*CB[2]
        ) - (2*A*CB[0]*CB[2]*G) + G*G*(CB[0]*CB[0]+CB[1]*CB[1]
        ) - (2*B*CB[1])*(A*CB[0]+CB[2]*G);

    double X=((B*B*CB[0]*F)-(A*B*CB[1]*F)+(F*G)*(-A*CB[2]+CB[0]*G)+cons
        )/denom;

    double Y,Z;
    if ((B==0 or CB[2]==0) && (CB[1]==0 or G==0))
    {
        double const1=sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(l-X)*(l+X)));
        Y= ((-A*B*X)+const1)/(B*B+G*G);
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G));
    }
    else
    {
        Y= ((A*A*CB[1]*F)*(B*CB[2]-CB[1]*G)+ G*(-F*pow(B*CB[2]-CB[1]*G,2
            ) + CB[0]*cons) - A*( B*B*CB[0]*CB[2]*F- B*CB[0]*CB[1]*F*G + \
            CB[2]*cons)) / ((B*CB[2]-CB[1]*G)*denom);
        Z= ((A*A*CB[2]*F)*(B*CB[2]-CB[1]*G) + (B*F)*pow(B*CB[2]-CB[1]*G,2
            ) + (A*CB[0]*F*G)*(-B*CB[2]+CB[1]*G) - B*CB[0]*cons + \
            A*CB[1]*cons) / ((B*CB[2]-CB[1]*G)*denom);
    }

    
    /* GET THE NEW VECTOR from the orgin */
    tmp_array[0]=X;
    tmp_array[1]=Y;
    tmp_array[2]=Z;
    vectorsum(tmp_array, refC, D);

    double angle=di-rad2deg(Points2Dihedral(refA, refB, refC, D));
    CoordinateRotation(D,refC,refB, angle, tmp_array);
    D[0]=tmp_array[0];
    D[1]=tmp_array[1];
    D[2]=tmp_array[2];
    return 0;
}

#endif
