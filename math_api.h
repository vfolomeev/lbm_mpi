#pragma once
#ifndef MATH_API_H
#define MATH_API_H

#endif // MATH_API_H
#include <math.h>
//#include "QTextStream"
//#include "QDebug"
#include <iostream>
inline int ind(int i,int j,int k){
    return (i+1)*9+(j+1)*3+(k+1);
}

class Vector3d{
public:
    double x;
    double y;
    double z;
    Vector3d();
    Vector3d(double a);
    Vector3d(double a,double b,double c);

    double mag();
    Vector3d  operator % (Vector3d a);

    friend Vector3d operator+(Vector3d a,Vector3d b);

    friend Vector3d operator*(double a,Vector3d b);

    friend Vector3d operator*(Vector3d a,double b);

    friend double operator*(Vector3d a,Vector3d b);

    friend Vector3d operator>=(Vector3d a,Vector3d b);

    friend Vector3d operator<=(Vector3d a,Vector3d b);

};
//QDebug operator<< (QDebug out,const Vector3d &r);
//QTextStream &operator <<(QTextStream &s, const Vector3d &r);
std::ostream& operator <<(std::ostream& cout_,const  Vector3d & r);
