
#include "math_api.h"


//Vector3d





Vector3d operator>=(Vector3d a,Vector3d b){
    return Vector3d(a.x>=b.x,a.y>=b.y,a.z>=b.z);
}

Vector3d operator<=(Vector3d a,Vector3d b){
    return Vector3d(a.x<=b.x,a.y<=b.y,a.z<=b.z);
}


Vector3d::Vector3d(){
    x=0;y=0;z=0;
}

Vector3d::Vector3d(double a){
    x=a;y=a;z=a;
}

Vector3d::Vector3d(double a,double b,double c){
    x=a;y=b;z=c;
}

double Vector3d::mag(){
    return sqrt(x*x+y*y+z*z);
}

std::ostream& operator <<(std::ostream& cout_,const Vector3d & r){
    cout_<<" "<<r.x<<" "<<r.y<<"  "<<r.z;
    return cout_;
}





