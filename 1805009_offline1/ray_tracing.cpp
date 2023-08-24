//implement ray tracing
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "bitmap_image.hpp"
using namespace std;
class Vector
{
    public:
    double x,y,z;
    Vector(double x,double y,double z)
    {
        this->x=x;
        this->y=y;
        this->z=z;
    }
    Vector()
    {
        this->x=0;
        this->y=0;
        this->z=0;
    }
    Vector (const Vector &p)
    {
        this->x=p.x;
        this->y=p.y;
        this->z=p.z;
    }
    Vector operator+(Vector p)
    {
        Vector temp;
        temp.x=this->x+p.x;
        temp.y=this->y+p.y;
        temp.z=this->z+p.z;
        return temp;
    }
    Vector operator-(Vector p)
    {
        Vector temp;
        temp.x=this->x-p.x;
        temp.y=this->y-p.y;
        temp.z=this->z-p.z;
        return temp;
    }
    Vector operator*(double d)
    {
        Vector temp;
        temp.x=this->x*d;
        temp.y=this->y*d;
        temp.z=this->z*d;
        return temp;
    }
    Vector operator/(double d)
    {
        Vector temp;
        temp.x=this->x/d;
        temp.y=this->y/d;
        temp.z=this->z/d;
        return temp;
    }
    double dot(Vector p)
    {
        return this->x*p.x+this->y*p.y+this->z*p.z;
    }
    Vector operator*(Vector p)
    {
        Vector temp;
        temp.x=this->y*p.z-this->z*p.y;
        temp.y=this->z*p.x-this->x*p.z;
        temp.z=this->x*p.y-this->y*p.x;
        return temp;
    }
    double length()
    {
        return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
    }
    Vector normalize()
    {
        Vector temp;
        double l=this->length();
        temp.x=this->x/l;
        temp.y=this->y/l;
        temp.z=this->z/l;
        return temp;
    }

};
class point
{
public:
    double x, y, z;
    point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;

    }
    point (const point &p)
    {
        this->x = p.x;
        this->y = p.y;
        this->z = p.z;
    }
    double length()
    {
        return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
    }
    point normalize()
    {
        point temp(0,0,0);
        double l = this->length();
        temp.x = this->x / l;
        temp.y = this->y / l;
        temp.z = this->z / l;
        return temp;
    }
    point operator+(Vector v)
    {
        point temp(0,0,0);
        temp.x = this->x + v.x;
        temp.y = this->y + v.y;
        temp.z = this->z + v.z;
        return temp;
    }
    point operator-(Vector v)
    {
        point temp(0,0,0);
        temp.x = this->x - v.x;
        temp.y = this->y - v.y;
        temp.z = this->z - v.z;
        return temp;
    }
    Vector operator-(point p)
    {
        Vector temp(0,0,0);
        temp.x = this->x - p.x;
        temp.y = this->y - p.y;
        temp.z = this->z - p.z;
        return temp;
    }
};

class Plane
{
public:
    point P;
    Vector N;
    Plane()
    {
        this->P = point(0, 0, 0);
        this->N = Vector(0, 0, 0);
    }
    Plane(point P, Vector N)
    {
        this->P = P;
        this->N = N;
    }

};
class cube
{
public:
    point bottom_lower_left;
    double side;
    double red,green,blue;
    cube(point bottom_lower_left, double side)
    {
        this->bottom_lower_left = bottom_lower_left;
        this->side = side;
    }
    void set_color(double red,double green,double blue)
    {
        this->red=red;
        this->green=green;
        this->blue=blue;
    }
};
class triangle
{
public:
    point A, B, C;
    triangle(point A, point B, point C)
    {
        this->A = A;
        this->B = B;
        this->C = C;
    }
};
class pyramid
{
    public:
    point bottom_center;
    double side;
    double height;
    double red,green,blue;
    pyramid(point bottom_center,double side,double height)
    {
        this->bottom_center=bottom_center;
        this->side=side;
        this->height=height;
    }
    void set_color(double red,double green,double blue)
    {
        this->red=red;
        this->green=green;
        this->blue=blue;
    }
};
double find_intersection_point_with_plane(point R_o, Vector ray, Plane plane)
{
    double t = (plane.P - R_o).dot(plane.N) / ray.dot(plane.N);
    return t;
}
double find_intersection_point_with_cube(point R_o, Vector ray, cube c)
{
    Plane planes[6];
    double t;
    planes[0] = Plane(c.bottom_lower_left, Vector(0, 0, 1));
    t=find_intersection_point_with_plane(R_o,ray,planes[0]);
    double min_t=-10;
    if(t>0)
    {
        point intersection_point=R_o+ray*t;
        if(intersection_point.x>=c.bottom_lower_left.x && intersection_point.x<=c.bottom_lower_left.x+c.side && intersection_point.y>=c.bottom_lower_left.y && intersection_point.y<=c.bottom_lower_left.y+c.side)
        {
           min_t=t; 
        }
    }
    planes[1] = Plane(c.bottom_lower_left, Vector(0, 1, 0));
    t=find_intersection_point_with_plane(R_o,ray,planes[1]);
    if(t>0)
    {
        point intersection_point=R_o+ray*t;
        if(intersection_point.x>=c.bottom_lower_left.x && intersection_point.x<=c.bottom_lower_left.x+c.side && intersection_point.z>=c.bottom_lower_left.z && intersection_point.z<=c.bottom_lower_left.z+c.side)
        {
           if (min_t==-10 || t<min_t)
           min_t=t; 
        }
    }
    planes[2] = Plane(c.bottom_lower_left, Vector(1, 0, 0));
    t=find_intersection_point_with_plane(R_o,ray,planes[2]);
    if(t>0)
    {
        point intersection_point=R_o+ray*t;
        if(intersection_point.y>=c.bottom_lower_left.y && intersection_point.y<=c.bottom_lower_left.y+c.side && intersection_point.z>=c.bottom_lower_left.z && intersection_point.z<=c.bottom_lower_left.z+c.side)
        {
           if (min_t==-10 || t<min_t)
           min_t=t; 
        }
    }
    planes[3] = Plane(c.bottom_lower_left + Vector(0, 0, c.side), Vector(0, 0, -1));
    t=find_intersection_point_with_plane(R_o,ray,planes[3]);
    if(t>0)
    {
        point intersection_point=R_o+ray*t;
        if(intersection_point.x>=c.bottom_lower_left.x && intersection_point.x<=c.bottom_lower_left.x+c.side && intersection_point.y>=c.bottom_lower_left.y && intersection_point.y<=c.bottom_lower_left.y+c.side)
        {
           if (min_t==-10 || t<min_t)
           min_t=t; 
        }
    }
    planes[4] = Plane(c.bottom_lower_left + Vector(0, c.side, 0), Vector(0, -1, 0));
    t=find_intersection_point_with_plane(R_o,ray,planes[4]);
    if(t>0)
    {
        point intersection_point=R_o+ray*t;
        if(intersection_point.x>=c.bottom_lower_left.x && intersection_point.x<=c.bottom_lower_left.x+c.side && intersection_point.z>=c.bottom_lower_left.z && intersection_point.z<=c.bottom_lower_left.z+c.side)
        {
           if (min_t==-10 || t<min_t)
           min_t=t; 
        }
    }
    planes[5] = Plane(c.bottom_lower_left + Vector(c.side, 0, 0), Vector(-1, 0, 0));
    t=find_intersection_point_with_plane(R_o,ray,planes[5]);
    if(t>0)
    {
        point intersection_point=R_o+ray*t;
        if(intersection_point.y>=c.bottom_lower_left.y && intersection_point.y<=c.bottom_lower_left.y+c.side && intersection_point.z>=c.bottom_lower_left.z && intersection_point.z<=c.bottom_lower_left.z+c.side)
        {
           if (min_t==-10 || t<min_t)
           min_t=t; 
        }
    }
    return min_t;
    
}
class Sphere
{
public:
    point center;
    double radius;
    double red,green,blue;
    Sphere(point center, double radius)
    {
        this->center = center;
        this->radius = radius;
    }
    void set_color(double red,double green,double blue)
    {
        this->red=red;
        this->green=green;
        this->blue=blue;
    }
};
class checkerboard
{   
    public:
    double side;
    checkerboard(double side)
    {
        this->side=side;
    }
};

class Matrix
{
public:
    double mat[3][3];
    Matrix()
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                if (i == j)
                    mat[i][j] = 1;
                else
                    mat[i][j] = 0;
            }
    }
    Matrix(double mat[3][3])
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                this->mat[i][j] = mat[i][j];
            }
    }
    Matrix operator*(const Matrix &other) const
    {
        Matrix ret;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                ret.mat[i][j] = 0;
                for (int k = 0; k < 3; k++)
                    ret.mat[i][j] += mat[i][k] * other.mat[k][j];
            }
        return ret;
    }
    point operator*(const point &p) const
    {
        point ret;
        ret.x = mat[0][0] * p.x + mat[0][1] * p.y + mat[0][2] * p.z + mat[0][3];
        ret.y = mat[1][0] * p.x + mat[1][1] * p.y + mat[1][2] * p.z + mat[1][3];
        ret.z = mat[2][0] * p.x + mat[2][1] * p.y + mat[2][2] * p.z + mat[2][3];
        return ret;
    }
    double determinant()
    {
        double det = 0;
        for (int i = 0; i < 3; i++)
        {
            det += mat[0][i] * (mat[1][(i + 1) % 3] * mat[2][(i + 2) % 3] - mat[1][(i + 2) % 3] * mat[2][(i + 1) % 3]);
        }
        return det;
    }
};
double find_intersection_point_with_triangle(point R_o, Vector ray, triangle t)
{
    Matrix A;
    A.mat[0][0] = t.A.x - t.B.x;
    A.mat[0][1] = t.A.x - t.C.x;
    A.mat[0][2] = ray.x;
    A.mat[1][0] = t.A.y - t.B.y;
    A.mat[1][1] = t.A.y - t.C.y;
    A.mat[1][2] = ray.y;
    A.mat[2][0] = t.A.z - t.B.z;
    A.mat[2][1] = t.A.z - t.C.z;
    A.mat[2][2] = ray.z;
    Matrix beta;
    beta.mat[0][0] = t.A.x - R_o.x;
    beta.mat[0][1] = t.A.x - t.C.x;
    beta.mat[0][2] = ray.x;
    beta.mat[1][0] = t.A.y - R_o.y;
    beta.mat[1][1] = t.A.y - t.C.y;
    beta.mat[1][2] = ray.y;
    beta.mat[2][0] = t.A.z - R_o.z;
    beta.mat[2][1] = t.A.z - t.C.z;
    beta.mat[2][2] = ray.z;
    Matrix gamma;
    gamma.mat[0][0] = t.A.x - t.B.x;
    gamma.mat[0][1] = t.A.x - R_o.x;
    gamma.mat[0][2] = ray.x;
    gamma.mat[1][0] = t.A.y - t.B.y;
    gamma.mat[1][1] = t.A.y - R_o.y;
    gamma.mat[1][2] = ray.y;
    gamma.mat[2][0] = t.A.z - t.B.z;
    gamma.mat[2][1] = t.A.z - R_o.z;
    gamma.mat[2][2] = ray.z;
    Matrix t_mat;
    t_mat.mat[0][0] = t.A.x - t.B.x;
    t_mat.mat[0][1] = t.A.x - t.C.x;
    t_mat.mat[0][2] = t.A.x - R_o.x;
    t_mat.mat[1][0] = t.A.y - t.B.y;
    t_mat.mat[1][1] = t.A.y - t.C.y;
    t_mat.mat[1][2] = t.A.y - R_o.y;
    t_mat.mat[2][0] = t.A.z - t.B.z;
    t_mat.mat[2][1] = t.A.z - t.C.z;
    t_mat.mat[2][2] = t.A.z - R_o.z;
    double det_A = A.determinant();
    double beta_val = beta.determinant() / det_A;
    double gamma_val = gamma.determinant() / det_A;
    double t_val = t_mat.determinant() / det_A;
    if (beta_val > 0 && gamma_val > 0 && beta_val + gamma_val < 1 && t_val > 0)
        return t_val;
    else
        return -1;
    
}
double find_intersection_point_with_pyramid(point R_o, Vector ray, pyramid p)
{
    triangle t1(p.bottom_center + Vector(-p.side / 2, 0, -p.side / 2), p.bottom_center + Vector(p.side / 2, 0, -p.side / 2), p.bottom_center + Vector(0, p.height, 0));
    triangle t2(p.bottom_center + Vector(p.side / 2, 0, -p.side / 2), p.bottom_center + Vector(p.side / 2, 0, p.side / 2), p.bottom_center + Vector(0, p.height, 0));
    triangle t3(p.bottom_center + Vector(p.side / 2, 0, p.side / 2), p.bottom_center + Vector(-p.side / 2, 0, p.side / 2), p.bottom_center + Vector(0, p.height, 0));
    triangle t4(p.bottom_center + Vector(-p.side / 2, 0, p.side / 2), p.bottom_center + Vector(-p.side / 2, 0, -p.side / 2), p.bottom_center + Vector(0, p.height, 0));
    double t1_val = find_intersection_point_with_triangle(R_o, ray, t1);
    double t2_val = find_intersection_point_with_triangle(R_o, ray, t2);
    double t3_val = find_intersection_point_with_triangle(R_o, ray, t3);
    double t4_val = find_intersection_point_with_triangle(R_o, ray, t4);
    double t5_val = find_intersection_point_with_plane(R_o, ray, Plane(p.bottom_center, Vector(0, 1, 0)));
    double min_t = -10;
    if (t1_val > 0)
    {
        min_t = t1_val;
    }
    if (t2_val > 0)
    {
        if (min_t == -10 || t2_val < min_t)
            min_t = t2_val;
    }
    if (t3_val > 0)
    {
        if (min_t == -10 || t3_val < min_t)
            min_t = t3_val;
    }
    if (t4_val > 0)
    {
        if (min_t == -10 || t4_val < min_t)
            min_t = t4_val;
    }
    if (t5_val > 0)
    {
        point temp=R_o+ray*t5_val;
        if(temp.x>=p.bottom_center.x-p.side/2 && temp.x<=p.bottom_center.x+p.side/2 && temp.z>=p.bottom_center.z-p.side/2 && temp.z<=p.bottom_center.z+p.side/2)
            if (min_t == -10 || t5_val < min_t)
                min_t = t5_val;
    }
    return min_t;
}
double find_intersection_point_with_sphere(point R_o, Vector ray, Sphere s)
{
    double a = (R_o-s.center).dot(R_o-s.center);
    double tp=-ray.dot(R_o-s.center);
    if (a>0 && tp<0)
    {
        return -10;
    }
    double d2=a-tp*tp;
    if(d2>s.radius*s.radius)
    {
        return -10;
    }
    double t=sqrt(s.radius*s.radius-d2);

    if (a>s.radius*s.radius)
    {
        return tp-t;
    }
    else if (a<s.radius*s.radius)
    {
        return tp+t;
    }
    else
    {
        return tp+0.001;
    }

}
double find_intersection_point_with_checkerboard(point R_o, Vector ray, checkerboard c)
{
    Plane plane(point(0,0,0),Vector(0,1,0));
    double t=find_intersection_point_with_plane(R_o,ray,plane);
    if(t>0)
    {
            return t;
    }
    return -10;
}
void create_pointBuffer(double near_plane_d,double far_plane_d,double fovy,double aspect_ratio,Vector l,Vector r,Vector u,point camera,int image_width,int image_height)
{
    // double near_plane_d=1;
    // double far_plane_d=100;
    fovy=fovy*acos(-1)/180;
    // double aspect_ratio=1;
    double near_plane_h=2*near_plane_d*tan(fovy/2);
    double fovx=aspect_ratio*fovy;
    double near_plane_w=2*near_plane_d*tan(fovx/2);
    // Vector l=Vector(0,0,-1);
    // Vector r=Vector(1,0,0);
    // Vector u=Vector(0,1,0);
    // point camera=point(0,1,5);
    // int image_width=512;
    // int image_height=512;
    point midpoint=camera+l*near_plane_d;
    double w_per_pixel=near_plane_w/image_width;
    double h_per_pixel=near_plane_h/image_height;
    point topleft=midpoint+u*(near_plane_h/2-h_per_pixel/2)-r*(near_plane_w/2-w_per_pixel/2);
    // point topright=midpoint+u*near_plane_h/2+l*near_plane_w/2;
    // point bottomleft=midpoint-u*near_plane_h/2-l*near_plane_w/2;
    // point bottomright=midpoint-u*near_plane_h/2+l*near_plane_w/2;
    point **pointBuffer=new point*[image_width];
    for(int i=0;i<image_width;i++)
    {
        pointBuffer[i]=new point[image_height];
    }
    for(int i=0;i<image_height;i++)
    {
        for(int j=0;j<image_width;j++)
        {
            pointBuffer[j][i]=topleft+u*(-h_per_pixel*i)+r*(w_per_pixel*j);
        }
    }
    checkerboard cb=checkerboard(1);
    cube c=cube(point(0,1,1),1);
    Sphere s=Sphere(point(-1,1,0),1);
    pyramid p=pyramid(point(1,1,0),1,1);
    c.set_color(0,1,1);
    s.set_color(1,0,0);
    p.set_color(1,1,0.5);
    bitmap_image* image;
    image = new bitmap_image(image_width, image_height);
    for (int i = 0; i < image_height; i++)
        for (int j = 0; j < image_width; j++)
            image->set_pixel(j, i, 0, 0, 0);
    for(int i=0;i<image_height;i++)
    {
        for(int j=0;j<image_width;j++)
        {
            Vector ray=pointBuffer[j][i]-camera;
            ray=ray.normalize();
            double t1=find_intersection_point_with_cube(pointBuffer[j][i],ray,c);
            double t2=find_intersection_point_with_sphere(pointBuffer[j][i],ray,s);
            double t3=find_intersection_point_with_pyramid(pointBuffer[j][i],ray,p);
            double t4=find_intersection_point_with_checkerboard(pointBuffer[j][i],ray,cb);
            double min_t=-10;
            if(t1>0)
            {
                min_t=t1;
                image->set_pixel(j, i, c.red*255, c.green*255, c.blue*255);
            }
            if(t2>0)
            {
                if(min_t==-10 || t2<min_t)
                {min_t=t2;
                image->set_pixel(j, i, s.red*255, s.green*255, s.blue*255);
                }
            }
            if(t3>0)
            {
                if(min_t==-10 || t3<min_t)
                {min_t=t3;
                image->set_pixel(j, i, p.red*255, p.green*255, p.blue*255);
                }
            }
            if(t4>0)
            {
                if(min_t==-10 || t4<min_t)
                {min_t=t4;
                point temp=pointBuffer[j][i]+ray*min_t;
                int x=abs(temp.x)/cb.side;
                int z=abs(temp.z)/cb.side;
                if((x+z)%2==0)
                {
                if (temp.x*temp.z>0)
                image->set_pixel(j, i, 0, 0, 0);
                else
                image->set_pixel(j, i, 255, 255, 255);
                }
                else
                {
                if (temp.x*temp.z>0)
                image->set_pixel(j, i, 255, 255, 255);
                else
                image->set_pixel(j, i, 0, 0, 0);
                }
                }
            }
        }
    }
    image->save_image("output.bmp");


    // int c=0;
    // for(int i=0;i<image_height;i++)
    // {
    //     for(int j=0;j<image_width;j++)
    //     {
    //         cout<<pointBuffer[j][i].x<<" "<<pointBuffer[j][i].y<<" "<<pointBuffer[j][i].z<<",";
    //     }
    //     cout<<endl;
    // }


}

// int main()
// {
//     create_pointBuffer();
//     return 0;
// }
