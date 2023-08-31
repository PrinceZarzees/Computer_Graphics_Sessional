#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "bitmap_image.hpp"
using namespace std;

bool texture = false;
bitmap_image *texture_w;
bitmap_image *texture_b;
double ***texture_w_array;
double ***texture_b_array;

class Vector
{
public:
    double x, y, z;
    Vector(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vector()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    Vector(const Vector &p)
    {
        this->x = p.x;
        this->y = p.y;
        this->z = p.z;
    }
    Vector operator+(Vector p)
    {
        Vector temp;
        temp.x = this->x + p.x;
        temp.y = this->y + p.y;
        temp.z = this->z + p.z;
        return temp;
    }
    Vector operator-(Vector p)
    {
        Vector temp;
        temp.x = this->x - p.x;
        temp.y = this->y - p.y;
        temp.z = this->z - p.z;
        return temp;
    }
    Vector operator*(double d)
    {
        Vector temp;
        temp.x = this->x * d;
        temp.y = this->y * d;
        temp.z = this->z * d;
        return temp;
    }
    Vector operator/(double d)
    {
        Vector temp;
        temp.x = this->x / d;
        temp.y = this->y / d;
        temp.z = this->z / d;
        return temp;
    }
    double dot(Vector p)
    {
        return this->x * p.x + this->y * p.y + this->z * p.z;
    }
    Vector operator*(Vector p)
    {
        Vector temp;
        temp.x = this->y * p.z - this->z * p.y;
        temp.y = this->z * p.x - this->x * p.z;
        temp.z = this->x * p.y - this->y * p.x;
        return temp;
    }
    double length()
    {
        return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
    }
    Vector normalize()
    {
        Vector temp;
        double l = this->length();
        temp.x = this->x / l;
        temp.y = this->y / l;
        temp.z = this->z / l;
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
    point(const point &p)
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
        point temp(0, 0, 0);
        double l = this->length();
        temp.x = this->x / l;
        temp.y = this->y / l;
        temp.z = this->z / l;
        return temp;
    }
    point operator+(Vector v)
    {
        point temp(0, 0, 0);
        temp.x = this->x + v.x;
        temp.y = this->y + v.y;
        temp.z = this->z + v.z;
        return temp;
    }
    point operator-(Vector v)
    {
        point temp(0, 0, 0);
        temp.x = this->x - v.x;
        temp.y = this->y - v.y;
        temp.z = this->z - v.z;
        return temp;
    }
    Vector operator-(point p)
    {
        Vector temp(0, 0, 0);
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
    Plane(const Plane &p)
    {
        this->P = p.P;
        this->N = p.N;
    }
    Plane(point P, Vector N)
    {
        this->P = P;
        this->N = N;
    }
};

class intersection_point_property
{
public:
    point P;
    Vector N;
    double t;
    double red, green, blue;
    double ambient, diffuse, specular, reflection, shininess;
    intersection_point_property(const intersection_point_property &p)
    {
        this->P = p.P;
        this->N = p.N;
        this->t = p.t;
        this->red = p.red;
        this->green = p.green;
        this->blue = p.blue;
        this->ambient = p.ambient;
        this->diffuse = p.diffuse;
        this->specular = p.specular;
        this->reflection = p.reflection;
        this->shininess = p.shininess;
    }
    intersection_point_property(point P, Vector N, double t)
    {
        this->P = P;
        this->N = N;
        this->t = t;
        this->red = 0;
        this->green = 0;
        this->blue = 0;
        this->ambient = 0;
        this->diffuse = 0;
        this->specular = 0;
        this->reflection = 0;
        this->shininess = 0;
    }
};

class cube
{
public:
    point bottom_lower_left;
    double side;
    double ambient, diffuse, specular, reflection, shininess;
    double red, green, blue;
    cube(point bottom_lower_left, double side)
    {
        this->bottom_lower_left = bottom_lower_left;
        this->side = side;
    }
    cube(const cube &c)
    {
        this->bottom_lower_left = c.bottom_lower_left;
        this->side = c.side;
        this->ambient = c.ambient;
        this->diffuse = c.diffuse;
        this->specular = c.specular;
        this->reflection = c.reflection;
        this->shininess = c.shininess;
        this->red = c.red;
        this->green = c.green;
        this->blue = c.blue;
    }
    void set_color(double red, double green, double blue)
    {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }
    void set_color_component(double ambient, double diffuse, double specular, double reflection, double shininess)
    {
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->reflection = reflection;
        this->shininess = shininess;
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
    double ambient, diffuse, specular, reflection, shininess;
    double red, green, blue;
    pyramid(point bottom_center, double side, double height)
    {
        this->bottom_center = bottom_center;
        this->side = side;
        this->height = height;
    }
    void set_color(double red, double green, double blue)
    {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }
    void set_color_component(double ambient, double diffuse, double specular, double reflection, double shininess)
    {
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->reflection = reflection;
        this->shininess = shininess;
    }
};

double find_intersection_point_with_plane(point R_o, Vector ray, Plane plane)
{
    double t = (plane.P - R_o).dot(plane.N) / ray.dot(plane.N);
    return t;
}

intersection_point_property find_intersection_point_with_cube(point R_o, Vector ray, cube c)
{
    Plane planes[6];
    double t;
    planes[0] = Plane(c.bottom_lower_left, Vector(0, 0, -1));
    t = find_intersection_point_with_plane(R_o, ray, planes[0]);
    // double min_t = -10;
    intersection_point_property ip = intersection_point_property(point(0, 0, 0), Vector(0, 0, 0), -10);
    if (t > 0)
    {
        point intersection_point = R_o + ray * t;
        if (intersection_point.x >= c.bottom_lower_left.x && intersection_point.x <= c.bottom_lower_left.x + c.side && intersection_point.y >= c.bottom_lower_left.y && intersection_point.y <= c.bottom_lower_left.y + c.side)
        {
            // if (min_t == -10 || t < min_t)
            //     min_t = t;
            if (ip.t == -10 || t < ip.t)
            {
                ip.P = intersection_point;
                ip.N = planes[0].N;
                ip.t = t;
            }
        }
    }
    planes[1] = Plane(c.bottom_lower_left, Vector(0, -1, 0));
    t = find_intersection_point_with_plane(R_o, ray, planes[1]);
    if (t > 0)
    {
        point intersection_point = R_o + ray * t;
        if (intersection_point.x >= c.bottom_lower_left.x && intersection_point.x <= c.bottom_lower_left.x + c.side && intersection_point.z >= c.bottom_lower_left.z && intersection_point.z <= c.bottom_lower_left.z + c.side)
        {
            if (ip.t == -10 || t < ip.t)
            {
                ip.P = intersection_point;
                ip.N = planes[1].N;
                ip.t = t;
            }
        }
    }
    planes[2] = Plane(c.bottom_lower_left, Vector(-1, 0, 0));
    t = find_intersection_point_with_plane(R_o, ray, planes[2]);
    if (t > 0)
    {
        point intersection_point = R_o + ray * t;
        if (intersection_point.y >= c.bottom_lower_left.y && intersection_point.y <= c.bottom_lower_left.y + c.side && intersection_point.z >= c.bottom_lower_left.z && intersection_point.z <= c.bottom_lower_left.z + c.side)
        {
            // if (min_t == -10 || t < min_t)
            //     min_t = t;
            if (ip.t == -10 || t < ip.t)
            {
                ip.P = intersection_point;
                ip.N = planes[2].N;
                ip.t = t;
            }
        }
    }
    planes[3] = Plane(c.bottom_lower_left + Vector(0, 0, c.side), Vector(0, 0, 1));
    t = find_intersection_point_with_plane(R_o, ray, planes[3]);
    if (t > 0)
    {
        point intersection_point = R_o + ray * t;
        if (intersection_point.x >= c.bottom_lower_left.x && intersection_point.x <= c.bottom_lower_left.x + c.side && intersection_point.y >= c.bottom_lower_left.y && intersection_point.y <= c.bottom_lower_left.y + c.side)
        {
            // if (min_t == -10 || t < min_t)
            //     min_t = t;
            if (ip.t == -10 || t < ip.t)
            {
                ip.P = intersection_point;
                ip.N = planes[3].N;
                ip.t = t;
            }
        }
    }
    planes[4] = Plane(c.bottom_lower_left + Vector(0, c.side, 0), Vector(0, 1, 0));
    t = find_intersection_point_with_plane(R_o, ray, planes[4]);
    if (t > 0)
    {
        point intersection_point = R_o + ray * t;
        if (intersection_point.x >= c.bottom_lower_left.x && intersection_point.x <= c.bottom_lower_left.x + c.side && intersection_point.z >= c.bottom_lower_left.z && intersection_point.z <= c.bottom_lower_left.z + c.side)
        {
            // if (min_t == -10 || t < min_t)
            //     min_t = t;
            if (ip.t == -10 || t < ip.t)
            {
                ip.P = intersection_point;
                ip.N = planes[4].N;
                ip.t = t;
            }
        }
    }
    planes[5] = Plane(c.bottom_lower_left + Vector(c.side, 0, 0), Vector(1, 0, 0));
    t = find_intersection_point_with_plane(R_o, ray, planes[5]);
    if (t > 0)
    {
        point intersection_point = R_o + ray * t;
        if (intersection_point.y >= c.bottom_lower_left.y && intersection_point.y <= c.bottom_lower_left.y + c.side && intersection_point.z >= c.bottom_lower_left.z && intersection_point.z <= c.bottom_lower_left.z + c.side)
        {
            // if (min_t == -10 || t < min_t)
            //     min_t = t;
            if (ip.t == -10 || t < ip.t)
            {
                ip.P = intersection_point;
                ip.N = planes[5].N;
                ip.t = t;
            }
        }
    }
    if (ip.t != -10)
    {
        ip.red = c.red;
        ip.green = c.green;
        ip.blue = c.blue;
        ip.ambient = c.ambient;
        ip.diffuse = c.diffuse;
        ip.specular = c.specular;
        ip.reflection = c.reflection;
        ip.shininess = c.shininess;
    }
    return ip;
}

class Sphere
{
public:
    point center;
    double radius;
    double ambient, diffuse, specular, reflection, shininess;
    double red, green, blue;
    Sphere(point center, double radius)
    {
        this->center = center;
        this->radius = radius;
    }
    void set_color(double red, double green, double blue)
    {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }
    void set_color_component(double ambient, double diffuse, double specular, double reflection, double shininess)
    {
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->reflection = reflection;
        this->shininess = shininess;
    }
};

class checkerboard
{
public:
    double side;
    double ambient, diffuse, specular, reflection, shininess;
    double red, green, blue;
    checkerboard()
    {
        this->side = 0;
        this->specular = 0;
        this->shininess = 0;
    }
    checkerboard(checkerboard &c)
    {
        this->side = c.side;
        this->ambient = c.ambient;
        this->diffuse = c.diffuse;
        this->reflection = c.reflection;
        this->specular = c.specular;
        this->shininess = c.shininess;
        this->red = c.red;
        this->green = c.green;
        this->blue = c.blue;
    }
    checkerboard(double side)
    {
        this->side = side;
        this->specular = 0;
        this->shininess = 0;
        this->red = 0;
        this->green = 0;
        this->blue = 0;
    }
    void set_color_component(double ambient, double diffuse, double reflection)
    {
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->reflection = reflection;
        this->specular = 0;
        this->shininess = 0;
    }
    void set_color(point temp)
    {
        int x = abs(temp.x) / side;
        int y = abs(temp.y) / side;
        double x_dist, y_dist;
        if ((x + y) % 2 == 0)
        {
            if (temp.x * temp.y > 0)
            {
                this->red = 1, this->green = 1, this->blue = 1;
                if (temp.x > 0)
                {
                    x_dist = temp.x - x * side;
                    y_dist = side - (temp.y - y * side);
                }
                else
                {
                    x_dist = side - (-x * side - temp.x);
                    y_dist = -y * side - temp.y;
                }
            }
            else
            {
                this->red = 0, this->green = 0, this->blue = 0;
                if (temp.x > 0)
                {
                    x_dist = temp.x - x * side;
                    y_dist = -y * side - temp.y;
                }
                else
                {
                    x_dist = side - (-x * side - temp.x);
                    y_dist = side - (temp.y - y * side);
                }
            }
        }
        else
        {
            if (temp.x * temp.y > 0)
            {
                this->red = 0, this->green = 0, this->blue = 0;
                if (temp.x > 0)
                {
                    x_dist = temp.x - x * side;
                    y_dist = side - (temp.y - y * side);
                }
                else
                {
                    x_dist = side - (-x * side - temp.x);
                    y_dist = -y * side - temp.y;
                }
            }
            else
            {
                this->red = 1, this->green = 1, this->blue = 1;
                if (temp.x > 0)
                {
                    x_dist = temp.x - x * side;
                    y_dist = -y * side - temp.y;
                }
                else
                {
                    x_dist = side - (-x * side - temp.x);
                    y_dist = side - (temp.y - y * side);
                }
            }
        }
        if (this->red == 0 && texture)
        {
            int x_pixel = int(texture_b->width() / side * x_dist);
            int y_pixel = int(texture_b->height() / side * y_dist);
            this->red = texture_b_array[y_pixel][x_pixel][0];
            this->green = texture_b_array[y_pixel][x_pixel][1];
            this->blue = texture_b_array[y_pixel][x_pixel][2];
        }
        else if (this->red == 1 && texture)
        {
            int x_pixel = int(texture_w->width() / side * x_dist);
            int y_pixel = int(texture_w->height() / side * y_dist);
            this->red = texture_w_array[y_pixel][x_pixel][0];
            this->green = texture_w_array[y_pixel][x_pixel][1];
            this->blue = texture_w_array[y_pixel][x_pixel][2];
        }
    }
};

class light_source
{
public:
    point p;
    point direction;
    double angle;
    double fall_off;
    bool normal;
    light_source(point p, point direction, double angle, double fall_off, bool normal)
    {
        this->p = p;
        this->direction = direction;
        this->angle = angle;
        this->fall_off = fall_off;
        this->normal = normal;
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

class objects
{
public:
    vector<cube> cubes;
    vector<pyramid> pyramids;
    vector<Sphere> spheres;
    checkerboard cb;
    objects()
    {
        cubes.clear();
        pyramids.clear();
        spheres.clear();
    }
    void add_cube(cube c)
    {
        cubes.push_back(c);
    }
    void add_pyramid(pyramid p)
    {
        pyramids.push_back(p);
    }
    void add_sphere(Sphere s)
    {
        spheres.push_back(s);
    }
    void add_checkerboard(checkerboard c)
    {
        cb = c;
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

intersection_point_property find_intersection_point_with_pyramid(point R_o, Vector ray, pyramid p)
{
    triangle t1(p.bottom_center + Vector(-p.side / 2, -p.side / 2, 0), p.bottom_center + Vector(p.side / 2, -p.side / 2, 0), p.bottom_center + Vector(0, 0, p.height));
    triangle t2(p.bottom_center + Vector(p.side / 2, -p.side / 2, 0), p.bottom_center + Vector(p.side / 2, p.side / 2, 0), p.bottom_center + Vector(0, 0, p.height));
    triangle t3(p.bottom_center + Vector(p.side / 2, p.side / 2, 0), p.bottom_center + Vector(-p.side / 2, p.side / 2, 0), p.bottom_center + Vector(0, 0, p.height));
    triangle t4(p.bottom_center + Vector(-p.side / 2, p.side / 2, 0), p.bottom_center + Vector(-p.side / 2, -p.side / 2, 0), p.bottom_center + Vector(0, 0, p.height));
    double t1_val = find_intersection_point_with_triangle(R_o, ray, t1);
    double t2_val = find_intersection_point_with_triangle(R_o, ray, t2);
    double t3_val = find_intersection_point_with_triangle(R_o, ray, t3);
    double t4_val = find_intersection_point_with_triangle(R_o, ray, t4);
    double t5_val = find_intersection_point_with_plane(R_o, ray, Plane(p.bottom_center, Vector(0, 0, -1)));
    // double min_t = -10;
    intersection_point_property ip = intersection_point_property(point(0, 0, 0), Vector(0, 0, 0), -10);
    if (t1_val > 0)
    {
        // min_t = t1_val;
        if (ip.t == -10 || t1_val < ip.t)
        {
            ip.P = R_o + ray * t1_val;
            ip.N = (t1.A - t1.B) * (t1.C - t1.B) * (-1);
            ip.N = (ip.N.normalize());
            ip.t = t1_val;
        }
    }
    if (t2_val > 0)
    {
        // if (min_t == -10 || t2_val < min_t)
        //     min_t = t2_val;
        if (ip.t == -10 || t2_val < ip.t)
        {
            ip.P = R_o + ray * t2_val;
            ip.N = (t2.A - t2.B) * (t2.C - t2.B) * (-1);
            ip.N = (ip.N.normalize());
            ip.t = t2_val;
        }
    }
    if (t3_val > 0)
    {
        // if (min_t == -10 || t3_val < min_t)
        //     min_t = t3_val;
        if (ip.t == -10 || t3_val < ip.t)
        {
            ip.P = R_o + ray * t3_val;
            ip.N = (t3.A - t3.B) * (t3.C - t3.B) * (-1);
            ip.N = (ip.N.normalize());
            ip.t = t3_val;
        }
    }
    if (t4_val > 0)
    {
        // if (min_t == -10 || t4_val < min_t)
        //     min_t = t4_val;
        if (ip.t == -10 || t4_val < ip.t)
        {
            ip.P = R_o + ray * t4_val;
            ip.N = (t4.A - t4.B) * (t4.C - t4.B) * (-1);
            ip.N = (ip.N.normalize());
            ip.t = t4_val;
        }
    }
    if (t5_val > 0)
    {
        point temp = R_o + ray * t5_val;
        if (temp.x >= p.bottom_center.x - p.side / 2 && temp.x <= p.bottom_center.x + p.side / 2 && temp.y >= p.bottom_center.y - p.side / 2 && temp.y <= p.bottom_center.y + p.side / 2)
            // if (min_t == -10 || t5_val < min_t)
            //     min_t = t5_val;
            if (ip.t == -10 || t5_val < ip.t)
            {
                ip.P = temp;
                ip.N = Vector(0, 0, -1);
                ip.t = t5_val;
            }
    }
    if (ip.t != -10)
    {
        ip.red = p.red;
        ip.green = p.green;
        ip.blue = p.blue;
        ip.ambient = p.ambient;
        ip.diffuse = p.diffuse;
        ip.specular = p.specular;
        ip.reflection = p.reflection;
        ip.shininess = p.shininess;
    }
    return ip;
}

intersection_point_property find_intersection_point_with_sphere(point R_o, Vector ray, Sphere s)
{
    double a = (R_o - s.center).dot(R_o - s.center);
    double tp = -ray.dot(R_o - s.center);
    intersection_point_property ip(point(0, 0, 0), Vector(0, 0, 0), -10);
    if (a > 0 && tp < 0)
    {
        return intersection_point_property(point(0, 0, 0), Vector(0, 0, 0), -10);
    }
    double d2 = a - tp * tp;
    if (d2 > s.radius * s.radius)
    {
        return intersection_point_property(point(0, 0, 0), Vector(0, 0, 0), -10);
    }
    double t = sqrt(s.radius * s.radius - d2);

    if (a > s.radius * s.radius)
    {
        // return tp - t;
        point temp = R_o + ray * (tp - t);
        Vector temp2 = temp - s.center;
        temp2 = temp2.normalize();
        ip = intersection_point_property(temp, temp2, tp - t);
        ip.red = s.red;
        ip.green = s.green;
        ip.blue = s.blue;
        ip.ambient = s.ambient;
        ip.diffuse = s.diffuse;
        ip.specular = s.specular;
        ip.reflection = s.reflection;
        ip.shininess = s.shininess;
    }
    else if (a < s.radius * s.radius)
    {
        // return tp + t;
        point temp = R_o + ray * (tp + t);
        Vector temp2 = temp - s.center;
        temp2 = temp2.normalize();
        ip = intersection_point_property(temp, temp2, tp + t);
        ip.red = s.red;
        ip.green = s.green;
        ip.blue = s.blue;
        ip.ambient = s.ambient;
        ip.diffuse = s.diffuse;
        ip.specular = s.specular;
        ip.reflection = s.reflection;
        ip.shininess = s.shininess;
    }
    else
    {
        // return tp + 0.001;
        point temp = R_o + ray * (tp);
        Vector temp2 = temp - s.center;
        temp2 = temp2.normalize();
        ip = intersection_point_property(temp, temp2, tp);
        ip.red = s.red;
        ip.green = s.green;
        ip.blue = s.blue;
        ip.ambient = s.ambient;
        ip.diffuse = s.diffuse;
        ip.specular = s.specular;
        ip.reflection = s.reflection;
        ip.shininess = s.shininess;
    }
    return ip;
}

intersection_point_property find_intersection_point_with_checkerboard(point R_o, Vector ray, checkerboard c)
{
    Plane plane(point(0, 0, 0), Vector(0, 0, 1));
    double t = find_intersection_point_with_plane(R_o, ray, plane);
    intersection_point_property ip(point(0, 0, 0), Vector(0, 0, 0), -10);
    if (t > 0)
    {
        // return t;
        point temp = R_o + ray * t;
        c.set_color(temp);
        if (R_o.z < 0)
            ip = intersection_point_property(temp, plane.N * (-1), t);
        else
            ip = intersection_point_property(temp, plane.N, t);
        ip.red = c.red;
        ip.green = c.green;
        ip.blue = c.blue;
        ip.ambient = c.ambient;
        ip.diffuse = c.diffuse;
        ip.reflection = c.reflection;
        ip.specular = c.specular;
        ip.shininess = c.shininess;
    }
    return ip;
}

intersection_point_property find_diffuse_specular_component(vector<light_source> ls_list, objects objs, intersection_point_property ip, point camera)
{
    double lambert = 0, phong = 0;
    for (light_source ls : ls_list)
    {
        Vector L = ls.p - ip.P;
        L = L.normalize();
        ip.P = ip.P + ip.N * 0.0001;
        bool flag = false;
        for (cube c : objs.cubes)
        {
            intersection_point_property ip2 = find_intersection_point_with_cube(ip.P, L, c);
            if (ip2.t == -10 || ip2.t < 0)
            {
                continue;
            }
            else if (ip2.t < (ls.p - ip.P).length())
            {
                flag = true;
                break;
            }
            else
            {
                continue;
            }
        }
        if (flag)
            continue;
        for (pyramid p : objs.pyramids)
        {
            intersection_point_property ip2 = find_intersection_point_with_pyramid(ip.P, L, p);
            if (ip2.t == -10 || ip2.t < 0)
            {
                continue;
            }
            else if (ip2.t < (ls.p - ip.P).length())
            {
                flag = true;
                break;
            }
            else
            {
                continue;
            }
        }
        if (flag)
            continue;
        for (Sphere s : objs.spheres)
        {
            intersection_point_property ip2 = find_intersection_point_with_sphere(ip.P, L, s);
            if (ip2.t == -10 || ip2.t < 0)
            {
                // cout<<ip.P.x<<" "<<ip.P.y<<" "<<ip.P.z<<endl;

                continue;
            }
            else if (ip2.t < (ls.p - ip.P).length())
            {
                flag = true;
                break;
            }
            else
            {
                continue;
            }
        }
        if (flag)
            continue;
        intersection_point_property ip2 = find_intersection_point_with_checkerboard(ip.P, L, objs.cb);
        if (ip2.t == -10 || ip2.t < 0)
        {
            flag = false;
        }
        else if (ip2.t < (ls.p - ip.P).length())
        {
            flag = true;
        }
        else
        {
            flag = false;
        }
        if (flag)
            continue;
        if (ls.normal == false)
        {
            Vector V1 = ip.P - ls.p;
            V1 = V1.normalize();
            Vector V2 = ls.direction - ls.p;
            V2 = V2.normalize();
            double angle = acos(V1.dot(V2));
            if (angle > ls.angle)
            {
                continue;
            }
            else
            {
                // cout<<angle<<endl;
            }
        }
        double distance = L.length();
        double scaling_factor = exp(-ls.fall_off * distance * distance);
        if (L.dot(ip.N) > 0)
        {
            lambert += scaling_factor * (L.dot(ip.N) > 0 ? L.dot(ip.N) : 0);
            Vector R = (ip.P - camera) - ip.N * (2 * (ip.P - camera).dot(ip.N));
            R = R.normalize();
            phong += scaling_factor * pow((R.dot(L) > 0 ? R.dot(L) : 0), ip.shininess);
        }
        else
        {
        }
    }
    ip.red = ip.red * (ip.ambient + ip.diffuse * lambert + ip.specular * phong);
    ip.green = ip.green * (ip.ambient + ip.diffuse * lambert + ip.specular * phong);
    ip.blue = ip.blue * (ip.ambient + ip.diffuse * lambert + ip.specular * phong);

    return ip;
}

intersection_point_property find_color(Vector ray, objects objs, point Ro, vector<light_source> ls_list, int level)
{
    intersection_point_property ip = intersection_point_property(point(0, 0, 0), Vector(0, 0, 0), -10);

    if (level == 0)
    {
        ip.red = 0;
        ip.green = 0;
        ip.blue = 0;
        return ip;
    }
    for (cube c : objs.cubes)
    {
        intersection_point_property t = find_intersection_point_with_cube(Ro, ray, c);
        if (t.t > 0)
        {
            if (ip.t == -10 || t.t < ip.t)
            {
                ip = t;
            }
        }
    }
    for (Sphere s : objs.spheres)
    {
        intersection_point_property t = find_intersection_point_with_sphere(Ro, ray, s);
        if (t.t > 0)
        {
            if (ip.t == -10 || t.t < ip.t)
            {
                ip = t;
            }
        }
    }
    for (pyramid p : objs.pyramids)
    {
        intersection_point_property t = find_intersection_point_with_pyramid(Ro, ray, p);
        if (t.t > 0)
        {
            if (ip.t == -10 || t.t < ip.t)
            {
                ip = t;
            }
        }
    }
    intersection_point_property t = find_intersection_point_with_checkerboard(Ro, ray, objs.cb);
    if (t.t > 0)
    {
        if (ip.t == -10 || t.t < ip.t)
        {
            ip = t;
        }
    }

    if (ip.t == -10)
    {
        ip.red = 0;
        ip.green = 0;
        ip.blue = 0;
        return ip;
    }
    if (ip.t != -10)
    {
        ip = find_diffuse_specular_component(ls_list, objs, ip, Ro);
        Vector R = (ip.P - Ro) - ip.N * (2 * (ip.P - Ro).dot(ip.N));
        R = R.normalize();
        ip.P = ip.P + ip.N * 0.0001;
        intersection_point_property ip2 = find_color(R, objs, ip.P, ls_list, level - 1);
        ip.red += ip2.red * ip.reflection;
        ip.green += ip2.green * ip.reflection;
        ip.blue += ip2.blue * ip.reflection;
    }
    return ip;
}

void create_pointBuffer(double near_plane_d, double far_plane_d, double fovy, double aspect_ratio, Vector l, Vector r, Vector u, point camera, int image_width, int image_height, int recursion_level, objects objs, vector<light_source> ls_list)
{
    fovy = fovy * acos(-1) / 180;
    double near_plane_h = 2 * near_plane_d * tan(fovy / 2);
    double fovx = aspect_ratio * fovy;
    double near_plane_w = 2 * near_plane_d * tan(fovx / 2);
    point midpoint = camera + l * near_plane_d;
    double w_per_pixel = near_plane_w / image_width;
    double h_per_pixel = near_plane_h / image_height;
    point topleft = midpoint + u * (near_plane_h / 2 - h_per_pixel / 2) - r * (near_plane_w / 2 - w_per_pixel / 2);
    point **pointBuffer = new point *[image_width];
    
    for (int i = 0; i < image_width; i++)
    {
        pointBuffer[i] = new point[image_height];
    }
    for (int i = 0; i < image_height; i++)
    {
        for (int j = 0; j < image_width; j++)
        {
            pointBuffer[j][i] = topleft + u * (-h_per_pixel * i) + r * (w_per_pixel * j);
        }
    }
    
    if (texture)
    {
        texture_w = new bitmap_image("texture_w.bmp");
        texture_b = new bitmap_image("texture_b.bmp");

        // load above images to 2 2d arrays
        texture_w_array = new double **[texture_w->height() + 1];
        for (int i = 0; i <= texture_w->height(); i++)
        {
            texture_w_array[i] = new double *[texture_w->width() + 1];
            for (int j = 0; j <= texture_w->width(); j++)
            {
                texture_w_array[i][j] = new double[3];
            }
        }
        for (int i = 0; i <= texture_w->height(); i++)
        {
            for (int j = 0; j <= texture_w->width(); j++)
            {
                texture_w_array[i][j][0] = texture_w->red_channel(j, i) / 255.0;
                texture_w_array[i][j][1] = texture_w->green_channel(j, i) / 255.0;
                texture_w_array[i][j][2] = texture_w->blue_channel(j, i) / 255.0;
            }
        }
        texture_b_array = new double **[texture_b->height() + 1];
        for (int i = 0; i <= texture_b->height(); i++)
        {
            texture_b_array[i] = new double *[texture_b->width() + 1];
            for (int j = 0; j <= texture_b->width(); j++)
            {
                texture_b_array[i][j] = new double[3];
            }
        }
        for (int i = 0; i <= texture_b->height(); i++)
        {
            for (int j = 0; j <= texture_b->width(); j++)
            {
                texture_b_array[i][j][0] = texture_b->red_channel(j, i) / 255.0;
                texture_b_array[i][j][1] = texture_b->green_channel(j, i) / 255.0;
                texture_b_array[i][j][2] = texture_b->blue_channel(j, i) / 255.0;
            }
        }
    }

    bitmap_image *image;
    image = new bitmap_image(image_width, image_height);
    for (int i = 0; i < image_height; i++)
        for (int j = 0; j < image_width; j++)
            image->set_pixel(j, i, 0, 0, 0);
    for (int i = 0; i < image_height; i++)
    {
        for (int j = 0; j < image_width; j++)
        {
            Vector ray = pointBuffer[j][i] - camera;
            ray = ray.normalize();
            intersection_point_property ip = find_color(ray, objs, pointBuffer[j][i], ls_list, recursion_level);
            if (ip.t != -10)
            {
                if (ip.t * ray.dot(l) <= (far_plane_d - near_plane_d))
                    image->set_pixel(j, i, ip.red * 255, ip.green * 255, ip.blue * 255);
            }
        }
    }

    cout << "done" << endl;
    image->save_image("output.bmp");

    // clear pointBuffer
    for (int i = 0; i < image_width; i++)
    {
        delete[] pointBuffer[i];
    }
    delete[] pointBuffer;
    // clear image pointer
    delete image;
    // clear texture pointer
    // clear texture array
    if (texture)
    {
        for (int i = 0; i <= texture_b->height(); i++)
        {
            for (int j = 0; j <= texture_b->width(); j++)
            {
                delete[] texture_b_array[i][j];
            }
            delete[] texture_b_array[i];
        }
        delete[] texture_b_array;
        for (int i = 0; i <= texture_w->height(); i++)
        {
            for (int j = 0; j <= texture_w->width(); j++)
            {
                delete[] texture_w_array[i][j];
            }
            delete[] texture_w_array[i];
        }
        delete[] texture_w_array;

        delete texture_w;
        delete texture_b;
    }
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
