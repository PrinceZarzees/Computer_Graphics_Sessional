#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stack>
#include <cmath>
#include <iomanip>
#include <vector>
#include "bitmap_image.hpp"
using namespace std;
class Point
{
public:
    double x, y, z;
    Point()
    {
        x = y = z = 0;
    }
    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};
class Vector
{
public:
    double x, y, z;
    Vector()
    {
        x = y = z = 0;
    }
    Vector(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void normalize()
    {
        double length = sqrt(x * x + y * y + z * z);
        x /= length;
        y /= length;
        z /= length;
    }
    Vector operator*(const Vector &other) const
    {
        Vector ret;
        ret.x = y * other.z - z * other.y;
        ret.y = z * other.x - x * other.z;
        ret.z = x * other.y - y * other.x;
        return ret;
    }
    Vector operator*(const double &other) const
    {
        Vector ret;
        ret.x = x * other;
        ret.y = y * other;
        ret.z = z * other;
        return ret;
    }
    Vector operator+(const Vector &other) const
    {
        Vector ret;
        ret.x = x + other.x;
        ret.y = y + other.y;
        ret.z = z + other.z;
        return ret;
    }
    Vector operator-(const Vector &other) const
    {
        Vector ret;
        ret.x = x - other.x;
        ret.y = y - other.y;
        ret.z = z - other.z;
        return ret;
    }
    double dot(const Vector &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }
};
Vector R(Vector x, Vector a, double theta)
{
    Vector ret = x * cos(theta) + a * (a.dot(x)) * (1 - cos(theta)) + (a * x) * sin(theta);
    return ret;
}
class Matrix
{
public:
    double mat[4][4];
    Matrix()
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                if (i == j)
                    mat[i][j] = 1;
                else
                    mat[i][j] = 0;
            }
    }
    Matrix(double mat[4][4])
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                this->mat[i][j] = mat[i][j];
            }
    }
    Matrix operator*(const Matrix &other) const
    {
        Matrix ret;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                ret.mat[i][j] = 0;
                for (int k = 0; k < 4; k++)
                    ret.mat[i][j] += mat[i][k] * other.mat[k][j];
            }
        return ret;
    }
    Point operator*(const Point &p) const
    {
        Point ret;
        ret.x = mat[0][0] * p.x + mat[0][1] * p.y + mat[0][2] * p.z + mat[0][3];
        ret.y = mat[1][0] * p.x + mat[1][1] * p.y + mat[1][2] * p.z + mat[1][3];
        ret.z = mat[2][0] * p.x + mat[2][1] * p.y + mat[2][2] * p.z + mat[2][3];
        double w = mat[3][0] * p.x + mat[3][1] * p.y + mat[3][2] * p.z + mat[3][3];
        ret.x /= w;
        ret.y /= w;
        ret.z /= w;
        return ret;
    }
};
Matrix translation_matrix(double tx, double ty, double tz)
    {
        double mat[4][4] = {{1, 0, 0, tx}, {0, 1, 0, ty}, {0, 0, 1, tz}, {0, 0, 0, 1}};
        return Matrix(mat);
    }
Matrix rotation_matrix(double angle, double ax, double ay, double az)
    {
        Vector a(ax, ay, az);
        a.normalize();
        Vector c1 = R(Vector(1, 0, 0), a, angle);
        Vector c2 = R(Vector(0, 1, 0), a, angle);
        Vector c3 = R(Vector(0, 0, 1), a, angle);
        double mat[4][4] = {{c1.x, c2.x, c3.x, 0}, {c1.y, c2.y, c3.y, 0}, {c1.z, c2.z, c3.z, 0}, {0, 0, 0, 1}};
        return Matrix(mat);
    }
Matrix scale_matrix(double sx, double sy, double sz)
    {
        double mat[4][4] = {{sx, 0, 0, 0}, {0, sy, 0, 0}, {0, 0, sz, 0}, {0, 0, 0, 1}};
        return Matrix(mat);
    }


void stage1(ifstream &inputFile)
{
    stack<Matrix> st;
    Matrix M;
    string line;
    ofstream outputFile("stage1.txt");

    while (1)
    {
        getline(inputFile, line);
        istringstream iss(line);
        string command;
        iss >> command;
        if (command == "triangle")
        {
            Point p[3];
            for(int i=0;i<3;i++)
                {getline(inputFile, line);
                istringstream iss(line);
                iss >> p[i].x >> p[i].y >> p[i].z;
                }
            for (int i=0;i<3;i++)
                {
                    Point p1 = M * p[i];
                    outputFile <<fixed<<setprecision(7);
                    outputFile << p1.x << " " << p1.y << " " << p1.z << endl;
                }
            outputFile<<endl;
        }
        else if (command == "translate")
        {
            getline(inputFile, line);
            istringstream iss(line);
            double tx, ty, tz;
            iss >> tx >> ty >> tz;
            Matrix T = translation_matrix(tx, ty, tz);
            M = M * T;
        }
        else if (command == "scale")
        {
            getline(inputFile, line);
            istringstream iss(line);
            double sx, sy, sz;
            iss >> sx >> sy >> sz;
            Matrix S = scale_matrix(sx, sy, sz);
            M = M * S;
        }
        else if (command == "rotate")
        {
            getline(inputFile, line);
            istringstream iss(line);
            double angle, ax, ay, az;
            iss >> angle >> ax >> ay >> az;
            angle=angle*acos(-1)/180;
            Matrix R = rotation_matrix(angle, ax, ay, az);
            M = M * R;
        }
        else if (command == "push")
        {
            st.push(M);
        }
        else if (command == "pop")
        {
            M = st.top();
            st.pop();
        }
        else if (command == "end")
        {
            break;
        }
    }
    outputFile.close();
}
void stage2(ifstream &inputFile, Vector eye, Vector look, Vector up)
{
    string line;
    Vector l = look - eye;
    l.normalize();
    Vector r = l * up;
    r.normalize();
    Vector u = r * l;
    Matrix T = translation_matrix(-eye.x,-eye.y,-eye.z);
    double mat[4][4]= {{r.x, r.y, r.z, 0}, {u.x, u.y, u.z, 0}, {-l.x, -l.y, -l.z, 0}, {0, 0, 0, 1}};
    Matrix R = Matrix(mat);
    Matrix V = R * T;

    ofstream outputFile("stage2.txt");
    while (getline(inputFile, line))
    {
        if (line.length() == 0)
        {    
           outputFile << endl;
           continue;
        }
    istringstream iss(line);
    Point p;
    iss >> p.x >> p.y >> p.z;
    Point p1 = V * p;
    outputFile <<fixed<<setprecision(7);
    outputFile << p1.x << " " << p1.y << " " << p1.z << endl;


    }

    outputFile.close();


}
void stage3(ifstream &inputFile, double fovY, double aspectRatio, double near, double far)
{
    double fovX = fovY * aspectRatio;
    double t = near * tan(fovY*acos(-1)/360);
    double r = near * tan(fovX*acos(-1)/360);
    Matrix P = Matrix();
    P.mat[0][0] = near / r;
    P.mat[1][1] = near / t;
    P.mat[2][2] = -(far + near) / (far - near);
    P.mat[2][3] = -(2 * far * near) / (far - near);
    P.mat[3][2] = -1;
    P.mat[3][3] = 0;
    string line;
    ofstream outputFile("stage3.txt");
    while (getline(inputFile, line))
    {
        if (line.length() == 0)
        {
            outputFile << endl;
            continue;
        }
        istringstream iss(line);
        Point p;
        iss >> p.x >> p.y >> p.z;
        Point p1 = P * p;
        outputFile <<fixed<<setprecision(7);
        outputFile << p1.x << " " << p1.y << " " << p1.z << endl;
    }
    outputFile.close();
}

static unsigned long int g_seed = 1;
inline int random()
{
 g_seed = (214013 * g_seed + 2531011);
 return (g_seed >> 16) & 0x7FFF;
}

class triangle
{
public:
    Point p[3];
    double red, green, blue;
    triangle(Point p1, Point p2, Point p3)
    {
        p[0] = p1;
        p[1] = p2;
        p[2] = p3;
    }

};
void z_buffer(ifstream &inputFile)
{
    //Read the config.txt file and store the values as Screen_Width, Screen_Height
    int Screen_Width, Screen_Height;
    ifstream config_file("1/config.txt");
    string line;
    getline(config_file, line);
    istringstream iss(line);
    iss >> Screen_Width >> Screen_Height;
    config_file.close();
    //Read stage3.txt and Use a suitable data structure to hold this information
    vector<triangle> triangles;
    while (getline(inputFile, line))
    {
        if (line.length() == 0)
            continue;
        Point p[3];
        for (int i=0;i<3;i++)
          { 
            istringstream iss(line);
            iss >> p[i].x >> p[i].y >> p[i].z;
            if (i==2)
                break;
            getline(inputFile, line);
          }
        triangle t(p[0], p[1], p[2]);
        t.red = random()%255;
        t.green = random()%255;
        t.blue = random()%255;
        triangles.push_back(t);
    //             for (int j=0;j<3;j++)
    //                 {
    //                     Point p1 = triangles[triangles.size()-1].p[j];
    //                     cout<<p1.x<<" "<<p1.y<<" "<<p1.z<<endl;
    //                 }
    //             cout<<triangles[triangles.size()-1].blue<<" "<<triangles[triangles.size()-1].green<<" "<<triangles[triangles.size()-1].red<<endl;
    //
    }
    inputFile.close();
    double right_limit=1.0;
    double left_limit=-1.0;
    double top_limit=1.0;
    double bottom_limit=-1.0;
    double dx=(right_limit-left_limit)/Screen_Width;
    double dy=(top_limit-bottom_limit)/Screen_Height;
    double Top_Y=top_limit-dy/2;
    double Bottom_Y=bottom_limit+dy/2;
    double Left_X=left_limit+dx/2;
    double Right_X=right_limit-dx/2;
    //cout<<dx<<" "<<dy<<endl;
    //cout<<Top_Y<<" "<<Bottom_Y<<" "<<Left_X<<" "<<Right_X<<endl;
    //Initialize the Z-buffer
    double **Z_buffer = new double *[Screen_Height];
    double z_max=1;
    double z_front_limit=-1;
    for (int i = 0; i < Screen_Height; i++)
    {
        Z_buffer[i] = new double[Screen_Width];
        for (int j = 0; j < Screen_Width; j++)
            Z_buffer[i][j] = z_max;
    }
    bitmap_image image(Screen_Width,Screen_Height);
    for (int i=0;i<Screen_Height;i++)
        for (int j=0;j<Screen_Width;j++)
            image.set_pixel(j,i,0,0,0);

    for (triangle t: triangles)
        {
            double top_scanline,bottom_scanline;
            //compute max y & min_y
            double min_y=t.p[0].y;
            double max_y=t.p[0].y;
            for (int i=1;i<3;i++)
                {if (t.p[i].y>max_y)
                    max_y=t.p[i].y;
                if (t.p[i].y<min_y)
                    min_y=t.p[i].y;
                }
            //cout<<min_y<<" "<<max_y<<endl;
            
            if (max_y>Top_Y)
                top_scanline=0;
            else
                {
                    double diff=top_limit-max_y;
                    int row_no=diff/dy;
                    top_scanline=row_no;

                }
            if (min_y<Bottom_Y)
                bottom_scanline=Screen_Height-1;
            else
                {
                    double diff=top_limit-min_y;
                    int row_no=diff/dy;
                    bottom_scanline=row_no;
                }
             //cout<<top_scanline<<" "<<bottom_scanline<<endl;
            for (int row_no=top_scanline+1;row_no<bottom_scanline;row_no++)
            {
                double ys=Top_Y-row_no*dy;
                double min_x=t.p[0].x;
                double max_x=t.p[0].x;
                double za=t.p[0].z;
                double zb=t.p[0].z;
                for(int i=0;i<2;i++)
                    for(int j=i+1;j<3;j++)
                    {
                        double denom=(t.p[j].y-t.p[i].y);
                        if (denom==0)
                            continue;
                        double temp=(ys-t.p[i].y)/(t.p[j].y-t.p[i].y);
                        if (temp<0 || temp>1)
                            continue;
                        double x=(ys-t.p[i].y)*(t.p[j].x-t.p[i].x)/(t.p[j].y-t.p[i].y)+t.p[i].x;
                        double z=(ys-t.p[i].y)*(t.p[j].z-t.p[i].z)/(t.p[j].y-t.p[i].y)+t.p[i].z;
                        if (x<=min_x)
                            {min_x=x;
                            za=z;
                            }
                        else if (x>=max_x)
                            {max_x=x;
                            zb=z;
                            }
                    }
                    //cout<<za<<" "<<zb<<endl;
                double left_scanline,right_scanline;
                if (min_x<Left_X)
                    left_scanline=0;
                else
                    {
                        double diff=min_x-left_limit;
                        int col_no=diff/dx;
                        left_scanline=col_no;
                    }
                if (max_x>Right_X)
                    right_scanline=Screen_Width-1;
                else
                    {
                        double diff=max_x-left_limit;
                        int col_no=diff/dx;
                        right_scanline=col_no;
                    }
                // cout<<left_scanline<<" "<<right_scanline<<endl;
                for (int col_no=left_scanline;col_no<=right_scanline;col_no++)
                    {
                        double xp=Left_X+col_no*dx;
                        double xa=Left_X+left_scanline*dx;
                        double xb=Left_X+right_scanline*dx;
                        double denom=(xa-xb);
                        double zp;
                        if (denom==0)
                            {
                              zp=za;
                            }
                        else
                            zp=zb+(za-zb)*(xp-xb)/(xa-xb);
                        if (zp<z_front_limit)
                            {
                                continue;
                            }
                        if (zp<Z_buffer[row_no][col_no])
                            {
                                Z_buffer[row_no][col_no]=zp;
                                image.set_pixel(col_no,row_no,t.red,t.green,t.blue);
                            }

                    }
                

            }




        }

        //print the Z-buffer
        ofstream outputFile("z_buffer.txt");
        for (int i = 0; i < Screen_Height; i++)
        {
            for (int j = 0; j < Screen_Width; j++)
                if (Z_buffer[i][j]>=z_max)
                    outputFile << "";
                else
                    outputFile << Z_buffer[i][j] << " ";
            outputFile << endl;
        }
        outputFile.close();
        image.save_image("5/out_.bmp");



}
int main()
{
    ifstream inputFile("5/scene.txt");
    string line;
    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspectRatio, near, far;

    if (inputFile.is_open())
    {
        // Read line 1
        if (getline(inputFile, line))
        {
            istringstream iss(line);
            iss >> eyeX >> eyeY >> eyeZ;
        }

        // Read line 2
        if (getline(inputFile, line))
        {
            istringstream iss(line);
            iss >> lookX >> lookY >> lookZ;
        }

        // Read line 3
        if (getline(inputFile, line))
        {
            istringstream iss(line);
            iss >> upX >> upY >> upZ;
        }

        // Read line 4
        if (getline(inputFile, line))
        {
            istringstream iss(line);
            iss >> fovY >> aspectRatio >> near >> far;
        }

        stage1(inputFile);
        inputFile.close();
        inputFile= ifstream("stage1.txt");
        stage2(inputFile,Vector(eyeX,eyeY,eyeZ),Vector(lookX,lookY,lookZ),Vector(upX,upY,upZ));
        inputFile.close();
        inputFile= ifstream("stage2.txt");
        stage3(inputFile,fovY,aspectRatio,near,far);
        inputFile.close();
        inputFile= ifstream("stage3.txt");
        z_buffer(inputFile);
        inputFile.close();

    }
    else
    {
        cerr << "Failed to open the file." << endl;
    }

    //compare 2 text files
    // ifstream file1("stage2.txt");
    // ifstream file2("1/stage2.txt");
    // string str1, str2;
    // bool flag = true;
    // while (getline(file1, str1) && getline(file2, str2))
    // {
    //     if (str1 != str2)
    //     {
    //         flag = false;
    //         break;
    //     }
    // }
    // if (flag)
    //     cout << "Same" << endl;
    // else
    //     cout << "Not Same" << endl;

    return 0;
}
