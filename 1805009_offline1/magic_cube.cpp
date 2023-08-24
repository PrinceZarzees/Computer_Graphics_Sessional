#include <GL/glut.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "ray_tracing.cpp"
//GLdouble scale=1;
//GLdouble rotateAngle=0;
GLdouble param=0.5;
// struct point {
//     GLfloat x, y, z;
// };
// struct point pos;   // position of the eye
// struct point l;     // look/forward direction
// struct point r;     // right direction
// struct point u;     // up direction
Vector l;
Vector r;
Vector u;
point pos;
int i_width;
int i_height;
double aspect_ratio;
double fovy;
double znear;
double zfar;
// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
std::vector<float> buildUnitPositiveX(int subdivision)
{
    const float DEG2RAD = acos(-1) / 180.0f;

    std::vector<float> vertices;
    float n1[3];        // normal of longitudinal plane rotating along Y-axis
    float n2[3];        // normal of latitudinal plane rotating along Z-axis
    float v[3];         // direction vector intersecting 2 planes, n1 x n2
    float a1;           // longitudinal angle along Y-axis
    float a2;           // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scale = 1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;

            // add a vertex into array
            vertices.push_back(v[0]);
            vertices.push_back(v[1]);
            vertices.push_back(v[2]);
        }
    }
    return vertices;
}
// void drawCylindricalFace(GLdouble angle, GLdouble height,GLdouble radius,GLdouble segments)
// {
//     double tempx = 0, tempy = radius;
//     double currx, curry;
//     glBegin(GL_QUADS);
//         for (int i = 1; i <= segments; i++) {
//             double theta = i * angle / (2*segments);
//             currx = radius * sin(theta);
//             curry = radius * cos(theta);
//             glVertex3f(currx, curry, height/2);
//             glVertex3f(currx, curry, -height/2);

//             glVertex3f(tempx, tempy, -height/2);
//             glVertex3f(tempx, tempy, height/2);

//             tempx = currx;
//             tempy = curry;
//         }
//         tempx = 0;
//         tempy = radius;
//          for (int i = 1; i <= segments; i++) {
//             double theta = -i * angle / (2*segments);
//             currx = radius * sin(theta);
//             curry = radius * cos(theta);
//             glVertex3f(currx, curry, height/2);
//             glVertex3f(currx, curry, -height/2);

//             glVertex3f(tempx, tempy, -height/2);
//             glVertex3f(tempx, tempy, height/2);

//             tempx = currx;
//             tempy = curry;
//         }
//     glEnd();
// }
void drawFace(const std::vector<float>& vertices, int pointsPerRow)
{
    glBegin(GL_TRIANGLES);

    for (int i = 0; i < pointsPerRow - 1; ++i) {
        for (int j = 0; j < pointsPerRow - 1; ++j) {
            // Compute indices for the current quad
            int idx1 = i * pointsPerRow + j;
            int idx2 = (i + 1) * pointsPerRow + j;
            int idx3 = i * pointsPerRow + (j + 1);
            int idx4 = (i + 1) * pointsPerRow + (j + 1);

            // Draw the first triangle
            glVertex3f(vertices[idx1 * 3], vertices[idx1 * 3 + 1], vertices[idx1 * 3 + 2]);
            glVertex3f(vertices[idx2 * 3], vertices[idx2 * 3 + 1], vertices[idx2 * 3 + 2]);
            glVertex3f(vertices[idx3 * 3], vertices[idx3 * 3 + 1], vertices[idx3 * 3 + 2]);

            // Draw the second triangle
            glVertex3f(vertices[idx3 * 3], vertices[idx3 * 3 + 1], vertices[idx3 * 3 + 2]);
            glVertex3f(vertices[idx2 * 3], vertices[idx2 * 3 + 1], vertices[idx2 * 3 + 2]);
            glVertex3f(vertices[idx4 * 3], vertices[idx4 * 3 + 1], vertices[idx4 * 3 + 2]);
        }
    }

    glEnd();
}
void drawSphere(GLdouble radius)
{
    int subdivisionLevel = 3; // Example subdivision level
std::vector<float> vertices = buildUnitPositiveX(subdivisionLevel);
int pointsPerRow = std::pow(2, subdivisionLevel) + 1;
// Scale the vertices by the radius
for (int i = 0; i < vertices.size(); i += 3) {
    vertices[i] *= radius;
    vertices[i + 1] *= radius;
    vertices[i + 2] *= radius;
}
for(int i=0;i<4;i++)
{
glPushMatrix();
glRotatef(90*i,0,0,1);
drawFace(vertices, pointsPerRow);
glPopMatrix();
}
for(int i=1;i<4;i+=2)
{
glPushMatrix();
glRotatef(90*i,0,1,0);
drawFace(vertices, pointsPerRow);
glPopMatrix();
}
}
// void drawTraingularFace()
// {
//     glBegin(GL_TRIANGLES);
//     glVertex3f(0,0,scale);
//     glVertex3f(0,scale,0);
//     glVertex3f(scale,0,0);
//     glEnd();
// }
void drawPyramid(GLdouble width, GLdouble height)
{
    glBegin(GL_TRIANGLES);
    glVertex3f(width/2,0,width/2);
    glVertex3f(0,height,0);
    glVertex3f(width/2,0,-width/2);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3f(-width/2,0,width/2);
    glVertex3f(0,height,0);
    glVertex3f(width/2,0,width/2);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3f(-width/2,0,-width/2);
    glVertex3f(0,height,0);
    glVertex3f(-width/2,0,width/2);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3f(width/2,0,-width/2);
    glVertex3f(0,height,0);
    glVertex3f(-width/2,0,-width/2);
    glEnd();
    glBegin(GL_QUADS);
    glVertex3f(width/2,0,width/2);
    glVertex3f(-width/2,0,width/2);
    glVertex3f(-width/2,0,-width/2);
    glVertex3f(width/2,0,-width/2);
    glEnd();
}
void drawSquareFace(GLdouble width)
{
    glBegin(GL_QUADS);
    glVertex3f(0,0,0);
    glVertex3f(width,0,0);
    glVertex3f(width,0,width);
    glVertex3f(0,0,width);
    glEnd();
}
void drawCube(GLdouble side)
{
    glPushMatrix();
    glTranslatef(0,side,0);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glRotatef(-90,1,0,0);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0,0,side);
    glRotatef(-90,1,0,0);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glRotatef(90,0,0,1);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(side,0,0);
    glRotatef(90,0,0,1);
    drawSquareFace(side);
    glPopMatrix();
}
void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    // gluPerspective(45.0f, aspect, 0.1f, 100.0f);
    i_width=width;
    i_height=height;
    aspect_ratio=aspect;
    znear=1;
    zfar=1000;
    fovy=80;
    gluPerspective(80.0f, aspect, 1.0f, 1000.0f);
}
/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y) {
   GLdouble rate = 0.01;
   //GLdouble tempx,tempy,tempz,t;
    switch (key) {
    // case ',':
    //     scale-=0.1;
    //     if (scale<=0)
    //         scale=0;
    //     break;
    // case '.':
    //     scale+=0.1;
    //     if (scale>=1)
    //         scale=1;
    //     break;
    // case 'a':
    //     rotateAngle-=2;
    //     break;
    // case 'd':
    //     rotateAngle+=2;
    //     break;
    case '0':
            create_pointBuffer(znear,zfar,fovy,aspect_ratio,l,r,u,pos,i_width,i_height);
            break;
    case '1':
			r.x = r.x*cos(rate)+l.x*sin(rate);
			r.y = r.y*cos(rate)+l.y*sin(rate);
			r.z = r.z*cos(rate)+l.z*sin(rate);

			l.x = l.x*cos(rate)-r.x*sin(rate);
			l.y = l.y*cos(rate)-r.y*sin(rate);
			l.z = l.z*cos(rate)-r.z*sin(rate);
			break;

        case '2':
			r.x = r.x*cos(-rate)+l.x*sin(-rate);
			r.y = r.y*cos(-rate)+l.y*sin(-rate);
			r.z = r.z*cos(-rate)+l.z*sin(-rate);

			l.x = l.x*cos(-rate)-r.x*sin(-rate);
			l.y = l.y*cos(-rate)-r.y*sin(-rate);
			l.z = l.z*cos(-rate)-r.z*sin(-rate);
			break;

        case '3':
			l.x = l.x*cos(rate)+u.x*sin(rate);
			l.y = l.y*cos(rate)+u.y*sin(rate);
			l.z = l.z*cos(rate)+u.z*sin(rate);

			u.x = u.x*cos(rate)-l.x*sin(rate);
			u.y = u.y*cos(rate)-l.y*sin(rate);
			u.z = u.z*cos(rate)-l.z*sin(rate);
			break;

        case '4':
			l.x = l.x*cos(-rate)+u.x*sin(-rate);
			l.y = l.y*cos(-rate)+u.y*sin(-rate);
			l.z = l.z*cos(-rate)+u.z*sin(-rate);

			u.x = u.x*cos(-rate)-l.x*sin(-rate);
			u.y = u.y*cos(-rate)-l.y*sin(-rate);
			u.z = u.z*cos(-rate)-l.z*sin(-rate);
			break;

        case '6':
			u.x = u.x*cos(rate)-r.x*sin(rate);
			u.y = u.y*cos(rate)-r.y*sin(rate);
			u.z = u.z*cos(rate)-r.z*sin(rate);

			r.x = r.x*cos(rate)+u.x*sin(rate);
			r.y = r.y*cos(rate)+u.y*sin(rate);
			r.z = r.z*cos(rate)+u.z*sin(rate);
			break;

        case '5':

			u.x = u.x*cos(-rate)-r.x*sin(-rate);
			u.y = u.y*cos(-rate)-r.y*sin(-rate);
			u.z = u.z*cos(-rate)-r.z*sin(-rate);

			r.x = r.x*cos(-rate)+u.x*sin(-rate);
			r.y = r.y*cos(-rate)+u.y*sin(-rate);
			r.z = r.z*cos(-rate)+u.z*sin(-rate);
			break;
        // case 'w':
        // a.x=l.z;a.y=0;a.z=-l.x;
        // a.x=a.x/sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        // a.y=a.y/sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        // a.z=a.z/sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        // tempx=cos(-rate)*l.x+(1-cos(-rate))*(a.x*a.x*l.x+a.x*a.y*l.y+a.x*a.z*l.z)+sin(-rate)*(a.y*l.z-a.z*l.y);
        // tempy=cos(-rate)*l.y+(1-cos(-rate))*(a.x*a.y*l.x+a.y*a.y*l.y+a.y*a.z*l.z)+sin(-rate)*(a.z*l.x-a.x*l.z);
        // tempz=cos(-rate)*l.z+(1-cos(-rate))*(a.x*a.z*l.x+a.y*a.z*l.y+a.z*a.z*l.z)+sin(-rate)*(a.x*l.y-a.y*l.x);
        // l.x=tempx;l.y=tempy;l.z=tempz;
        // tempx=cos(-rate)*u.x+(1-cos(-rate))*(a.x*a.x*u.x+a.x*a.y*u.y+a.x*a.z*u.z)+sin(-rate)*(a.y*u.z-a.z*u.y);
        // tempy=cos(-rate)*u.y+(1-cos(-rate))*(a.x*a.y*u.x+a.y*a.y*u.y+a.y*a.z*u.z)+sin(-rate)*(a.z*u.x-a.x*u.z);
        // tempz=cos(-rate)*u.z+(1-cos(-rate))*(a.x*a.z*u.x+a.y*a.z*u.y+a.z*a.z*u.z)+sin(-rate)*(a.x*u.y-a.y*u.x);
        // u.x=tempx;u.y=tempy;u.z=tempz;
        // tempx=cos(-rate)*r.x+(1-cos(-rate))*(a.x*a.x*r.x+a.x*a.y*r.y+a.x*a.z*r.z)+sin(-rate)*(a.y*r.z-a.z*r.y);
        // tempy=cos(-rate)*r.y+(1-cos(-rate))*(a.x*a.y*r.x+a.y*a.y*r.y+a.y*a.z*r.z)+sin(-rate)*(a.z*r.x-a.x*r.z);
        // tempz=cos(-rate)*r.z+(1-cos(-rate))*(a.x*a.z*r.x+a.y*a.z*r.y+a.z*a.z*r.z)+sin(-rate)*(a.x*r.y-a.y*r.x);
        // r.x=tempx;r.y=tempy;r.z=tempz;
        // // move up without changing reference point
        // t=sqrt(pos.x*pos.x+pos.z*pos.z)*tan(atan(pos.y/sqrt(pos.x*pos.x+pos.z*pos.z))+rate)-pos.y;
        // pos.y+=t;
        // break;

        // case 's':
        // a.x=l.z;a.y=0;a.z=-l.x;
        // a.x=a.x/sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        // a.y=a.y/sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        // a.z=a.z/sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        // tempx=cos(rate)*l.x+(1-cos(rate))*(a.x*a.x*l.x+a.x*a.y*l.y+a.x*a.z*l.z)+sin(rate)*(a.y*l.z-a.z*l.y);
        // tempy=cos(rate)*l.y+(1-cos(rate))*(a.x*a.y*l.x+a.y*a.y*l.y+a.y*a.z*l.z)+sin(rate)*(a.z*l.x-a.x*l.z);
        // tempz=cos(rate)*l.z+(1-cos(rate))*(a.x*a.z*l.x+a.y*a.z*l.y+a.z*a.z*l.z)+sin(rate)*(a.x*l.y-a.y*l.x);
        // l.x=tempx;l.y=tempy;l.z=tempz;
        // tempx=cos(rate)*u.x+(1-cos(rate))*(a.x*a.x*u.x+a.x*a.y*u.y+a.x*a.z*u.z)+sin(rate)*(a.y*u.z-a.z*u.y);
        // tempy=cos(rate)*u.y+(1-cos(rate))*(a.x*a.y*u.x+a.y*a.y*u.y+a.y*a.z*u.z)+sin(rate)*(a.z*u.x-a.x*u.z);
        // tempz=cos(rate)*u.z+(1-cos(rate))*(a.x*a.z*u.x+a.y*a.z*u.y+a.z*a.z*u.z)+sin(rate)*(a.x*u.y-a.y*u.x);
        // u.x=tempx;u.y=tempy;u.z=tempz;
        // tempx=cos(rate)*r.x+(1-cos(rate))*(a.x*a.x*r.x+a.x*a.y*r.y+a.x*a.z*r.z)+sin(rate)*(a.y*r.z-a.z*r.y);
        // tempy=cos(rate)*r.y+(1-cos(rate))*(a.x*a.y*r.x+a.y*a.y*r.y+a.y*a.z*r.z)+sin(rate)*(a.z*r.x-a.x*r.z);
        // tempz=cos(rate)*r.z+(1-cos(rate))*(a.x*a.z*r.x+a.y*a.z*r.y+a.z*a.z*r.z)+sin(rate)*(a.x*r.y-a.y*r.x);
        // r.x=tempx;r.y=tempy;r.z=tempz;
        // // move down without changing reference point
        // t=-sqrt(pos.x*pos.x+pos.z*pos.z)*tan(atan(pos.y/sqrt(pos.x*pos.x+pos.z*pos.z))-rate)+pos.y;
        // pos.y-=t;      
        //     break;

    }
    glutPostRedisplay();    // Post a paint request to activate display()
}
/* Callback handler for special-key event */
void specialKeyListener(int key, int x,int y) {

    switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x-=param*l.x;
			pos.y-=param*l.y;
			pos.z-=param*l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
			pos.x+=param*l.x;
			pos.y+=param*l.y;
			pos.z+=param*l.z;
			break;

		case GLUT_KEY_RIGHT:
			pos.x+=param*r.x;
			pos.y+=param*r.y;
			pos.z+=param*r.z;
			break;
		case GLUT_KEY_LEFT :
			pos.x-=param*r.x;
			pos.y-=param*r.y;
			pos.z-=param*r.z;
			break;
        case GLUT_KEY_PAGE_UP:
		    pos.x+=param*u.x;
			pos.y+=param*u.y;
			pos.z+=param*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
            pos.x-=param*u.x;
			pos.y-=param*u.y;
			pos.z-=param*u.z;
			break;

		default:
			break;
    }

    glutPostRedisplay(); // Post a paint request to activate display()
}
  // Post a paint request to activate display()
void displayFunc() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
gluLookAt(pos.x,pos.y,pos.z,
              pos.x+l.x,pos.y+l.y,pos.z+l.z,
              u.x,u.y,u.z);
// draw axis
glLineWidth(3);
glBegin(GL_LINES);
glColor3f(0,1,0);
glVertex3f(pos.x-100,0,0);
glVertex3f(pos.x+100,0,0);
glColor3f(0,1,0);
glVertex3f(0,pos.y-100,0);
glVertex3f(0,pos.y+100,0);
glColor3f(0,1,0);
glVertex3f(0,0,pos.z-100);
glVertex3f(0,0,pos.z+100);
glEnd();

// draw a checkerboard with certain width of each cell
glColor3f(0,0,0);
GLdouble checkerboard_width = 1;
drawSquareFace(checkerboard_width);
int checker_x = (pos.x / checkerboard_width)*checkerboard_width;
int checker_z = (pos.z / checkerboard_width)*checkerboard_width;
for (int i=checker_x-100;i<=checker_x+100;i++)
{
    for (int j=checker_z-100;j<=checker_z+100;j++)
    {
            glPushMatrix();
            glColor3f(abs(i+j)%2,abs(i+j)%2,abs(i+j)%2);
            glTranslatef(i*checkerboard_width,0,j*checkerboard_width);
            drawSquareFace(checkerboard_width);
            glPopMatrix();
        
    }
}

//draw a sphere
glColor3f(1.0f, 0.0f, 0.0f);
GLdouble radius = 1;
GLdouble sphere_x=-1,sphere_y=1,sphere_z=0;
glPushMatrix();
glTranslatef(sphere_x,sphere_y,sphere_z);
drawSphere(radius);
glPopMatrix();

//draw a pyramid
GLdouble pyramid_width=1,pyramid_height=1;
GLdouble pyramid_x=1,pyramid_y=1,pyramid_z=0;
glColor3f(1.0f, 1.0f, 0.5f);
glPushMatrix();
glTranslatef(pyramid_x,pyramid_y,pyramid_z);
drawPyramid(pyramid_width,pyramid_height);
glPopMatrix();

//draw a cube
GLdouble cube_side=1;
GLdouble cube_x=0,cube_y=1,cube_z=1;
glColor3f(0.0f, 1.0f, 1.0f);
glPushMatrix();
glTranslatef(cube_x,cube_y,cube_z);
drawCube(cube_side);
glPopMatrix();





// for(int i=0;i<4;i++)
// {
//     glPushMatrix();
//     if (i%2==0)
//     glColor3f(0,1,1);
//     else
//     glColor3f(1,0,1);
//     glRotatef(90*i,0,0,1);
//     glTranslatef(1/3.0-scale/3.0,1/3.0-scale/3.0,1/3.0-scale/3.0);
//     drawTraingularFace();
//     glPopMatrix();
// }

// for(int i=0;i<4;i++)
// {
// glPushMatrix();
// if  (i%2==0)
//     glColor3f(1,0,1);
// else
//     glColor3f(0,1,1);
// glRotatef(90,0,1,0);
// glRotatef(90*i,1,0,0);
// glTranslatef(1/3.0-scale/3.0,1/3.0-scale/3.0,1/3.0-scale/3.0);
// drawTraingularFace();
// glPopMatrix();

// }

// GLdouble angle=acos(1/3.0);
// GLdouble height=sqrt(2)*scale;
// GLdouble radius=1/sqrt(3)*(1-scale);
// GLdouble d=1/sqrt(2)*scale;
// for(int i=0;i<4;i++)
// {glPushMatrix();
// glColor3f(1,1,0);
// glRotatef(90,0,1,0);
// glRotatef(45+i*90,1,0,0);
// glTranslatef(0,d,0);
// drawCylindricalFace(angle,height,radius,50);
// glPopMatrix();
// }
// for(int i=0;i<4;i++)
// {glPushMatrix();
// glColor3f(1,1,0);
// glRotatef(45+i*90,1,0,0);
// glTranslatef(0,d,0);
// drawCylindricalFace(angle,height,radius,50);
// glPopMatrix();
// }
// for(int i=0;i<4;i++)
// {glPushMatrix();
// glColor3f(1,1,0);
// glRotatef(-90,0,0,1);
// glRotatef(45+i*90,1,0,0);
// glTranslatef(0,d,0);
// drawCylindricalFace(angle,height,radius,50);
// glPopMatrix();
// }


//     // Set color to blue

    glutSwapBuffers();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(768, 768);
    glutCreateWindow("Magic Cube");
     pos.x=0;pos.y=1;pos.z=5;
    
    l.x=0;l.y=0;l.z=-1;
    u.x=0;u.y=1;u.z=0;
    r.x=1;r.y=0;r.z=0;
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
     
    glutDisplayFunc(displayFunc);
     glutReshapeFunc(reshapeListener);           // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);         // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();
    return 0;
}
