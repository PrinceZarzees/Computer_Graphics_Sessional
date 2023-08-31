#include <GL/glut.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "ray_tracing.cpp"

GLdouble param = 5;
objects objs;
vector<light_source> ls_list;
double near_plane_d, far_plane_d, fovy, aspect_ratio;
int recursion_level, image_resolution;
Vector l;
Vector r;
Vector u;
point pos;
// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
std::vector<float> buildUnitPositiveX(int subdivision)
{
    const float DEG2RAD = acos(-1) / 180.0f;

    std::vector<float> vertices;
    float n1[3]; // normal of longitudinal plane rotating along Y-axis
    float n2[3]; // normal of latitudinal plane rotating along Z-axis
    float v[3];  // direction vector intersecting 2 planes, n1 x n2
    float a1;    // longitudinal angle along Y-axis
    float a2;    // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for (unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for (unsigned int j = 0; j < pointsPerRow; ++j)
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
            float scale = 1 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
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
void drawFace(const std::vector<float> &vertices, int pointsPerRow)
{
    glBegin(GL_TRIANGLES);

    for (int i = 0; i < pointsPerRow - 1; ++i)
    {
        for (int j = 0; j < pointsPerRow - 1; ++j)
        {
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
    for (int i = 0; i < vertices.size(); i += 3)
    {
        vertices[i] *= radius;
        vertices[i + 1] *= radius;
        vertices[i + 2] *= radius;
    }
    for (int i = 0; i < 4; i++)
    {
        glPushMatrix();
        glRotatef(90 * i, 0, 0, 1);
        drawFace(vertices, pointsPerRow);
        glPopMatrix();
    }
    for (int i = 1; i < 4; i += 2)
    {
        glPushMatrix();
        glRotatef(90 * i, 0, 1, 0);
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
    glVertex3f(width / 2, width / 2, 0);
    glVertex3f(0, 0, height);
    glVertex3f(width / 2, -width / 2, 0);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3f(-width / 2, width / 2, 0);
    glVertex3f(0, 0, height);
    glVertex3f(width / 2, width / 2, 0);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3f(-width / 2, -width / 2, 0);
    glVertex3f(0, 0, height);
    glVertex3f(-width / 2, width / 2, 0);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3f(width / 2, -width / 2, 0);
    glVertex3f(0, 0, height);
    glVertex3f(-width / 2, -width / 2, 0);
    glEnd();
    glBegin(GL_QUADS);
    glVertex3f(width / 2, width / 2, 0);
    glVertex3f(-width / 2, width / 2, 0);
    glVertex3f(-width / 2, -width / 2, 0);
    glVertex3f(width / 2, -width / 2, 0);
    glEnd();
}
void drawSquareFace(GLdouble width)
{
    glBegin(GL_QUADS);
    glVertex3f(0, 0, 0);
    glVertex3f(width, 0, 0);
    glVertex3f(width, width, 0);
    glVertex3f(0, width, 0);
    glEnd();
}
void drawCube(GLdouble side)
{
    glPushMatrix();
    glTranslatef(0, 0, side);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glRotatef(90, 1, 0, 0);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0, side, 0);
    glRotatef(90, 1, 0, 0);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glRotatef(-90, 0, 1, 0);
    drawSquareFace(side);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(side, 0, 0);
    glRotatef(-90, 0, 1, 0);
    drawSquareFace(side);
    glPopMatrix();
}
void reshapeListener(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    //cout<<width<<" "<<height<<endl;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();            // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    // gluPerspective(45.0f, aspect, 0.1f, 100.0f);
    gluPerspective(fovy, aspect, near_plane_d, far_plane_d);
}
/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y)
{
    GLdouble rate = 0.01;
    switch (key)
    {
    case '0':
        create_pointBuffer(near_plane_d, far_plane_d, fovy, aspect_ratio, l, r, u, pos, image_resolution,image_resolution,recursion_level, objs, ls_list);
        break;
    case '1':
        r.x = r.x * cos(rate) + l.x * sin(rate);
        r.y = r.y * cos(rate) + l.y * sin(rate);
        r.z = r.z * cos(rate) + l.z * sin(rate);

        l.x = l.x * cos(rate) - r.x * sin(rate);
        l.y = l.y * cos(rate) - r.y * sin(rate);
        l.z = l.z * cos(rate) - r.z * sin(rate);
        break;

    case '2':
        r.x = r.x * cos(-rate) + l.x * sin(-rate);
        r.y = r.y * cos(-rate) + l.y * sin(-rate);
        r.z = r.z * cos(-rate) + l.z * sin(-rate);

        l.x = l.x * cos(-rate) - r.x * sin(-rate);
        l.y = l.y * cos(-rate) - r.y * sin(-rate);
        l.z = l.z * cos(-rate) - r.z * sin(-rate);
        break;

    case '3':
        l.x = l.x * cos(rate) + u.x * sin(rate);
        l.y = l.y * cos(rate) + u.y * sin(rate);
        l.z = l.z * cos(rate) + u.z * sin(rate);

        u.x = u.x * cos(rate) - l.x * sin(rate);
        u.y = u.y * cos(rate) - l.y * sin(rate);
        u.z = u.z * cos(rate) - l.z * sin(rate);
        break;

    case '4':
        l.x = l.x * cos(-rate) + u.x * sin(-rate);
        l.y = l.y * cos(-rate) + u.y * sin(-rate);
        l.z = l.z * cos(-rate) + u.z * sin(-rate);

        u.x = u.x * cos(-rate) - l.x * sin(-rate);
        u.y = u.y * cos(-rate) - l.y * sin(-rate);
        u.z = u.z * cos(-rate) - l.z * sin(-rate);
        break;

    case '6':
        u.x = u.x * cos(rate) - r.x * sin(rate);
        u.y = u.y * cos(rate) - r.y * sin(rate);
        u.z = u.z * cos(rate) - r.z * sin(rate);

        r.x = r.x * cos(rate) + u.x * sin(rate);
        r.y = r.y * cos(rate) + u.y * sin(rate);
        r.z = r.z * cos(rate) + u.z * sin(rate);
        break;

    case '5':

        u.x = u.x * cos(-rate) - r.x * sin(-rate);
        u.y = u.y * cos(-rate) - r.y * sin(-rate);
        u.z = u.z * cos(-rate) - r.z * sin(-rate);

        r.x = r.x * cos(-rate) + u.x * sin(-rate);
        r.y = r.y * cos(-rate) + u.y * sin(-rate);
        r.z = r.z * cos(-rate) + u.z * sin(-rate);
        break;
    //pressing <space>
    case 32:
        texture=!texture;
        break;
    }
    glutPostRedisplay(); // Post a paint request to activate display()
}
/* Callback handler for special-key event */
void specialKeyListener(int key, int x, int y)
{

    switch (key)
    {


    case GLUT_KEY_DOWN: // down arrow key
        pos.x -= param * l.x;
        pos.y -= param * l.y;
        pos.z -= param * l.z;
        break;
    case GLUT_KEY_UP: // up arrow key
        pos.x += param * l.x;
        pos.y += param * l.y;
        pos.z += param * l.z;
        break;

    case GLUT_KEY_RIGHT:
        pos.x += param * r.x;
        pos.y += param * r.y;
        pos.z += param * r.z;
        break;
    case GLUT_KEY_LEFT:
        pos.x -= param * r.x;
        pos.y -= param * r.y;
        pos.z -= param * r.z;
        break;
    case GLUT_KEY_PAGE_UP:
        pos.x += param * u.x;
        pos.y += param * u.y;
        pos.z += param * u.z;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos.x -= param * u.x;
        pos.y -= param * u.y;
        pos.z -= param * u.z;
        break;

    default:
        break;
    }

    glutPostRedisplay(); // Post a paint request to activate display()
}

int process_input(ifstream &infile)
{
    if (!infile)
    {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }
    double checkerboard_width, checkerboard_ambient, checkerboard_diffuse, checkerboard_reflection;
    int num_objects, num_lights, num_spot_lights;

    infile >> near_plane_d >> far_plane_d >> fovy >> aspect_ratio;
    infile >> recursion_level;
    infile >> image_resolution;
    objs = objects();
    infile >> checkerboard_width;
    infile >> checkerboard_ambient >> checkerboard_diffuse >> checkerboard_reflection;
    objs.add_checkerboard(checkerboard(checkerboard_width));
    objs.cb.set_color_component(checkerboard_ambient, checkerboard_diffuse, checkerboard_reflection);
    infile >> num_objects;
    for (int i = 0; i < num_objects; ++i)
    {
        string type;
        infile >> type;
        if (type == "sphere")
        {
            double c_x, c_y, c_z, radius, color_r, color_g, color_b, ambient, diffuse, specular, reflection, shininess;
            infile >> c_x >> c_y >> c_z;
            infile >> radius;
            infile >> color_r >> color_g >> color_b;
            infile >> ambient >> diffuse >> specular >> reflection;
            infile >> shininess;
            objs.add_sphere(Sphere(point(c_x, c_y, c_z), radius));
            objs.spheres[objs.spheres.size() - 1].set_color(color_r, color_g, color_b);
            objs.spheres[objs.spheres.size() - 1].set_color_component(ambient, diffuse, specular, reflection, shininess);
        }
        if (type == "pyramid")
        {
            double l_x, l_y, l_z, width, height, color_r, color_g, color_b, ambient, diffuse, specular, reflection, shininess;
            infile >> l_x >> l_y >> l_z;
            infile >> width >> height;
            infile >> color_r >> color_g >> color_b;
            infile >> ambient >> diffuse >> specular >> reflection;
            infile >> shininess;
            objs.add_pyramid(pyramid(point(l_x+width/2, l_y+width/2, l_z), width, height));
            objs.pyramids[objs.pyramids.size() - 1].set_color(color_r, color_g, color_b);
            objs.pyramids[objs.pyramids.size() - 1].set_color_component(ambient, diffuse, specular, reflection, shininess);
        }
        if (type == "cube")
        {
            double b_x, b_y, b_z, side, color_r, color_g, color_b, ambient, diffuse, specular, reflection, shininess;
            infile >> b_x >> b_y >> b_z;
            infile >> side;
            infile >> color_r >> color_g >> color_b;
            infile >> ambient >> diffuse >> specular >> reflection;
            infile >> shininess;
            objs.add_cube(cube(point(b_x, b_y, b_z), side));
            objs.cubes[objs.cubes.size() - 1].set_color(color_r, color_g, color_b);
            objs.cubes[objs.cubes.size() - 1].set_color_component(ambient, diffuse, specular, reflection, shininess);
        }
    }

    infile >> num_lights;
    for (int i = 0; i < num_lights; ++i)
    {
        double x, y, z, falloff;
        infile >> x >> y >> z >> falloff;
        ls_list.push_back(light_source(point(x, y, z), point(0, 0, 0), -1, falloff, true));
    }

    infile >> num_spot_lights;
    for (int i = 0; i < num_spot_lights; ++i)
    {
        double x, y, z, falloff, d_x, d_y, d_z, cutoff_angle;
        infile >> x >> y >> z >> falloff;
        infile >> d_x >> d_y >> d_z;
        infile >> cutoff_angle;
        cutoff_angle = cutoff_angle * M_PI / 180.0;

        ls_list.push_back(light_source(point(x, y, z), point(d_x, d_y, d_z), cutoff_angle, falloff, false));
    }

    // Now you have parsed the data and stored it in appropriate structures
    // You can use the data for further processing

    return 0;
}
// Post a paint request to activate display()
void displayFunc()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(pos.x, pos.y, pos.z,
              pos.x + l.x, pos.y + l.y, pos.z + l.z,
              u.x, u.y, u.z);
    // draw axis
    glLineWidth(3);
    glBegin(GL_LINES);
    glColor3f(0, 1, 0);
    glVertex3f(pos.x - 1000, 0, 0);
    glVertex3f(pos.x + 1000, 0, 0);
    glColor3f(0, 1, 0);
    glVertex3f(0, pos.y - 1000, 0);
    glVertex3f(0, pos.y + 1000, 0);
    glColor3f(0, 1, 0);
    glVertex3f(0, 0, pos.z - 1000);
    glVertex3f(0, 0, pos.z + 1000);
    glEnd();

    // draw a checkerboard with certain width of each cell
    glColor3f(0, 0, 0);
    //GLdouble checkerboard_width = 50;
    // drawSquareFace(checkerboard_width);
    int checker_x = (pos.x / objs.cb.side);
    int checker_y = (pos.y / objs.cb.side);
    for (int i = checker_x - 100; i <= checker_x + 100; i++)
    {
        for (int j = checker_y - 100; j <= checker_y + 100; j++)
        {
            glPushMatrix();
            glColor3f(abs(i + j + 1) % 2, abs(i + j + 1) % 2, abs(i + j + 1) % 2);
            glTranslatef(i * objs.cb.side, j * objs.cb.side, 0);
            drawSquareFace(objs.cb.side);
            // std::cout<<i*checkerboard_width<<","<<j*checkerboard_width<<" "<<pos.x<<","<<pos.y<<std::endl;
            glPopMatrix();
        }
    }

    for(int i=0;i<objs.spheres.size();i++)
    {
         // draw a sphere
    glColor3f(objs.spheres[i].red, objs.spheres[i].green, objs.spheres[i].blue);
    glPushMatrix();
    glTranslatef(objs.spheres[i].center.x, objs.spheres[i].center.y, objs.spheres[i].center.z);
    drawSphere(objs.spheres[i].radius);
    glPopMatrix();
    }

    for(int i=0;i<objs.pyramids.size();i++)
    {
    // draw a pyramid
    glColor3f(objs.pyramids[i].red, objs.pyramids[i].green, objs.pyramids[i].blue);
    glPushMatrix();
    glTranslatef(objs.pyramids[i].bottom_center.x, objs.pyramids[i].bottom_center.y, objs.pyramids[i].bottom_center.z);
    drawPyramid(objs.pyramids[i].side, objs.pyramids[i].height);
    glPopMatrix();
    }

    for(int i=0;i<objs.cubes.size();i++)
    {
    glColor3f(objs.cubes[i].red, objs.cubes[i].green, objs.cubes[i].blue);
    glPushMatrix();
    glTranslatef(objs.cubes[i].bottom_lower_left.x, objs.cubes[i].bottom_lower_left.y, objs.cubes[i].bottom_lower_left.z);
    drawCube(objs.cubes[i].side);
    glPopMatrix();
    }

   

    // draw a sphere
    // glColor3f(1.0f, 0.0f, 1.0f);
    // radius = 15;
    // sphere_x = -20, sphere_y = -20, sphere_z = 20;
    // glPushMatrix();
    // glTranslatef(sphere_x, sphere_y, sphere_z);
    // drawSphere(radius);
    // glPopMatrix();

    // draw a pyramid
    // GLdouble pyramid_width = 30, pyramid_height = 40;
    // GLdouble pyramid_x = -40, pyramid_y = 0, pyramid_z = 5;
    // glColor3f(1.0f, 0.0f, 0.0f);
    // glPushMatrix();
    // glTranslatef(pyramid_x + pyramid_width / 2, pyramid_y + pyramid_width / 2, pyramid_z);
    // drawPyramid(pyramid_width, pyramid_height);
    // glPopMatrix();

    // draw a cube
    // GLdouble cube_side = 40;
    // GLdouble cube_x = -100, cube_y = -100, cube_z = 10;
    // glColor3f(0.0f, 0.5f, 1.0f);
    // glPushMatrix();
    // glTranslatef(cube_x, cube_y, cube_z);
    // drawCube(cube_side);
    // glPopMatrix();

    glutSwapBuffers();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    ifstream infile("description.txt");
    process_input(infile);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(700,700);
    glutCreateWindow("Magic Cube");
    pos.x = 0;
    pos.y = -150;
    pos.z = 50;

    l.x = 0;
    l.y = 1;
    l.z = 0;
    u.x = 0;
    u.y = 0;
    u.z = 1;
    r.x = 1;
    r.y = 0;
    r.z = 0;
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(displayFunc);
    glutReshapeFunc(reshapeListener);   // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener); // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();
    return 0;
}
