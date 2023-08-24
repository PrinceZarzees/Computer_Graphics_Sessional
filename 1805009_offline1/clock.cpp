/*
 * GL04ModelTransform.cpp: Model Transform - Translation and Rotation
 * Transform primitives from their model spaces to world space.
 */
#define USE_MATH_DEFINES
#include <windows.h> // for MS Windows
#include <GL/glut.h> // GLUT, include glu.h and gl.h
#include <math.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

std::time_t currentTime;

string timeString;

int hours, minutes, seconds;
GLdouble hangle, mangle, sangle, pangle, pangle_max;
GLdouble p_t = 0;
GLdouble p_l=0.5;
GLdouble innerRadius=0.40;
GLdouble outerRadius=0.525;
/* Initialize OpenGL Graphics */
void initGL()
{
   // Set "clearing" or background color
   pangle_max = 16;
   pangle = pangle_max;
   currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
   timeString = ctime(&currentTime);
   sscanf(timeString.c_str(), "%*s %*s %*d %d:%d:%d", &hours, &minutes, &seconds);
   if (hours == 0)
      hours = 12;
   else if (hours > 12)
      hours = hours - 12;
   sangle = seconds * 6;
   mangle = minutes * 6 + seconds * 0.1;
   hangle = hours * 30 + minutes * 0.5 + seconds * (1 / 120.0);
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

void drawRectangle(GLdouble height, GLdouble width)
{
   glBegin(GL_QUADS);
   glVertex2f(-width / 2, -height / 2);
   glVertex2f(width / 2, -height / 2);
   glVertex2f(width / 2, height / 2);
   glVertex2f(-width / 2, height / 2);
   glEnd();
}

void drawTriangle()
{
   glBegin(GL_TRIANGLES);
   glVertex2f(-0.05, 0.1);
   glVertex2f(0.05, 0.1);
   glVertex2f(0, 0.15);
   glEnd();
}
void drawLine(GLdouble length)
{
   glBegin(GL_LINES);
   glVertex2f(0, 0);
   glVertex2f(0, length);
   glEnd();
}
void drawhorizontalline(GLdouble length)
{
   glBegin(GL_LINES);
   glVertex2f(-length/2, 0);
   glVertex2f(length/2, 0);
   glEnd();

}
void drawClockHand(GLdouble length, GLdouble width)
{
   glBegin(GL_POLYGON);
   glVertex2f(-width / 2, length / 3);
   glVertex2f(0, length);
   glVertex2f(width / 2, length / 3);
   glVertex2f(0, 0);
   glEnd();
}
void drawCircle(GLdouble radius)
{
   glBegin(GL_LINE_LOOP);
   for (GLdouble i = 0; i < 360; i += 0.3)
   {
      glVertex2f(radius * cos(i * M_PI / 180.0), radius * sin(i * M_PI / 180.0));
   }
   glEnd();
}
void drawCircleFilled(GLdouble radius)
{
   glBegin(GL_POLYGON);
   for (GLdouble i = 0; i < 360; i += 0.5)
   {
      glVertex2f(radius * cos(i * M_PI / 180.0), radius * sin(i * M_PI / 180.0));
   }
   glEnd();
}
void drawClock()
{
   drawRectangle(outerRadius/14.0, outerRadius/14.0);
   drawCircle(outerRadius);
   glColor3f(0,0.5,0.5);
   drawCircle(innerRadius);
   for (int i = 0; i < 360; i += 30)
   {
      glPushMatrix();
      glRotatef(i, 0, 0, 1);
      if (i % 90 == 0)
      {
         glTranslatef(0, innerRadius-0.075, 0);
         drawLine(0.075);
      }
      else
      {
         glTranslatef(0,innerRadius-0.0375, 0);
         drawLine(0.0375);
      }
      glPopMatrix();
   }
   glPushMatrix();
   glRotatef(-hangle, 0, 0, 1);
   glColor3f(0,0,1);
   drawClockHand(innerRadius/2, innerRadius/20);
   glPopMatrix();
   glPushMatrix();
   glRotatef(-mangle, 0, 0, 1);
   glColor3f(0,1,0);
   drawClockHand(innerRadius*7/12.0, innerRadius/30);
   glPopMatrix();
   glPushMatrix();
   glRotatef(-sangle, 0, 0, 1);
   glColor3f(1,0,0);
   drawClockHand(innerRadius*2/3.0, innerRadius/40);
   glPopMatrix();
   glPushMatrix();
   glBegin(GL_LINE_LOOP);
   glColor3f(1.0f, 1.0f, 0.0f);
   glVertex2f(-outerRadius-0.1,outerRadius+0.1);
   glColor3f(1.0f, 0.5f, 0.0f);
   glVertex2f(outerRadius+0.1,outerRadius+0.1);
   glColor3f(0.5f, 0.5f, 0.0f);
   glVertex2f(outerRadius+0.1,-outerRadius-(p_l+0.1)/cos(pangle_max*M_PI/180.0)+(outerRadius+0.1)*tan(pangle_max*M_PI/180.0));
   glVertex2f(0,-outerRadius-(p_l+0.1)/cos(pangle_max*M_PI/180.0));
   glVertex2f(-outerRadius-0.1,-outerRadius-(p_l+0.1)/cos(pangle_max*M_PI/180.0)+(outerRadius+0.1)*tan(pangle_max*M_PI/180.0));
   glEnd();
   glPopMatrix();
}
void drawPendulum()
{

   glPushMatrix();
   glRotatef(pangle, 0, 0, 1);
   glTranslatef(0, -p_l, 0);
   glColor3f(0,0.6,0.7);
   drawCircleFilled(0.1);
   glPopMatrix();
   glPushMatrix();
   glRotatef(pangle, 0, 0, 1);
   glTranslatef(0, -p_l/2.0, 0);
   glColor3f(0,0.6,0.7);
   drawRectangle(p_l, 0.05);
   glPopMatrix();
}
/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display()
{
   glClear(GL_COLOR_BUFFER_BIT); // Clear the color buffer
   glMatrixMode(GL_MODELVIEW);   // To operate on Model-View matrix
   glLoadIdentity();             // Reset the model-view matrix

   glPushMatrix();
   glTranslatef(0, 0.3, 0);
   drawClock();
   glPopMatrix();
   glPushMatrix();
   glTranslatef(0, -outerRadius+0.3, 0);
   drawPendulum();
   glPopMatrix();
   glPushMatrix();
   glTranslatef(0,  -outerRadius+0.3, 0);
   drawCircleFilled(0.05);
   glPopMatrix();
   // glPushMatrix();
   // glTranslatef(0, -outerRadius+0.3, 0);
   // glRotatef(pangle_max,0,0,1);
   // glTranslatef(0, -p_l-0.1, 0);
   //  glBegin(GL_LINES);
   // glVertex2f((-p_l-0.1)*sin(pangle_max*M_PI/180.0)-0.01, 0);
   // glVertex2f(0.3, 0);
   // glEnd();
   // glPopMatrix();
   // glPushMatrix();
   // glTranslatef(0, -outerRadius+0.3, 0);
   // glRotatef(-pangle_max,0,0,1);
   // glTranslatef(0, -p_l-0.1, 0);
   // glBegin(GL_LINES);
   // glVertex2f((p_l+0.1)*sin(pangle_max*M_PI/180.0)+0.01, 0);
   // glVertex2f(-0.3, 0);
   // glEnd();
   // glPopMatrix();
   glFlush();
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
   // Compute aspect ratio of the new window
   if (height == 0)
      height = 1; // To prevent divide by 0
   GLfloat aspect = (GLfloat)width / (GLfloat)height;

   // Set the viewport to cover the new window
   glViewport(0, 0, width, height);

   // Set the aspect ratio of the clipping area to match the viewport
   glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
   glLoadIdentity();
   if (width >= height)
   {
      // aspect >= 1, set the height from -1 to 1, with larger width
      gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
   }
   else
   {
      // aspect < 1, set the width to -1 to 1, with larger height
      gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
   }
}
void update(int value)
{
   sangle += 6;
   mangle += 0.1;
   hangle += (1 / 120.0);
   if (sangle>=360)
   {
      sangle-=360;
   }
   if (mangle>=360)
      mangle-=360;
   if (hangle>=360)
      hangle-=360;
   glutPostRedisplay();
   glutTimerFunc(1000, update, 0);
}
void update_p(int value)
{

   pangle = pangle_max * cos((2 * M_PI / 2.0) * p_t);
   p_t += 0.01;
   if (p_t==2)
   {
      p_t=0;
   }
   glutPostRedisplay();
   glutTimerFunc(10, update_p, 0);
}
/* Main function: GLUT runs as a console application starting at main() */
int main(int argc, char **argv)
{
   glutInit(&argc, argv);               // Initialize GLUT
   glutInitWindowSize(640, 500);        // Set the window's initial width & height - non-square
   glutInitWindowPosition(50, 50);      // Position the window's initial top-left corner
   glutCreateWindow("Clock"); // Create window with the given title
   glutDisplayFunc(display);            // Register callback handler for window re-paint event
   glutReshapeFunc(reshape);            // Register callback handler for window re-size event
   initGL();                            // Our own OpenGL initialization
   glutTimerFunc(0, update, 0);
   glutTimerFunc(0, update_p, 0);
   glutMainLoop(); // Enter the infinite event-processing loop
   return 0;
}