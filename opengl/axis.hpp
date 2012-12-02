# ifndef AXIS_H
# define AXIS_H

# include "graphicsDefs.hpp"



point3 primaryColors[] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

point3 worldBasis[3] = {{ 1.0, 0.0, 0.0 },
                                                 { 0.0, 1.0, 0.0 },
                                                 { 0.0, 0.0, 1.0 }};
point3 worldOrigin = {0.0, 0.0, 0.0};


// draws an axis system at any frame given by 

void drawAxis(point3 basis[3], point3 origin) {
         int j, k;
         point3 endPoints[3];
         for (j = 0; j < 3; ++j)
                 for ( k = 0; k < 3; ++k)
                         endPoints[j][k] = basis[j][k] + origin[k];

         glLineWidth(3);
         glBegin(GL_LINES);
            for (j = 0; j < 3; ++j)
            {
                         glColor3fv(primaryColors[j]);
                         glVertex3fv(origin);
                         glVertex3fv(endPoints[j]);
            }
         glEnd();
}


void drawGrid(float left, float right, float back, float front, float interval) {
         float j;
         glColor3f(0.25, 0.25, 0.25);
         glLineWidth(0.5);
         glBegin(GL_LINES);
                 for (j = left; j <= right; j += interval)
                 {
                         glVertex3f(j, 0, back);
                         glVertex3f(j, 0, front);
                 }
                 for (j = back; j <= front; j += interval)
                 {
                         glVertex3f(left, 0, j);
                         glVertex3f(right, 0, j);
                 }
         glEnd();
}


#endif
