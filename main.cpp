#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <Crowd.hpp>
#include <App.hpp>

using namespace std;

Stopwatch stopwatch;
double    elapsedTime = 0.0f;

float             FPS = 0.0f;

int               width = 800,
                  height = 600,
                  mCurrentTimeStep = 0;

Crowd *cellularAutomata = NULL;
// Função de visualização
void cellularAutomataRender(void){

      double dX = 0.0f,
            dY = 0.0f,
            r  = 0.0f,
            g  = 0.0f,
            b  = 0.0f;
      int type = GL_LINE_LOOP;


      //glScalef(, mScaleY, 1.0f);

      for (int j = 0; j < cellularAutomata->getCellY(); j++){
          for (int i = 0; i < cellularAutomata->getCellX(); i++){
              glPushMatrix();
              glTranslatef(dX, dY, 0.0f);
              int value = cellularAutomata->getCellValue(i, j);
              if (value == Crowd::empty){
                  type = GL_LINE_LOOP;
                  r = 0.0f;
                  g = 0.0f;
                  b = 0.0f;


              }else if (value == Crowd::wall) {
                type = GL_QUADS;
                r = 1.0f;
                g = 1.0f;
                b = 1.0f;

              }else if (value == Crowd::door) {
                type = GL_QUADS;
                r = 0.0f;
                g = 1.0f;
                b = 0.0f;
              }else{ //if (mCurrentState[p] == person) {
                type = GL_QUADS;
                r = 1.0f;
                g = 1.0f;
                b = 0.0f;
              }
              glBegin(type); //GL_QUADS);GL_LINE_LOOP
              glColor3f(r, g, b);
              glVertex3f(0.0f, 0.0f, 0.0f);
              glVertex3f(0.0f, -cellularAutomata->getScaleY(), 0.0f);
              glVertex3f(cellularAutomata->getScaleX(), -cellularAutomata->getScaleY(), 0.0f);
              glVertex3f(cellularAutomata->getScaleX(), 0.0f, 0.0f);


              glEnd();
              glPopMatrix();
              dX += cellularAutomata->getScaleX();

          }
          dY -= cellularAutomata->getScaleY();
          dX = 0.0f;

      }


}
void render(void){
    glClear(GL_COLOR_BUFFER_BIT);


    cellularAutomataRender();

    float delta  = 0.0f;
    glColor3f(0.0f, 0.0f, 1.0f);
//    for (int i = 0; i < cellularAutomata->getCellX(); i++)

    for (int j = 0; j <= cellularAutomata->getCellY(); j++){
        glBegin(GL_LINES); //GL_QUADS);GL_LINE_LOOP
        glVertex3f(0.0f, -delta, 0.0f);
        glVertex3f(1.0f, -delta, 0.0f);
        glEnd();
        delta += cellularAutomata->getScaleY();
    }

    delta = 0.0f;
    for (int j = 0; j <= cellularAutomata->getCellX(); j++){
        glBegin(GL_LINES); //GL_QUADS);GL_LINE_LOOP
        glVertex3f(delta, 0.0f, 0.0f);
        glVertex3f(delta, -1.0f, 0.0f);
        glEnd();
        delta += cellularAutomata->getScaleX();
    }


    glutSwapBuffers();


}

// Função de inicialização de parâmetros (OpenGL e outros)
void init (void){
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
}

// Função de evento do teclado
void keyboardEvent(unsigned char key, int x, int y)
{
     glutPostRedisplay();
    switch (key) {
        case 'b':
        case 'B':

            break;

        case 'e':
        case 'E':

            break;
        case 'c':
        case 'C':
            cellularAutomata->clear();
            break;

        case 'r':
        case 'R':
            cellularAutomata->initialCondition();
            break;

        case 's':
        case 'S':
            cellularAutomata->applyRule();
            cellularAutomata->update();
            break;


        case 'q':
        case 'Q':
        case 27:
            exit (EXIT_SUCCESS);
            break;




        default:
            break;
   }
   glutPostRedisplay();
}

// Função de evento do mouse
void mouseEvent(int button, int state, int x, int y){

    if (button == GLUT_LEFT_BUTTON)
        if (state == GLUT_DOWN){
            float x1 = (static_cast<float>(x) / width) * static_cast<float>(cellularAutomata->getCellX());
            float y1 = (static_cast<float>(y) / height) * static_cast<float>(cellularAutomata->getCellY());
            int x2 = static_cast<int>(x1)  ;
            int y2 = static_cast<int>(y1) ;
            //x2 /=  cellularAutomata->getCellX();
            //y2 /=  cellularAutomata->getCellY();

            //cellularAutomata->changeState(x2, y2);

        }

    if (button == GLUT_RIGHT_BUTTON)
        if (state == GLUT_DOWN)
            cerr << "Right button" << endl;

   glutPostRedisplay();
}

//Viewport
void viewPort(int w, int h){

    if(h == 0) h = 1;


    width = w;
    height = h;
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho (-0.01f, 1.01f, -1.01f, 0.01f, -1.0f, 1.0f);
    glutPostRedisplay();
}

//Loop principal da visualização
void mainloop(void){
    glutPostRedisplay();

   char msg[1024];
   sprintf(msg, "Crowd\t %d/%d", cellularAutomata->getCount(), cellularAutomata->getTimeStep());
   glutSetWindowTitle(msg);


}


int main(int argc, char**argv){


    cellularAutomata = new Crowd("board3.txt");
    START_STOPWATCH(stopwatch);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(width, height);
    glutCreateWindow("Cellular Automata - Game of life");
    glutDisplayFunc(render);
    glutReshapeFunc(viewPort);
    glutMouseFunc(mouseEvent);
    glutKeyboardFunc(keyboardEvent);
    glutIdleFunc(mainloop);
    init();
    glutMainLoop();
    return EXIT_SUCCESS;
}
