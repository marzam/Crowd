#include <Crowd.hpp>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
using namespace std;

//------------------------------------------------------------------------------
Crowd::Crowd(std::string fileName):
mDistance(1),
mCurrentState(NULL),
mBoard(NULL),
mPeople(NULL),
mTimeStep(0),
mCount(0){

  FILE *ptr = NULL;
  ptr = fopen(fileName.c_str(), "r");
  assert(ptr!=NULL);

  fscanf(ptr, "%d", &mCellX);
  fscanf(ptr, "%d", &mCellY);
  fscanf(ptr, "%d", &mNPeople);
  cout << "(" << mCellX << "," << mCellY << " | " << mNPeople << ")" << endl;

  assert(posix_memalign((void**)&mCurrentState, ALIGN, mCellX * mCellY * sizeof(int)) == 0);
  assert(posix_memalign((void**)&mBoard, ALIGN,  mCellX * mCellY * sizeof(int)) == 0);
  assert(posix_memalign((void**)&mPeople, ALIGN,  mNPeople * sizeof(stEntity)) == 0);

  //assert(posix_memalign((void**)&mCrowd, ALIGN,  mCellX * mCellY * sizeof(stCrowd)) == 0);

  for (int j = 0; j < mCellY; j++)
      for (int i = 0 ; i < mCellX; i++){
          int v = -1,
              p = j * mCellX + i;

          fscanf(ptr, "%d", &v);

          mBoard[p] = v;

      }


  fclose(ptr);
  mScaleX = 1.0f/static_cast<double>(mCellX);
  mScaleY = 1.0f/static_cast<double>(mCellY);
  clear();
  for (int i = 0; i < mNPeople; i++){
    mPeople[i].alive = true;
  }

}
Crowd::~Crowd(){
  if (mCurrentState != NULL){
    free(mCurrentState);
    mCurrentState = NULL;

  }



  if (mBoard != NULL){
    free(mBoard);
    mBoard = NULL;
  }

  if (mPeople != NULL){
    free(mPeople);
    mPeople = NULL;
  }


}

void Crowd::initialCondition(void){

  /*
    for (int j = 0 ; j < mCellY; j++){
        for (int i = 0; i < mCellX; i++){
            double prob = static_cast <double> (rand() % 65535 + 1) / 65535.0f;
            int p = (mCellX * j) + i;
            if (prob < 0.025f)
                mCurrentState[p] = 1;
            else
                mCurrentState[p] = 0;
        }
    }
*/

};
/*
    Moore neighborhood with r = 1

            |nw|n|ne|
            |w | |e |
            |sw|s|se|

 */
double Crowd::getStaticDynamicField(const int i, const int j){
    double x1 = static_cast<double>(i);
    double y1 = static_cast<double>(j);
    double x2 = static_cast<double>(mDoor.x);
    double y2 = static_cast<double>(mDoor.y);
    double x =  fabs(x1 - x2);
    double y =  fabs(y1 - y2);

    double l = x + y + 1;
    double o = 0.0f;
    double d = 0.0f;
    double e = 0.0f;
    int   r = mDistance;
    double E = 0.0f;

    if ((mCurrentState[j * mCellX + i] == empty) || (mCurrentState[j * mCellX + i] == door))
      E = 1.0f;

    for (int rj = -r; rj <= r; rj++){
      int pj = j - rj;
      for (int ri = -r; ri <= r; ri++){
          int pi = i + ri;
          int p = (pj * mCellX) + pi;
          if ((pi >= 0) && (pi < mCellX) && (pj >= 0) && (pj < mCellY)){
            if ((ri != 0) || (rj != 0)){
              if (mCurrentState[p] == empty)
                d++;
              if (mCurrentState[p] != wall)
                o++;
              if (mCurrentState[p] != door)
                e++;
            }//end-if ((ri != 0) && (rj != 0)){
          }//end-  if ((pi >= 0) && (pi < mCellX) && (pj >= 0) && (pj < mCellY)){
          //Bounderies condition represent nothing


      }//end-for (int ri = -r; ri <= r; ri++){
    }//ebdfor (int rj = -r; rj <= r; rj++){
//    cout << i << "\t" << j << "\t --> " << l << endl;
//double result = pow(E, exp((K1 * l) + (K2 * o) + (K3 * d)));/

//double result = (1.0 / exp(l * 0.125f)) * E; //(1.0f / l) * E; //pow(E, exp(K1 * l));

//double result = exp(d * 0.125f) * E; //(1.0f / l) * E; //pow(E, exp(K1 * l));
double result = (exp(o + d + (-e * 5)) / l) * E; //(1.0f / l) * E; //pow(E, exp(K1 * l));

return result;
}

void Crowd::applyRule(void){
  /*
    mTimeInstant+=time;
    if (time >= 0.0f){
    if (mTimeInstant < mDeltaT)
        return;

   }
*/
  double p[3][3];
  double n[3][3];

  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
      p[j][i] = n[j][i] = 0.0f;

  for (int index = 0; index < mNPeople; index++){
    if (mPeople[index].alive){
      int u = mPeople[index].x;
      int v = mPeople[index].y;
      //finding the center cell. In this case, it is 3x3
      int aux_i = 1;
      int aux_j = 1;
      int r = mDistance;
      double sum = 0.0f;
      for (int rj = -r; rj <= r; rj++){
        int pj = v - rj;
        for (int ri = -r; ri <= r; ri++){
            int pi = u + ri;
            double myN = 0.0f;
            if ((pi >= 0) && (pi < mCellX) && (pj >= 0) && (pj < mCellY)){
              myN  = getStaticDynamicField(pi, pj);
              sum += myN;
            }

            n[(2 - (aux_j + rj))][(aux_i + ri)] = myN;
        }//end-for (int ri = -r; ri <= r; ri++){
      }//ebdfor (int rj = -r; rj <= r; rj++){

      double max = -1.0f;
      for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++){
          p[j][i] = n[j][i] / sum;
          if (isnan(p[j][i]) || isinf(p[j][i]))
            p[j][i] = 0.0f;
          if (p[j][i] > max){
            int k = j * 3 + i;
            switch (k){
              case 0:mPeople[index].vx = -1; mPeople[index].vy = -1; break;
              case 1:mPeople[index].vx =  0; mPeople[index].vy = -1; break;
              case 2:mPeople[index].vx =  1; mPeople[index].vy = -1; break;
              case 3:mPeople[index].vx = -1; mPeople[index].vy =  0; break;
              //case 4:mPeople[index].vx = -1; mPeople[index].vy = -1; break;
              case 5:mPeople[index].vx =  1; mPeople[index].vy =  0; break;
              case 6:mPeople[index].vx = -1; mPeople[index].vy = +1; break;
              case 7:mPeople[index].vx =  0; mPeople[index].vy = +1; break;
              case 8:mPeople[index].vx =  1; mPeople[index].vy = +1; break;
              default:mPeople[index].vx =  0; mPeople[index].vy = 0 ; break;
            }
            max = p[j][i];
          }//end-for (int i = 0; i < 3; i++){
        }




    /*    cout << "N = " << endl;
        for (int j = 0; j < 3; j++){
          for (int i = 0; i < 3; i++){
            cout << n[j][i] << "\t";
          }
          cout << endl;
        }*/
        /*
        cout << "P = " << endl;
        cout << fixed << setprecision(15);
        for (int j = 0; j < 3; j++){
          for (int i = 0; i < 3; i++){
            cout << p[j][i] << "\t";
          }
          cout << endl;
        }
        cout << "\tPosition(" << u << "," << v << ")" << endl;
        cout << "\tVelocity(" <<  mPeople[index].vx << "," << mPeople[index].vy << ")" << endl;
        cout << "----------------------------------------------------" << endl;
*/
    }


  }//end-for (int i = 0; i < mNPeople; i++){



};

void Crowd::update(void){
  mTimeStep++;
  memcpy(mCurrentState, mBoard, sizeof(int)*mCellX*mCellY);

  for (int index = 0; index < mNPeople; index++){
    mPeople[index].x += mPeople[index].vx;
    mPeople[index].y += mPeople[index].vy;
    mPeople[index].vx = 0;
    mPeople[index].vy = 0;
    //cout << index << "(" << mPeople[index].x << "," << mPeople[index].y << ") " << mDoor.x << "," << mDoor.y << endl;
    if ((mPeople[index].x == mDoor.x) && (mPeople[index].y == mDoor.y) && (mPeople[index].alive)){
      mPeople[index].alive = false;
      mCount++;
    }

    if (mPeople[index].alive){
      int p = mPeople[index].y * mCellX + mPeople[index].x;
      assert(mCurrentState[p] != person);
      mCurrentState[p] = person;
    }

  }//end-for (int i = 0; i < mNPeople; i++){





};
void Crowd::clear(void){
    bzero(mCurrentState, sizeof(int)*mCellX*mCellY);
    bzero(mPeople, sizeof(stEntity)*mNPeople);

    memcpy(mCurrentState, mBoard, sizeof(int)*mCellX*mCellY);

    int k = 0;
    for (int j = 0; j < mCellY; j++){
        for (int i = 0; i < mCellX; i++){
            int p = j * mCellX + i;
            if (mBoard[p] == door){
              mDoor.x = i; mDoor.y = j;
            }//end- if (mBoard[p] == door){

            if (mBoard[p] == person){
              assert(k < mNPeople);
              mPeople[k].x = i;  mPeople[k].y = j;
              mBoard[p] = 0;
              k++;
            }
        }//end-for (int i = 0; i < mCellX; i++){
    }//end-for (int j = 0; j < mCellY; j++){

};



void Crowd::render(void){

    double dX = 0.0f,
          dY = 0.0f,
          r  = 0.0f,
          g  = 0.0f,
          b  = 0.0f;
    int type = GL_LINE_LOOP;


    //glScalef(, mScaleY, 1.0f);

    for (int j = 0; j < mCellY; j++){
        for (int i = 0; i < mCellX; i++){
            glPushMatrix();
            glTranslatef(dX, dY, 0.0f);
            int p = (mCellX * j) + i;
            if (mCurrentState[p] == empty){
                type = GL_LINE_LOOP;
                r = 0.0f;
                g = 0.0f;
                b = 0.0f;


            }else if (mCurrentState[p] == wall) {
              type = GL_QUADS;
              r = 1.0f;
              g = 1.0f;
              b = 1.0f;

            }else if (mCurrentState[p] == door) {
              type = GL_QUADS;
              r = 0.0f;
              g = 1.0f;
              b = 0.0f;
            }else if (mCurrentState[p] == person) {
              type = GL_QUADS;
              r = 1.0f;
              g = 1.0f;
              b = 0.0f;
            }
            glBegin(type); //GL_QUADS);GL_LINE_LOOP
            glColor3f(r, g, b);
            glVertex3f(0.0f, 0.0f, 0.0f);
            glVertex3f(0.0f, -mScaleY, 0.0f);
            glVertex3f(mScaleX, -mScaleY, 0.0f);
            glVertex3f(mScaleX, 0.0f, 0.0f);


            glEnd();
            glPopMatrix();
            dX += mScaleX;

        }
        dY -= mScaleY;
        dX = 0.0f;

    }



};

//Constructor
