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
#include <climits>
#include <fstream>
using namespace std;
//------------------------------------------------------------------------------
double gamaFunction (double n){

    unsigned long in1 = static_cast <unsigned long> (n-1),
    iacc = in1;

    if (n <= 2.0f) return 1.0f;

    for (unsigned long i = in1-1 ; i > 0 ; i--)
        iacc *= i;

    return static_cast <double> (iacc);
}
//------------------------------------------------------------------------------
Crowd::Crowd(std::string fileName):
mDistance(1),
mBoard(NULL),
mPeople(NULL),
mTimeStep(0),
mCount(0),
mMesh(NULL),
mVarTime(0.0f),
mAveTime(0.0f){

  srand (time(NULL));
  FILE *ptr = NULL;
  ptr = fopen(fileName.c_str(), "r");
  assert(ptr!=NULL);

  fscanf(ptr, "%d", &mCellX);
  fscanf(ptr, "%d", &mCellY);
  fscanf(ptr, "%d", &mNPeople);
  cout << "(" << mCellX << "," << mCellY << " | " << mNPeople << ")" << endl;

  assert(posix_memalign((void**)&mMesh, ALIGN, mCellX * mCellY * sizeof(uintptr_t)) == 0);
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


  initialCondition();
  cout << "Ko = " << Ko << endl;
  cout << "Kd = " << Kd << endl;
  cout << "Ke = " << Ke << endl;
  cout << "Kl = " << Kl << endl;
  cout << "-----------------------------" << endl << endl;
}
Crowd::~Crowd(){
  if (mMesh != NULL){
    free(mMesh);
    mMesh = NULL;

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
      int id = 0;
      for (int i = 0; i < mNPeople; i++){
        mPeople[i].id = ++id;
        mPeople[i].alive = true;
        mPeople[i].stopET = walking;
        fprintf (stderr, "This is line %d.\n", __LINE__);

        mPeople[i].prob = myRand();
      }
      int k = 0;
      for (int j = 0; j < mCellY; j++){
          for (int i = 0; i < mCellX; i++){
              int p = j * mCellX + i;
              if (mBoard[p] == door){
                mDoor.x = i; mDoor.y = j;
                mBoard[p] = 0;
              }//end- if (mBoard[p] == door){

              if (mBoard[p] == person){
                assert(k < mNPeople);
                mPeople[k].x = i;  mPeople[k].y = j;
                mPeople[k].x0 = i;  mPeople[k].y0 = j;
                mBoard[p] = 0;
                k++;
              }

              mMesh[p] = mBoard[p];
          }//end-for (int i = 0; i < mCellX; i++){
      }//end-for (int j = 0; j < mCellY; j++){
      update();
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

    double l = x + y + 1E-5;
    double o = 0.0f;
    double d = 0.0f;
    double e = 0.0f;
    int   r = mDistance;
    double E = 0.0f;

    assert((i >= 0) && (i < mCellX) && (j >= 0) && (j < mCellY));
    if ((mMesh[j * mCellX + i] == empty) || (mMesh[j * mCellX + i] == door))
      E = 1.0f;

    for (int rj = -r; rj <= r; rj++){
      int pj = j - rj;
      for (int ri = -r; ri <= r; ri++){
          int pi = i + ri;
          int p = (pj * mCellX) + pi;
          if ((pi >= 0) && (pi < mCellX) && (pj >= 0) && (pj < mCellY)){
            if ((ri != 0) || (rj != 0)){
              bool isDoor = ((pi == mDoor.x) && (pj == mDoor.y));
              //if (isDoor)
              //cerr << "Door" << endl;
              if ((mMesh[p] == empty) && !isDoor)
                d++;
              else if ((mMesh[p] != empty) && !isDoor) //((mMesh[p] != wall)&& isDoor)
                o++;
              if (isDoor)// && !isDoor)
                e++;
            }//end-if ((ri != 0) && (rj != 0)){
          }//end-  if ((pi >= 0) && (pi < mCellX) && (pj >= 0) && (pj < mCellY)){
          //Bounderies condition represent nothing


      }//end-for (int ri = -r; ri <= r; ri++){
    }//ebdfor (int rj = -r; rj <= r; rj++){

  double alpha_o = 1.0f; //betaFunction(mParam);
  double alpha_d = 1.0f; //betaFunction(mParam);
  double alpha_e = 1.0f; //betaFunction(mParam);
  double alpha_l = 1.0f; //betaFunction(mParam);

  double in_o = (o * Ko) * alpha_o;
  double in_d = (d * Kd) * alpha_d;
  double in_e = (e * Ke) * alpha_e;//
  double in_l = (l * Kl) * alpha_l;
  //double ee   = exp(in_o + in_d + in_e + in_l);
  double ee   = exp(in_l + in_o + in_d);
  //cout << "\te = " << ee << endl;
  //cout << "\t_l = " << l << endl;
  //exp(in_o + in_d + in_e)
  //double ee;//result = (ee / l) * E; //(1.0f / l) * E; //pow(E, exp(K1 * l));
  double result = ee * E;
  return result;
}
double Crowd::fillDynamicMatrix(double *out_n,
                                        int u,
                                        int v,
                                        int l){

      if (l == 0) return 0;
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
            int k = ((2 - (aux_j + rj))) * 3 + ((aux_i + ri));
            out_n[k] = myN;
        }//end-for (int ri = -r; ri <= r; ri++){
      }//ebdfor (int rj = -r; rj <= r; rj++){

      return sum;
}

void Crowd::fillProbMatrix(double *out_p,
                              int u,
                              int v){
  //double p[9];
  double n[9];
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++){
      int k = j * 3 + i;
      out_p[k] = n[k] = 0.0f;
    }
    double sum = fillDynamicMatrix(&n[0], u, v, 1);
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++){
        int k = j * 3 + i;
        out_p[k] = n[k] / sum;
        if (isnan(out_p[k]) || isinf(out_p[k]))
          out_p[k] = 0.0f;
      }//end-for (int j = 0; j < 3; j++)

}

void Crowd::applyRule(void){
  /*
    mTimeInstant+=time;
    if (time >= 0.0f){
    if (mTimeInstant < mDeltaT)
        return;

   }
*/

  double p[9];
  for (int index = 0; index < mNPeople; index++){
    if ((mPeople[index].alive) &&  (mPeople[index].stopET == walking)){

      int u = mPeople[index].x;
      int v = mPeople[index].y;
     

      fillProbMatrix(&p[0], u, v);
      //finding the center cell. In this case, it is 3x3
      //#AQUI
      double prob[9] = {0.0f};
      int    posi[9] = {0};

      double max = -1.0f;
      int    idx = 0;

      for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++){
            int k = j * 3 + i;
          if (p[k] > max){
            max = p[k];
            double distance = abs(max - prob[idx]);
            if (distance < ER){
              idx++;
              assert(idx < 9);
            }//end-if (distance < ER){
            posi[idx] = k;
            prob[idx] = max;
          }//end-if (p[k] > max){
        }//end-for (int i = 0; i < 3; i++){

        for (int i = 0; i < idx - 1; i++){
          double distance = abs(prob[idx + 1] - prob[idx]);
          assert(distance < ER);
        }
        if (idx > 1){
          cerr << ".";
        }
        int l = posi[rand() % (idx + 1)];
        switch (l){
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
/*
      double max = -1.0f;
      for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++){
            int k = j * 3 + i;
          if (p[k] > max){
            int l = j * 3 + i;
            switch (l){
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
            max = p[k];
          }//end-for (int i = 0; i < 3; i++){
        }//end-for (int j = 0; j < 3; j++)
*/
    }//end-if ((mPeople[index].alive) &&  (mPeople[index].stopET == walking)){


  }//end-for (int i = 0; i < mNPeople; i++){



};

void Crowd::setPosition(stEntity *ptrPerson, int level){
  int p1 = ptrPerson->y * mCellX + ptrPerson->x;
  if (mMesh[p1] == empty){
      mMesh[p1] = reinterpret_cast<uintptr_t> (ptrPerson);

      if ((ptrPerson->x == mDoor.x) && (ptrPerson->y == mDoor.y)){
        if (ptrPerson->stopET == walking ){
          ptrPerson->stopET = stopped;

        }else if (ptrPerson->stopET == stopped){
          mMesh[p1] = empty;
          ptrPerson->stopET = passed;
          ptrPerson->alive = false;
          mCount++;
          mAveTime += static_cast<double>(mTimeStep);
          mVarTime += static_cast<double>(mTimeStep * mTimeStep);

        }
        //end-if (ptrPerson->stopET == walking ){
        //ptrPerson->alive = false;

      }//end-if ((ptrPerson->x == mDoor.x) && (ptrPerson->y == mDoor.y) && (ptrPerson->alive)){
      ptrPerson = NULL;
      return;
    }//end-if (mMesh[p1] == empty){

    if ((mMesh[p1] == wall)){ // || (p1 == door)){
      ptrPerson->x = ptrPerson->x0;
      ptrPerson->y = ptrPerson->y0;
      int p2 = ptrPerson->y * mCellX + ptrPerson->x;
      mMesh[p2] = reinterpret_cast<uintptr_t> (ptrPerson);
      ptrPerson = NULL;
      return;
    }/*else{
      cerr << "All your bases are belong to us!!!!"  << endl;
      exit(-1);
    }*/

    assert((mMesh[p1] != wall) && (mMesh[p1] != empty));
    stEntity *ptrPerson2 = reinterpret_cast<stEntity*>(mMesh[p1]);
    
    if (!ptrPerson2->alive){
          mMesh[p1] = reinterpret_cast<uintptr_t> (ptrPerson);
          return;
    }
    
    if (ptrPerson2->id == ptrPerson->id)
            return;
    
    if (level > 1000)
        cerr << "." << endl;
    if (ptrPerson->prob >= ptrPerson2->prob){
      p1 = ptrPerson->y * mCellX + ptrPerson->x;
      mMesh[p1] = reinterpret_cast<uintptr_t> (ptrPerson);

      ptrPerson2->x = ptrPerson2->x0;
      ptrPerson2->y = ptrPerson2->y0;
      ptrPerson2->prob *= myRand();
      ptrPerson->prob *= myRand();

      cerr << "402 L = " << level << " " << ptrPerson2->id << " C(" << ptrPerson2->x << "," << ptrPerson2->y << ")" << " v(" << ptrPerson2->vx << "," << ptrPerson2->vy << ")" << " C0(" << ptrPerson2->x0 << "," << ptrPerson2->y0 << ") P(" << ptrPerson2->prob << ")" << endl;
      cerr << "        " << level << " " << ptrPerson->id << " C(" << ptrPerson->x << "," << ptrPerson->y << ")" << " v(" << ptrPerson->vx << "," << ptrPerson->vy << ")"<< " C0(" << ptrPerson->x0 << "," << ptrPerson->y0 << ") P(" << ptrPerson->prob << ")" << endl;

      ptrPerson = NULL;

      setPosition(ptrPerson2, level+1);
    }else{
      p1 = ptrPerson2->y * mCellX + ptrPerson2->x;
      if (ptrPerson2->id == ptrPerson->id)
        return;

      mMesh[p1] = reinterpret_cast<uintptr_t> (ptrPerson2);

      ptrPerson->x = ptrPerson->x0;
      ptrPerson->y = ptrPerson->y0;
      ptrPerson->prob *= myRand();
      ptrPerson2->prob *= myRand();

      cerr << "2 > 1 L = " << level << " " << ptrPerson->id << " C(" << ptrPerson->x << "," << ptrPerson->y << ")" << " v(" << ptrPerson->vx << "," << ptrPerson->vy << ")" << " C0(" << ptrPerson->x0 << "," << ptrPerson->y0 << ") P(" << ptrPerson->prob << ")" << endl;
      cerr << "      *" << level << " " << ptrPerson2->id << " C(" << ptrPerson2->x << "," << ptrPerson2->y << ")" << " v(" << ptrPerson2->vx << "," << ptrPerson2->vy << ")" << " C0(" << ptrPerson2->x0 << "," << ptrPerson2->y0 << ") P(" << ptrPerson2->prob << ")" << endl;

      ptrPerson2 = NULL;
      setPosition(ptrPerson, level+1);
    }




}
void Crowd::update(void){
    mTimeStep++;
    for (int i = 0; i < mCellX*mCellY; i++){
      mMesh[i] = static_cast<uintptr_t>(mBoard[i]);
    }

    
    
    
    for (int index = 0; index < mNPeople; index++){
        stEntity *ptrPerson = &mPeople[index];
        if (ptrPerson->stopET == passed)
            ptrPerson->stopET = dieded;
        
        ptrPerson->x0  = ptrPerson->x;
        ptrPerson->y0  = ptrPerson->y;
        ptrPerson->x += ptrPerson->vx;
        ptrPerson->y += ptrPerson->vy;
      
       
      //end-if (ptrPerson->alive){
    }//end-for (int i = 0; i < mNPeople; i++){

    for (int index = 0; index < mNPeople; index++){
        stEntity *ptrPerson = &mPeople[index];
        if (ptrPerson->alive)
            setPosition(ptrPerson, 1);
      //end-if (ptrPerson->alive){
    }//end-for (int i = 0; i < mNPeople; i++){
    
    
    for (int index = 0; index < mNPeople; index++){
      stEntity *ptrPerson = &mPeople[index];
    
      if (ptrPerson->alive){
        ptrPerson->vx = 0;
        ptrPerson->vy = 0;
        fprintf (stderr, "This is line %d.\n", __LINE__);
        ptrPerson->prob = myRand();
      }//end-if (ptrPerson->alive){
   }//end-for (int i = 0; i < mNPeople; i++){

   saveLog();

};
void Crowd::clear(void){
    bzero(mMesh, sizeof(uintptr_t)*mCellX*mCellY);
    bzero(mPeople, sizeof(stEntity)*mNPeople);
};





float Crowd::myRand(void){
  return static_cast <float> (rand() % INT_MAX + 1) / static_cast <float> (INT_MAX);

};

void Crowd::printProb(int index){
  double p[3][3];
  double n[3][3];

  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
      p[j][i] = n[j][i] = 0.0f;


    if (mPeople[index].alive){
      int u = mPeople[index].x;
      int v = mPeople[index].y;
      fprintf (stderr, "This is line %d.\n", __LINE__);

      mPeople[index].prob = myRand();

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
      cout << "P = " << endl;
      for (int j = 0; j < 3; j++){
        for (int i = 0; i < 3; i++){
          p[j][i] = n[j][i] / sum;
          if (isnan(p[j][i]) || isinf(p[j][i]))
            p[j][i] = 0.0f;
          cout << setprecision(4) << setw(12) << fixed << p[j][i] ;
        }//end-for (int i = 0; i < 3; i++){
        cout << endl;
      }
    }
}
/*
void Crowd::statistic(void){

}
*/
void Crowd::saveLog(void){
  fstream output;
  string fileName = "log-8-4.txt";
  output.open(fileName, std::fstream::out | std::fstream::app);
   assert(output.is_open());
   for (int index = 0; index < mNPeople; index++){
     if ((!mPeople[index].alive) && (mPeople[index].stopET == passed)){
       output << setw(10) << mTimeStep  << " " << mCount << " " << mPeople[index].id << endl;
     }//end-if ((mPeople[index].alive) && (mPeople[index],stopET == passed)){
   }//end-for (int index = 0; index < mNPeople; index++){


   output.close();
};

void Crowd::setBeta(const double alpha, const double beta){
  mParam[0] = alpha;
  mParam[1] = beta;
  mParam[2] = gamaFunction(mParam[0]);
  mParam[3] = gamaFunction(mParam[1]);
  mParam[4] = gamaFunction(mParam[0]+mParam[1]);


/*
  mRules->getGrid()->getVehicleType(*nType)->param[2] = gamaFunction(mRules->getGrid()->getVehicleType(*nType)->param[0]);
  mRules->getGrid()->getVehicleType(*nType)->param[3] = gamaFunction(mRules->getGrid()->getVehicleType(*nType)->param[1]);
  mRules->getGrid()->getVehicleType(*nType)->param[4] = gamaFunction(mRules->getGrid()->getVehicleType(*nType)->param[0]+mRules->getGrid()->getVehicleType(*nType)->param[1]);
*/
};
double  Crowd::betaFunction (double gamaAlpha,
                                    double gamaBeta,
                                    double gamaAlphaBeta,
                                    double alpha,
                                    double beta,
                                    double x){


	double  xalpha = pow(x, (alpha-1.0f)),
                xbeta  = pow((1.0f - x), (beta-1.0f));

	double ret = gamaAlphaBeta / (gamaAlpha * gamaBeta);

    ret *= xalpha * xbeta;

    return ret;
};

double  Crowd::betaFunction (double *vet){

	double x, y, y1;
    do{
        x  = static_cast <double> (rand() % 65535 + 1) / 65535.0f;
        y1 = static_cast <double> (rand() % 65535 + 1) / 65535.0f;

        y = betaFunction( vet[2],
                         vet[3],
                         vet[4],
                         vet[0],
                         vet[1],
                         x);


    }while (y1 > y);

    if (x < 0.0f)
        x = 0.0f;

    if (x > 1.0f)
        x = 1.0f;

    return x;

};

void Crowd::debug(void){
  cout << "---------------------------------- ---------------------" << endl;
  for (int index = 0; index < mNPeople; index++){
    stEntity person = mPeople[index];
      cout << person.id << " (" << person.x << "," << person.y << ") Alive:" << person.alive << " P(" << person.prob << ")" << endl;
  }
  cout.flush();
};
