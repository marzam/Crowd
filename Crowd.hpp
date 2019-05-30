#ifndef _CROWD_HPP_
#define _CROWD_HPP_
#include <string>
#include <cstdint>
//#define DELTA_T 100.0f
double gamaFunction (double n);

struct stEntity{
  bool alive;
  int id,
      stopET,
      x0,
      y0,
      x,
      y,
      vx,
      vy;
  float prob;
};



class Crowd {
private:
  const double Ko = 1.0f;//1.0f; //-0.25f; //-1.0f;   //Obstacle repulsion
  const double Kd = 1.0f; //-1.0f; //  //Is empty
  const double Ke = 1.0f; //1.0f;   //Is door exit
  const double Kl = -5.0f ;//-1.0f;  //Distance
  const double ER = 1E-20;
  double getStaticDynamicField(const int, const int);
public:
  enum STATES{empty, wall, door, person};
  enum PERSON{walking, stopped, passed, dieded}; //states of person. In door cell

   Crowd(std::string );
    ~Crowd();
    int getCellX(void)    { return mCellX; };
    int getCellY(void)    { return mCellY; };
    float getScaleX(void) { return mScaleX; };
    float getScaleY(void) { return mScaleY; };

    void initialCondition(void);
    void applyRule(void);
    void update(void);
    void clear(void);
    //void render(void);

    void setDistance(int d) { this->mDistance = d; };

    int getCount(void)      { return this->mCount; }
    int getTimeStep(void)   { return this->mTimeStep; }
    int getCellValue(int i, int j)
                            { return static_cast<int>(this->mMesh[j * this->mCellX + i]);}
    float myRand(void);
    void printProb(int);
    void setBeta(const double alpha, const double beta);
    void debug(void);

//    void insertPerson(void);
//    void removePerson(void);
  protected:
    void setPosition(stEntity *ptrPerson, int level);
    void saveLog(void);
    double fillDynamicMatrix(double *out_n,
                                     int u,
                                     int v,
                                     int l);
    void fillProbMatrix(double *out_p,
                        int u,
                        int v);
    double betaFunction (double, double, double, double, double, double );
    double betaFunction (double *);

    int mCellX,
        mCellY,
       *mBoard,
        mDistance,
        mNPeople,
        mTimeStep,
        mCount;

    uintptr_t *mMesh;

    double mVarTime,
           mAveTime;

    float mScaleX,
          mScaleY,
          mCells;

    double     mParam[5];


    stEntity mDoor;
    stEntity *mPeople;

};

#endif
