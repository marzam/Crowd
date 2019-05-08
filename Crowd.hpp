#ifndef _CROWD_HPP_
#define _CROWD_HPP_
#include <string>
#include <cstdint>
//#define DELTA_T 100.0f


struct stEntity{

  bool alive;

  int id, stopET, x0, y0, x, y, vx, vy;
  float prob;
};



class Crowd {
private:
  const double K1 = 1.0f;
  const double K2 = 1.0f;
  const double K3 = 1.0f;

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
//    void insertPerson(void);
//    void removePerson(void);
  protected:
    void saveLog(void);
    int mCellX,
        mCellY,
       *mBoard,
        mDistance,
        mNPeople,
        mTimeStep,
        mCount;

    uintptr_t *mMesh;

    float mScaleX,
          mScaleY,
          mCells;


    stEntity mDoor;
    stEntity *mPeople;

};

#endif
