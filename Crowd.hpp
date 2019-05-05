#ifndef _CROWD_HPP_
#define _CROWD_HPP_
#include <string>

//#define DELTA_T 100.0f


struct stEntity{
  bool alive;
  int x, y, vx, vy;
};



class Crowd {
private:
  const double K1 = 1.0f;
  const double K2 = 1.0f;
  const double K3 = 1.0f;

  enum STATES{empty, wall, door, person};
  double getStaticDynamicField(const int, const int);
public:
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
    void render(void);

    void setDistance(int d) { this->mDistance = d; };

    int getCount(void)      { return this->mCount; }
    int getTimeStep(void)   { return this->mTimeStep; }
  protected:
    int mCellX,
        mCellY,
       *mCurrentState,
       *mBoard,
        mDistance,
        mNPeople,
        mTimeStep,
        mCount;


    float mScaleX,
          mScaleY,
          mCells;


    stEntity mDoor;
    stEntity *mPeople;

};

#endif
