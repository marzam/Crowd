#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>

#define PI      4.0f * atan(1.0f)
#define SQUARE  sqrt(2.0f * PI)


using namespace std;
double normal (const double x,
               const double ave,
               const double var){

   double a = (x - ave) / var;
   a *= a;
   return exp(-0.5f * a) / (var * SQUARE);

}
double normal (const double x0,
               const double x1,
               const double h,
               const double ave,
               const double var){
    double x  = x0;
    double h2 = h / 2.0f;
    double n  = normal(x, ave, var);
    double I1 = n * h2;


    do{
      cout << setprecision(5) << fixed <<  x << " " << n << " " << I1 << endl;
      x  += h;
      n   =  normal(x, ave, var);
      I1 += (2.0f * n * h2);
    }while (x < (x1 - h));
    n   =  normal(x, ave, var);
    I1 += (2.0f * n * h2);
    cout << x << "\t" << n << "\t" << I1 << endl;
    return I1;
}


int main (int ac, char **av){
    cout << "# normal" << endl;
    cerr << normal(0.0f, 1.0f, 0.00001f, 0.25f, 0.01f) << endl;
    return EXIT_SUCCESS;
}
