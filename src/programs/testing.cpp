#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>

#include "gnuplot-iostream.h"
#include "casadi/casadi.hpp"

#include <unistd.h>
inline void mysleep(unsigned millis) {
  ::usleep(millis * 1000);
}

void demo_animation() {

    Gnuplot gp;

    std::cout << "Press Ctrl-C to quit (closing gnuplot window doesn't quit)." << std::endl;

    gp << "set yrange [-1:1]\n";

    const int N = 1000;
    std::vector<double> pts(N);

    double theta = 0;
    while(1) {
        for(int i=0; i<N; i++) {
            double alpha = (static_cast<double>(i)/N-0.5) * 10;
            pts[i] = sin(alpha*8.0 + theta) * exp(-alpha*alpha/2.0);
        }

        gp << "plot '-' binary" << gp.binFmt1d(pts, "array") << "with lines notitle\n";
        gp.sendBinary1d(pts);
        gp.flush();

        theta += 0.2;
        mysleep(100);
    }
}



void demo_basic() {
  Gnuplot gp;
  // For debugging or manual editing of commands:
  //Gnuplot gp(std::fopen("plot.gnu", "w"));
  // or
  //Gnuplot gp("tee plot.gnu | gnuplot -persist");

  std::vector<std::pair<double, double>> xy_pts_A;
  for(double x=-2; x<2; x+=0.01) {
    double y = x*x*x;
    xy_pts_A.emplace_back(x, y);
  }

  std::vector<std::pair<double, double>> xy_pts_B;
  for(double alpha=0; alpha<1; alpha+=1.0/24.0) {
    double theta = alpha*2.0*3.14159;
    xy_pts_B.emplace_back(cos(theta), sin(theta));
  }

  gp << "set xrange [-2:2]\nset yrange [-2:2]\n";
  gp << "plot '-' with lines title 'cubic', '-' with points title 'circle'\n";
  gp.send1d(xy_pts_A);
  gp.send1d(xy_pts_B);

}


void demo_tmpfile() {
    Gnuplot gp;

    std::vector<std::pair<double, double>> xy_pts_A;
    for(double x=-2; x<2; x+=0.01) {
        double y = x*x*x;
        xy_pts_A.emplace_back(x, y);
    }

    std::vector<std::pair<double, double>> xy_pts_B;
    for(double alpha=0; alpha<1; alpha+=1.0/24.0) {
        double theta = alpha*2.0*3.14159;
        xy_pts_B.emplace_back(cos(theta), sin(theta));
    }

    gp << "set xrange [-2:2]\nset yrange [-2:2]\n";
    // Data will be sent via a temporary file.  These are erased when you call
    // gp.clearTmpfiles() or when gp goes out of scope.  If you pass a filename
    // (i.e. `gp.file1d(pts, "mydata.dat")`), then the named file will be created
    // and won't be deleted.
    //
    // Note: you need std::endl here in order to flush the buffer.  The send1d()
    // function flushes automatically, but we're not using that here.
    gp << "plot" << gp.file1d(xy_pts_A) << "with lines title 'cubic',"
        << gp.file1d(xy_pts_B) << "with points title 'circle'" << std::endl;

}



int main()
{
  std::cout << "Hello" << std::endl;

  demo_basic();
 // demo_tmpfile();
  demo_animation();

  return 0;
}
