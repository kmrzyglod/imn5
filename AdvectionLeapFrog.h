//
// Created by Kamil on 22.05.2016.
//

#ifndef IMN5_ADVECTIONLEAPFROG_H
#define IMN5_ADVECTIONLEAPFROG_H

#include "Advection.h"

class AdvectionLeapFrog: public Advection {
protected:
    void fillAdvectionMatrix() {
        for(int i=0;i< _xsize;i++) {
            for(int j=0;j < _ysize;j++) {
                double a = pow((i * _dx - 0.4), 2);
                double b = pow((j * _dy - 0.45), 2);
                _density0Matrix[i][j] = exp(-25.0*(a+b));
                a = pow((i * _dx - _uVelocity[i][j]*_deltaT - 0.4), 2);
                b = pow((j * _dy - _vVelocity[i][j]*_deltaT  - 0.45), 2);
                _density1Matrix[i][j] =  exp(-25.0*(a+b));

            }
        }
    }
    void makeAdvection() {
        for(int i=0;i< _xsize;i++) {
            for(int j=1;j < _ysize-1;j++) {
                if(i == 0) {
                    double a = _uVelocity[i][j] * ((_density0Matrix[i+1][j] - _density0Matrix[_xsize-1][j])/_dx);
                    double b = _vVelocity[i][j] * ((_density0Matrix[i][j+1] - _density0Matrix[i][j-1])/_dy);
                    _density2Matrix[i][j] = _density1Matrix[i][j] - _deltaT*(a+b);
                }
                else if(i == _xsize - 1) {
                    double a = _uVelocity[i][j] * ((_density0Matrix[0][j] - _density0Matrix[i-1][j])/_dx);
                    double b = _vVelocity[i][j] * ((_density0Matrix[i][j+1] - _density0Matrix[i][j-1])/_dy);
                    _density2Matrix[i][j] = _density1Matrix[i][j] - _deltaT*(a+b);
                }
                else {
                    double a = _uVelocity[i][j] * ((_density0Matrix[i + 1][j] - _density0Matrix[i - 1][j]) / _dx);
                    double b = _vVelocity[i][j] * ((_density0Matrix[i][j + 1] - _density0Matrix[i][j - 1]) / _dy);
                    _density2Matrix[i][j] = _density1Matrix[i][j] - _deltaT * (a + b);
                }
            }
        }
    }
    void SaveResults(char* filename) {
        imnd::write_data2D(filename, _xsize, _ysize, _dx, _dy);
    }
public:
    AdvectionLeapFrog(int xsize, int ysize, double dx, double dy, FlagMatrix* flagMatrix,  function < void(double** uMatrix, double** vMatrix) > fillVeloctyMatrixFnc)
            : Advection(xsize, ysize, dx, dy, flagMatrix, fillVeloctyMatrixFnc) {}
};

#endif //IMN5_ADVECTIONLEAPFROG_H
