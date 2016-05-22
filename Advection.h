//
// Created by Kamil on 22.05.2016.
//

#ifndef IMN5_ADVECTION_H
#define IMN5_ADVECTION_H

#include <functional>
#include "FlagMatrix.h"
class Advection {
protected:
    double** _density0Matrix;
    double** _density1Matrix;
    double** _density2Matrix;

    double** _uVelocity;
    double** _vVelocity;

    double _maxU, _maxV, _deltaT, _nowT;

    int _xsize, _ysize, _iter;
    double  _dx, _dy;
    FlagMatrix* _flagMatrix;
    //function < double(int x, int y) > _denisty0Fnc;
    function < void(double** uMatrix, double** vMatrix) > _fillVelocityMatrixFnc;

    Advection() {}

    virtual void makeAdvection()  =  0;
    virtual void fillAdvectionMatrix() = 0;

    void fillVelocityMatrixAndDeltaT() {
        _fillVelocityMatrixFnc(_uVelocity, _vVelocity);
        _maxU = _uVelocity[0][0];
        _maxV = _vVelocity[0][0];

        for(int i=0;i<_xsize;i++) {
            for(int j=0;j < _ysize;j++) {
                if(_maxU <= _uVelocity[i][j]) {
                    _maxU = _uVelocity[i][j];
                }
                if(_maxV <= _vVelocity[i][j]) {
                    _maxV = _vVelocity[i][j];
                }
            }
        }
        double max = sqrt(pow(_maxU, 2) +  pow(_maxV, 2));
        _deltaT = _dx/(4*max);
    }


    void setBoundaryConditions() {
    }

public:
    Advection(int xsize, int ysize, double dx, double dy,  FlagMatrix* flagMatrix, function < void(double** uMatrix, double** vMatrix) > fillVeloctyMatrixFnc):
        _xsize(xsize), _ysize(ysize), _flagMatrix(flagMatrix), _dx(dx), _dy(dy), _fillVelocityMatrixFnc(fillVeloctyMatrixFnc)  {
        _density0Matrix = imnd::matrix(_xsize, _ysize);
        _density1Matrix = imnd::matrix(_xsize, _ysize);
        _density2Matrix = imnd::matrix(_xsize, _ysize);

        _uVelocity = imnd::matrix(_xsize, _ysize);
        _vVelocity = imnd::matrix(_xsize, _ysize);
    }

    void Reset() {
        imnd::set_matrix(_density0Matrix, _xsize, _ysize, 0);
        imnd::set_matrix(_density1Matrix, _xsize, _ysize, 0);
        imnd::set_matrix(_density2Matrix, _xsize, _ysize, 0);
        imnd::set_matrix(_uVelocity, _xsize, _ysize, 0);
        imnd::set_matrix(_vVelocity, _xsize, _ysize, 0);
        //setBoundaryConditions();
        fillVelocityMatrixAndDeltaT();
        fillAdvectionMatrix();
        _nowT = 0;
        _iter = 0;
    }

    void NextTimestamp() {
        makeAdvection();
        //calculateIntegral();

        if(_iter++%int(100/_deltaT/50) == 0) {
            imnd::push_data2D(_density2Matrix, _xsize, _ysize);
        }
        imnd::copy_matrix(_density1Matrix, _density0Matrix,_xsize,_ysize);
        imnd::copy_matrix(_density2Matrix, _density1Matrix,_xsize,_ysize);
        _nowT+=_deltaT;

    }

    double GetTimestamp() {
        return _nowT;
    }

    void calculateIntegral() {

    }
    int GetIteration() {
        return _iter;
    }
    void SaveFlagsMatrixToPNGFile(const char *fileName) {
        imnd::plot_2d_system(fileName, _density0Matrix, _xsize, _ysize,  _dx, _dy);
    }

    virtual void SaveResults(char* filename) = 0;

    void PrintMatrixToFile(const char* fileName) {
        ofstream oFile;
        oFile.open (fileName);
        for(int i=0;i<_xsize;i++) {
            for(int j=0;j<_ysize;j++) {
                oFile << _density0Matrix[i][j] << ", ";
            }
            oFile << "\b\b\n";
        }
        oFile.close();
    }

};
#endif //IMN5_ADVECTION_H
