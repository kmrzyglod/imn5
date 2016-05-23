#include <iostream>
#include "imnmath.h"
#include "Advection.h"
#include "AdvectionLeapFrog.h"

using namespace std;
#define DX 0.02
#define DY 0.02
#define YMIN 0
#define XMIN 0

void fillVelocityLeapFrog(double** uVelocity, double** vVelocity) {
    int iMax = 301;
    int jMax = 41;
    double Q = -1;
    double mi = 1;


    for (int i=0;i<iMax;i++) {
        for (int j=0;j<jMax;j++) {
            uVelocity[i][j] = Q/(2*mi)*(DY * j - YMIN)*(DY * j - jMax*DY);
            vVelocity[i][j] = XMIN;
        }
    }
}

void fillVelocityFromFile(double** uVelocity, double** vVelocity) {
    int iMax = 301;
    int jMax = 41;

    imnd::uv_load(uVelocity, vVelocity,iMax,jMax,"predkosc.dat");
}


void zad1() {
    FlagMatrix flagMatrix = FlagMatrix(301, 41);
    flagMatrix.SetBorders();
    Advection * advection = new AdvectionLeapFrog(301, 41, DX, DY, &flagMatrix, fillVelocityLeapFrog);
    advection->Reset();


    while(advection->GetTimestamp() <= 100) {
        advection->NextTimestamp();
    }
    advection->SaveResults("zad1_rho.txt",  "zad1_integral.txt", "zad1_center.txt");

}

void zad2() {
    FlagMatrix flagMatrix = FlagMatrix(301, 41);
    flagMatrix.SetBorders();
    Advection * advection = new AdvectionLeapFrog(301, 41, DX, DY, &flagMatrix, fillVelocityFromFile);
    advection->Reset();

    while(advection->GetTimestamp() <= 100) {
        advection->NextTimestamp();
    }
    advection->SaveResults("zad2_rho.txt", "zad2_integral.txt", "zad2_center.txt");

}
int main() {
    zad1();
    zad2();
    system("gnuplot plot.plt");
    return 0;
}