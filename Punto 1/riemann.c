#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "riemann.h"
#define MIN(x,y) x<y?x:y

void Riemann(double *U1, double *U4, double *F) {

  const double gamma = 1.4;
  
  // variables densidad
  double rho1 = U1[0];
  double rho4 = U4[0];

  // variable primitivas (u, p, h)
  double u1 = U1[1] / rho1;
  double p1 = (gamma - 1.0) * (U1[2] - rho1 * u1 * u1 / 2.0);
  double h1 = (U1[2] + p1) / rho1;
  double u4 = U4[1] / rho4;
  double p4 = (gamma - 1.0) * (U4[2] - rho4 * u4 * u4 / 2.0);
  double h4 = (U4[2] + p4) / rho4;

  // Cantidades promedios de Roe rho14, u14, h14, a14

  double rho14 = sqrt(rho1 * rho4);
  double u14 = ( sqrt(rho4) * u4 + sqrt(rho1) * u1 ) / ( sqrt(rho4) + sqrt(rho1) );
  double h14 = ( sqrt(rho4) * h4 + sqrt(rho1) * h1 ) / ( sqrt(rho4) + sqrt(rho1) );
  double a14 = sqrt ( (gamma - 1.0) * (h14 - u14 * u14 / 2.0) );

  // Velocidades promedio de Roe lam1, lam2, lam3  

  double lam1 = u14;
  double lam2 = u14 + a14;
  double lam3 = u14 - a14;

  // deltas drho, dp, du  &  wave strengths dv1, dv2, dv3

  double drho = rho4 - rho1;
  double dp = p4 - p1;
  double du = u4 - u1;

  double dv1 = drho - dp / (a14 * a14) ;
  double dv2 = du + dp / (rho14 * a14) ;
  double dv3 = du - dp / (rho14 * a14) ;

  // Expresiones para el flujo
  
 F[0] = rho1 * u1 + (MIN(0.0,lam1)) * dv1 + (rho14 * (MIN(0.0,lam2)) * dv2) / (2.0 * a14) - (MIN(0.0,lam3)) * dv3;

  F[1] = rho1 * u1 * u1 + p1 + u14 * (MIN(0.0,lam1)) * dv1 + 
         (lam2 * (MIN(0.0,lam2)) * dv2 - lam3 * (MIN(0.0, lam3)) * dv3) * rho14 / (2.0 * a14); 

  F[2] = rho1 * h1 * u1 + u14 * u14 * (MIN(0.0, lam1)) * dv1 / 2.0 + 
    ((h14 + a14 * u14) * (MIN(0.0,lam2)) * dv2 - (h14 - a14 * u14) * (MIN(0.0,lam3)) * dv3 ) * rho14 / (2.0 * a14);
  
}
