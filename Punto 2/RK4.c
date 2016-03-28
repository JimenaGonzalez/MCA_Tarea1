#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void RK4_step(double delta_t, double t, double *q1, double *p1,double *q3,double *p3,double e);
double deriv_q1(double t, double q1, double p1, double q3, double p3,double e);
double deriv_p1(double t, double q1, double p1, double q3, double p3,double e);
double deriv_q3(double t, double q1, double p1, double q3, double p3,double e);
double deriv_p3(double t, double q1, double p1, double q3, double p3,double e);
double E(double *q1, double *p1, double *q2, double *p2, double e);
int main(int argc, char ** argv){
  double q1_RK4;
  double q3_RK4; 
  double p1_RK4;
  double p3_RK4;
  double p2_RK4;
  double q2_RK4;
  double e;
  int i; 
  int j;
  double t;
  double energia; 
  double T=2800;
  double delta_t=0.006;
  int n_step =300000;//(int)(T/delta_t);
  e=0.681;
  t=0.0;
  
  q1_RK4 =strtod(argv[1],NULL);
  p1_RK4 =strtod(argv[2], NULL);
  q3_RK4 =strtod(argv[3], NULL);
  p3_RK4 =strtod(argv[4], NULL);
  for(i=0;i<n_step;i++){  

  
   if(p1_RK4<0.001 && p1_RK4>-0.001){ 
    printf("%.4e %.4e %.4e %.4e\n",q3_RK4, p3_RK4,t,energia);
   }
    q2_RK4=-q1_RK4;
    p2_RK4=-p1_RK4;
    energia=E( &q1_RK4, &p1_RK4, &q2_RK4, &p2_RK4, e);
  
    RK4_step(delta_t, t, &q1_RK4, &p1_RK4, &q3_RK4,&p3_RK4,e);
    
    t += delta_t;
  }

  return 0;
}

double deriv_q1(double t, double q1, double p1, double q3, double p3,double e){
  return p1;
}

double deriv_p1(double t, double q1, double p1, double q3, double p3,double e){
  return -2*q1/ pow( (4*pow(q1,2) + pow(e,2)) ,3/2 ) ;
}

double deriv_q3(double t, double q1, double p1, double q3, double p3,double e){
  return  p3;
  
}

double deriv_p3(double t, double q1, double p1, double q3, double p3,double e){
  return (q1-q3)/ pow( (pow((q1-q3),2) + pow(e,2.0)/4) ,3/2 ) - (q1+q3)/pow(  (pow((q1+q3),2.0) + pow(e,2.0)/4),  (3/2)   );
}

void RK4_step(double delta_t, double t, double *q1, double *p1,double *q3,double *p3,double e){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double w1, w2, w3, w4;
  double s1, s2, s3, s4;

  double q1_in;
  double p1_in;
  double q3_in;
  double p3_in;
  q1_in = *q1;
  p1_in = *p1;
  q3_in = *q3;
  p3_in = *p3;

//////// para q1 y p1///////////

  k1 = deriv_q1(t, q1_in, p1_in,  q3_in,  p3_in, e);  
  l1 = deriv_p1(t, q1_in, p1_in,  q3_in,  p3_in, e);  
               

  k2 = deriv_q1(t+ delta_t*0.5, q1_in+ k1*delta_t*0.5, p1_in+l1*delta_t*0.5,  q3_in+w1*delta_t*0.5,  p3_in+s1*delta_t*0.5, e);
  l2 = deriv_p1(t+ delta_t*0.5, q1_in+ k1*delta_t*0.5, p1_in+l1*delta_t*0.5,  q3_in+w1*delta_t*0.5,  p3_in+s1*delta_t*0.5, e) ;

  k3 = deriv_q1(t+ delta_t*0.5, q1_in+ k2*delta_t*0.5, p1_in+l2*delta_t*0.5,  q3_in+w2*delta_t*0.5,  p3_in+s2*delta_t*0.5, e);
  l3 = deriv_p1(t+ delta_t*0.5, q1_in+ k2*delta_t*0.5, p1_in+l2*delta_t*0.5,  q3_in+w2*delta_t*0.5,  p3_in+s2*delta_t*0.5, e);

  k4 = deriv_q1(t+ delta_t, q1_in+ k3*delta_t, p1_in+l3*delta_t,  q3_in+w3*delta_t,  p3_in+s3*delta_t, e);
  l4 = deriv_p1(t+ delta_t, q1_in+ k3*delta_t, p1_in+l3*delta_t,  q3_in+w3*delta_t,  p3_in+s3*delta_t, e);
  q1_in += (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  p1_in += (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*delta_t;
  
//////// para q3 y p3///////////

  w1 = deriv_q3(t, q1_in, p1_in,  q3_in,  p3_in, e);  
  s1 = deriv_p3(t, q1_in, p1_in,  q3_in,  p3_in, e);                 

  w2 = deriv_q3(t+ delta_t*0.5, q1_in+ k1*delta_t*0.5, p1_in+l1*delta_t*0.5,  q3_in+w1*delta_t*0.5,  p3_in+s1*delta_t*0.5, e);
  s2 = deriv_p3(t+ delta_t*0.5, q1_in+ k1*delta_t*0.5, p1_in+l1*delta_t*0.5,  q3_in+w1*delta_t*0.5,  p3_in+s1*delta_t*0.5, e) ;

  w3 = deriv_q3(t+ delta_t*0.5, q1_in+ k2*delta_t*0.5, p1_in+l2*delta_t*0.5,  q3_in+w2*delta_t*0.5,  p3_in+s2*delta_t*0.5, e);
  s3 = deriv_p3(t+ delta_t*0.5, q1_in+ k2*delta_t*0.5, p1_in+l2*delta_t*0.5,  q3_in+w2*delta_t*0.5,  p3_in+s2*delta_t*0.5, e);

  w4 = deriv_q3(t+ delta_t, q1_in+ k3*delta_t, p1_in+l3*delta_t,  q3_in+w3*delta_t,  p3_in+s3*delta_t, e);
  s4 = deriv_p3(t+ delta_t, q1_in+ k3*delta_t, p1_in+l3*delta_t,  q3_in+w3*delta_t,  p3_in+s3*delta_t, e);
 
  q3_in += (w1/6.0 + w2/3.0 + w3/3.0 + w4/6.0)*delta_t;
  p3_in += (s1/6.0 + s2/3.0 + s3/3.0 + s4/6.0)*delta_t;

  *q1 = q1_in;
  *p1 = p1_in;
  *q3 = q3_in;
  *p3 = p3_in;


}
double E(double *q1, double *p1, double *q2, double *p2, double e) {
  double p1_in;
  double q1_in;
  double p2_in;
  double q2_in;
  double K = 0;
  double U = 0;
  p1_in = *p1;
  q1_in = *q1;
  p2_in = *p2;
  q2_in = *q2;
  
  K += 0.5*p1_in*p1_in;
  K += 0.5*p2_in*p2_in;
  U += -1.0/(2.0*sqrt(4*q1_in*q1_in + e));
  U += -1.0/(2.0*sqrt(4*q2_in*q2_in + e));
return K + U;
}


