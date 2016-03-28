#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void leapfrog_step(double delta_t, double t, double *q1, double *p1,double *q3,double *p3, double e);
double deriv_p1(double t, double q1, double p1, double q3, double p3,double e);
double deriv_q3(double t, double q1, double p1, double q3, double p3,double e);
double deriv_p3(double t, double q1, double p1, double q3, double p3,double e);
double E(double *q1, double *p1, double *q2, double *p2, double e);

int main(int argc, char ** argv){
  double q1_LF;
  double p1_LF;
  double q3_LF;
  double p3_LF;
  double q2_LF;
  double p2_LF;
  double e;
  int i; 
  int j;
  double t;

  double T=2800;
  double delta_t=0.006;
  int n_step =300000;//(int)(T/delta_t);
  e=0.681;
  
  t=0.0;  

  double energia; 
  q1_LF =strtod(argv[1],NULL);
  p1_LF =strtod(argv[2], NULL);
  q3_LF =strtod(argv[3], NULL);
  p3_LF =strtod(argv[4], NULL);
  
  
  for(i=0;i<n_step;i++){  

   if(p1_LF<0.001 && p1_LF>-0.001){ 
    printf("%.4e %.4e %.4e %.4e\n",q3_LF, p3_LF,t,energia);
   }
    q2_LF=-q1_LF;
    p2_LF=-p1_LF;
    energia=E( &q1_LF, &p1_LF, &q2_LF, &p2_LF, e);
    leapfrog_step(delta_t, t, &q1_LF, &p1_LF, &q3_LF,&p3_LF,e);
    
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

void leapfrog_step(double delta_t, double t, double *q1, double *p1,double *q3,double *p3,double e){
  double p1_in;
  double q1_in;
  double p3_in;
  double q3_in;
  p1_in = *p1;
  q1_in = *q1;
  p3_in = *p3;
  q3_in = *q3;
  double theta;
  theta=1/(2-pow(2,1/3));
  
/////para q1 y p1///////


  /*drift*/
   q1_in += theta*p1_in*delta_t;

  /*kick*/
   p1_in += theta * deriv_p1(t, q1_in,p1_in, q3_in,p3_in, e) * delta_t;
 
  /*drift*/
   q1_in +=  0.5*(1-theta)*p1_in * delta_t;
  
  /*kick*/
   p1_in += (1-2*theta) * deriv_p1(t, q1_in,p1_in, q3_in,  p3_in, e) * delta_t;
  
  /*drift*/
   q1_in +=  0.5*(1-theta)*p1_in * delta_t; 

  /*kick*/
   p1_in += theta * deriv_p1(t, q1_in,p1_in, q3_in,p3_in, e) * delta_t;

  /*drift*/
   q1_in += theta*p1_in*delta_t;



//////para q3 y p3//////////


  /*drift*/
   q3_in += theta*p3_in*delta_t;

  /*kick*/
   p3_in += theta * deriv_p3(t, q1_in,p1_in, q3_in,p3_in, e) * delta_t;
 
  /*drift*/
   q3_in +=  0.5*(1-theta)*p3_in * delta_t;
  
  /*kick*/
   p3_in += (1-2*theta) * deriv_p3(t, q1_in,p1_in, q3_in,  p3_in, e) * delta_t;
  
  /*drift*/
   q3_in +=  0.5*(1-theta)*p3_in * delta_t; 

  /*kick*/
   p3_in += theta * deriv_p3(t, q1_in,p1_in, q3_in,p3_in, e) * delta_t;

  /*drift*/
   q3_in += theta*p3_in*delta_t;


 



  *p1 = p1_in;
  *q1 = q1_in;
  *p3 = p3_in;
  *q3 = q3_in;
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




