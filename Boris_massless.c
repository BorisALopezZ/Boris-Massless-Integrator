/***************************************************************************
 *
 * Boris-Massless-Integrator
 *
 * Structure-preserving Boris-type integrator for massless charged
 * particle dynamics.
 *
 * Author:
 * Boris Alberto Lopez Zubieta
 *
 * Manuscript:
 * "A Structure-preserving Boris-type integrator for massless charged
 * particle dynamics"
 * Submitted to Computer Physics Communications (2026).
 *
 * Description:
 * This program implements a Boris-type structure-preserving integrator
 * for massless charged particles propagating in prescribed
 * electromagnetic fields. The algorithm preserves the null constraint
 *
 *      ε² - p² = 0
 *
 * through a geometric formulation on the unit sphere S².
 * The implementation reproduces the numerical examples presented
 * in the accompanying manuscript.
 *
 * Repository:
 * https://github.com/BorisALopezZ/Boris-Massless-Integrator
 *
 * License:
 * MIT License
 *
 * Version:
 * 1.0
 *
 * Date:
 * June 2026
 *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


void BORIS_MASSLESS(
double n_old[3],
double x_old[3],
double epsilon,
double E[3],
double B[3],
double Dt,
double q,
double n_new[3],
double x_new[3],
double *epsilon_new)
{

    double n_minus[3], n_plus[3];
    double t[3], s[3];
    double cross1[3], cross2[3];
    double ndotE = n_old[0]*E[0]+n_old[1]*E[1]+n_old[2]*E[2];

    double eps_minus = epsilon + 0.5*q*Dt*ndotE;

    /* electric half step */

      for(int i=0;i<3;i++)
        /*n_minus[i]=n_old[i]+0.5*q*Dt/epsilon*(E[i]-ndotE*n_old[i]);  */
        n_minus[i]=n_old[i]+0.5*q*Dt/eps_minus*(E[i]-ndotE*n_old[i]);

    /* magnetic rotation */

      for(int i=0;i<3;i++)
        t[i]=q*B[i]*Dt/(2.0*epsilon);

    double T=t[0]*t[0]+t[1]*t[1]+t[2]*t[2];

    for(int i=0;i<3;i++)
      s[i]=2*t[i]/(1+T);

    /* cross products */

    cross1[0]=n_minus[1]*t[2]-n_minus[2]*t[1];
    cross1[1]=n_minus[2]*t[0]-n_minus[0]*t[2];
    cross1[2]=n_minus[0]*t[1]-n_minus[1]*t[0];

    double n_prime[3];

      for(int i=0;i<3;i++)
        n_prime[i]=n_minus[i]+cross1[i];

    cross2[0]=n_prime[1]*s[2]-n_prime[2]*s[1];
    cross2[1]=n_prime[2]*s[0]-n_prime[0]*s[2];
    cross2[2]=n_prime[0]*s[1]-n_prime[1]*s[0];

    for(int i=0;i<3;i++)
    n_plus[i]=n_minus[i]+cross2[i];

    /* second electric step */

    double ndotE2=n_plus[0]*E[0]+n_plus[1]*E[1]+n_plus[2]*E[2];

    *epsilon_new=eps_minus+0.5*q*Dt*ndotE2;

    for(int i=0;i<3;i++)
    n_new[i]=n_plus[i]+0.5*q*Dt/eps_minus*(E[i]-ndotE2*n_plus[i]);

    /* renormalize */

    double norm=sqrt(n_new[0]*n_new[0]+n_new[1]*n_new[1]+n_new[2]*n_new[2]);

      for(int i=0;i<3;i++)
        n_new[i]/=norm;

    /* position */
    
         for(int i=0;i<3;i++)
          x_new[i]=x_old[i] + 0.5*Dt*n_old[i]+0.5*Dt*n_new[i];
          //x_new[i]=x_old[i]+Dt*n_new[i];
}

FILE *output;      /*nombre del fichero de salida*/

int main(void){
                  /*tiempo de ejecucion del programa*/
			time_t tiempo=time(0);
			struct tm *tlocal=localtime(&tiempo);
			char output[128];
			strftime(output,128,"%d/%m/%y %H:%M:%S",tlocal);
			printf("%s\n",output);

	/*inicializamos la posición inicial*/   double x_old[3]={0,1,0};  // X, Y, Z
	/*inicializamos la posición final*/	double x_new[3];
	
	/*inicializamos el campo electrico*/    double E[3];
	/*inicializamos el campo magnetico*/    double B[3];
	/*vector potencial magnetico*/	        double A[3];
	                                               A[2]=0;
	/*inicializamos el momento final*/      double n_new[3];
	/*energia*/                             double epsilon=1.5;
	/*energia final*/                       double epsilon_new;      

	                                      double n_old[3]={1,0,0};    //Px, Py, Pz
	
	/*ponemos el delta "temporal"*/         double Dt=0.01;//M_PI/10;   /*el valor de Dt minimo optimo es 0.000028, con valores menores ya salen complejos*/
	/*hasta donde queremos que vaya*/       double max=1.0e7;//50000.0;    //max=15.0;
	/*phi*/				        double phi,q0,q3;
	
	/*carga q*/				double q=1;
	
	/*Campo eléctrico*/     E[0]=0.2;
				E[1]=0.0;
				E[2]=0.0;		      

    	outputrk=fopen("mass_less3.txt","w");
    	
    			int flotante=0;		
    			
          long long N = (long long)(max/Dt);
          for(long long j=0;j<N;j++){
	
	               B[0] = 0;
                       B[1] = 0;
                       B[2] = 1;
                       
                       BORIS_MASSLESS(n_old,x_old,epsilon,E,B,Dt,q, n_new, x_new, &epsilon_new);
                       if(flotante==10000) /*loop para imprimir uno de cada 1000*/
                       {
                       //double cuadradoB=B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
                       double moduloB=1.0; //sqrt(cuadradoB);
                       double adiabatica=epsilon*epsilon*(1-n_new[2])*(1-n_new[2])/moduloB;
                       double norm = sqrt(n_new[0]*n_new[0]+n_new[1]*n_new[1]+n_new[2]*n_new[2]);
                      
                       /*energia*/
			double phy=-E[0]*x_new[0]-E[1]*x_new[1]-E[2]*x_new[2];  /*ojo con el signo*/
			double cuak=epsilon_new*n_new[0]*epsilon_new*n_new[0]+epsilon_new*n_new[1]*epsilon_new*n_new[1]+epsilon_new*n_new[2]*epsilon_new*n_new[2];
			double raiz=sqrt(cuak);
			double energy=raiz+q*phy;
			
			/*momento angular*/
			double vec_cua=x_new[0]*x_new[0]+x_new[1]*x_new[1]+x_new[2]*x_new[2];
			double cuadrado=sqrt(vec_cua);
			double resultado=cuadrado*cuadrado*cuadrado;
			A[0]=-x_new[1]/resultado;
			A[1]=x_new[0]/resultado;
			
			/*Momento Dipolar*/
			q3=x_new[0]*(epsilon_new*n_new[1]+q*A[1])-x_new[1]*(epsilon_new*n_new[0]+q*A[0]);  

                      
                     /*
                     lista de componentes (Columns):

                            1 time
                            2 px
                            3 py
                            4 pz
                            5 x
                            6 y
                            7 z
                            8 mu
                            9 |B|
                            10 epsilon
                            11 |n|
                            12 energy
                            13 canonical angular momentum
                      */
                      fprintf(outputrk,"%f\t%f\t%lf\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                         j*Dt,
                         epsilon_new*n_new[0],
                         epsilon_new*n_new[1],
                         epsilon_new*n_new[2],
                         x_new[0],
                         x_new[1],
                         x_new[2],
                         adiabatica,
                         moduloB,
                         epsilon_new,
                         norm,
                         energy,
                         q3);
                         
                         flotante=0;		
                       }
                       flotante=flotante+1;
                       //epsilon_new = epsilon;
                       epsilon=epsilon_new;
                       n_old[0]=n_new[0];
                       n_old[1]=n_new[1];
                       n_old[2]=n_new[2];
                       
                       x_old[0]=x_new[0];
                       x_old[1]=x_new[1];
                       x_old[2]=x_new[2];
                       
                       if(!isfinite(x_new[0]) || !isfinite(epsilon_new)){
                       printf("NaN detectado en t=%f\n", j*Dt);
                       break;
                        }		      	
        }
        
       fclose(outputrk);	
		/*terminamos con el tiempo*/
		time_t tiempo1=time(0);
		struct tm *tlocal1=localtime(&tiempo1);
		strftime(output,128,"%d/%m/%y %H:%M:%S",tlocal1);
		printf("%s\n",output);
	
	return 0; 
}
/*  gcc -O3 -march=native Boris_massless.c -lm -o boris_massless */

/* gcc -O3 -march=native -ffast-math Boris_massless.c -lm -o boris_massless  */
