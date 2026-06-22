/***************************************************************************
 *
 * Boris-Massless-Integrator
 *
 * Classical fourth-order Runge-Kutta integrator used for comparison
 * with the Boris-type structure-preserving method.
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
 * This program integrates the equations of motion of a massless charged
 * particle using a standard explicit fourth-order Runge-Kutta scheme.
 * It is employed as a benchmark for assessing the long-time behavior
 * of the proposed Boris-type scheme and its preservation of the null
 * constraint
 *
 *      ε² - p² = 0
 *
 * and related invariants.
 *
 * Repository:
 *
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
#define N  1    /*numero de ecuaciones de primer orden*/
#define H  0.01   /* incremento de tiempo*/
#define Max 1.0e7  /*valor maximo del tiempo*/
#define q 1      /*carga*/


FILE *outputrk;      /*nombre del fichero de salida*/

double *RungeKutta4(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double step);    /*funcion runge kutta de cuarto orden*/
double funcion0(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);                      /*funcion para las derivadas*/
double funcion1(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
double funcion2(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
double funcion3(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
double funcion4(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
double funcion5(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
double funcion6(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
double funcion7(double L, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
double *puntero;	
double vec[9];
int main(void)  {

                        /*tiempo de ejecucion del programa*/
			time_t tiempo=time(0);
			struct tm *tlocal=localtime(&tiempo);
			char output[128];
			strftime(output,128,"%d/%m/%y %H:%M:%S",tlocal);
			printf("%s\n",output);

	double L,x0,x1,x2,x3,x4,x5,x6,x7;     
	              /*
	                t=lambda;
	                x0=tiempo; 
	                (x1,x2,x3) = (x,y,z);
	                (x4,x5,x6)=(px,py,pz) 
	                En=x7 
	                */
	int j;
		int flotante=0;		
	outputrk=fopen("outputrk.txt","w");
	
	L=0.0;
	vec[0]=0.0;
	vec[1]=0.0;       /*condiciones iniciales*/
	vec[2]=1.0;     
	vec[3]=0.0;
	vec[4]=1.0;
	vec[5]=0.0;
	vec[6]=0.0;
	vec[7]=1.5;
	
		for(j=0;j*H <= Max;j++) {

			L=j*H;
			puntero=RungeKutta4(L, vec[0], vec[1], vec[2], vec[3], vec[4], vec[5], vec[6], vec[7], H);
			vec[0]=puntero[0];
			vec[1]=puntero[1];
			vec[2]=puntero[2];
			vec[3]=puntero[3];
			vec[4]=puntero[4];
			vec[5]=puntero[5];
			vec[6]=puntero[6];
			vec[7]=puntero[7];
			vec[8]=puntero[8];
			/*
			double norm=sqrt(vec[4]*vec[4]+vec[5]*vec[5]+vec[6]*vec[6]);
			vec[4]*=vec[7]/norm;
			vec[5]*=vec[7]/norm;
			vec[6]*=vec[7]/norm;
			*/
			
			if(flotante==10000) /*loop para imprimir uno de cada 1000*/
				{
			    fprintf(outputrk,"%15f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
			    vec[8], 
			    vec[0], 
			    vec[1],
			    vec[2],
			    vec[3],
			    vec[4],
			    vec[5],
			    vec[6],
			    vec[7]);
                            
                            flotante=0;
			}

			flotante=flotante+1;
		}
fclose(outputrk);
printf("Programa finalizado \n");

                /*terminamos con el tiempo*/  
		time_t tiempo1=time(0);
		struct tm *tlocal1=localtime(&tiempo1);
		strftime(output,128,"%d/%m/%y %H:%M:%S",tlocal1);
		printf("%s\n",output);
return 0;
}


double *RungeKutta4(double L, double x0, double x1, double x2, double x3, double x4,double x5, double x6, double x7, double step){
		static double vec1[9];
double h=step/2.0; /* mitad de paso */
		
	double t11, t12, t13, t14, t15, t16, t17, t18;
	double t21, t22, t23, t24, t25, t26, t27, t28;
	double t31, t32, t33, t34, t35, t36, t37, t38;
	double t41, t42, t43, t44, t45, t46, t47, t48;
	
	double k11, k12, k13, k14, k15, k16, k17, k18;
	double k21, k22, k23, k24, k25, k26, k27, k28;
	double k31, k32, k33, k34, k35, k36, k37, k38;
	double k41, k42, k43, k44, k45, k46, k47, k48;    /* variables temporales para Runge-Kutta */

	t11= x0 + 0.5*(k11=step*funcion0(L, x0, x1, x2, x3, x4, x5, x6, x7));
	t12= x1 + 0.5*(k12=step*funcion1(L, x0, x1, x2, x3, x4, x5, x6, x7));
	t13= x2 + 0.5*(k13=step*funcion2(L, x0, x1, x2, x3, x4, x5, x6, x7));
	t14= x3 + 0.5*(k14=step*funcion3(L, x0, x1, x2, x3, x4, x5, x6, x7));
	t15= x4 + 0.5*(k15=step*funcion4(L, x0, x1, x2, x3, x4, x5, x6, x7));
	t16= x5 + 0.5*(k16=step*funcion5(L, x0, x1, x2, x3, x4, x5, x6, x7));
	t17= x6 + 0.5*(k17=step*funcion6(L, x0, x1, x2, x3, x4, x5, x6, x7));
	t18= x7 + 0.5*(k18=step*funcion7(L, x0, x1, x2, x3, x4, x5, x6, x7));
		
	t21= x0 + 0.5*(k21=step*funcion0(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	t22= x1 + 0.5*(k22=step*funcion1(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	t23= x2 + 0.5*(k23=step*funcion2(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	t24= x3 + 0.5*(k24=step*funcion3(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	t25= x4 + 0.5*(k25=step*funcion4(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	t26= x5 + 0.5*(k26=step*funcion5(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	t27= x6 + 0.5*(k27=step*funcion6(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	t28= x7 + 0.5*(k28=step*funcion7(L+h, t11, t12, t13, t14, t15, t16, t17, t18));
	
	t31= x0 + 0.5*(k31=step*funcion0(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	t32= x1 + 0.5*(k32=step*funcion1(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	t33= x2 + 0.5*(k33=step*funcion2(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	t34= x3 + 0.5*(k34=step*funcion3(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	t35= x4 + 0.5*(k35=step*funcion4(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	t36= x5 + 0.5*(k36=step*funcion5(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	t37= x6 + 0.5*(k37=step*funcion6(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	t38= x7 + 0.5*(k38=step*funcion7(L+h, t21, t22, t23, t24, t25, t26, t27, t28));
	
		
	k41 = step*funcion0(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	k42 = step*funcion1(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	k43 = step*funcion2(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	k44 = step*funcion3(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	k45 = step*funcion4(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	k46 = step*funcion5(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	k47 = step*funcion6(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	k48 = step*funcion7(L+step, t31, t32, t33, t34, t35, t36, t37, t38);
	
	x0 += (k11+2*k21+2*k31+k41)/6.0;
	x1 += (k12+2*k22+2*k32+k42)/6.0;
	x2 += (k13+2*k23+2*k33+k43)/6.0;
	x3 += (k14+2*k24+2*k34+k44)/6.0;
	x4 += (k15+2*k25+2*k35+k45)/6.0;
	x5 += (k16+2*k26+2*k36+k46)/6.0;
	x6 += (k17+2*k27+2*k37+k47)/6.0;
	x7 += (k18+2*k28+2*k38+k48)/6.0;
	
	vec1[0]=x0;
	vec1[1]=x1;
	vec1[2]=x2;
	vec1[3]=x3;
	vec1[4]=x4;
	vec1[5]=x5;
	vec1[6]=x6;
	vec1[7]=x7;
	vec1[8]=L;

	
return vec1;
/*free(vec1);*/
}


 

double funcion0(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6, double x7){

return(1);
}


double funcion1(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6, double x7){

return(x4/x7);
}


double funcion2(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6, double x7){

return(x5/x7);
}


double funcion3(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6, double x7){
return(x6/x7);
}

double funcion4(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6, double x7){
	double B1=0.0;   /*campo magnetico en la dirección e_x*/
	double B2=0.0;   /*campo magnetico en la dirección e_y*/
	double B3=1.0;   /*campo magnetico en la direccion e_z*/
	
	double Ex=0.2;       /*campo eléctrico en dirección e_x */
	double Ey=0.0; /*campo eléctrico en dirección e_y*/
	double Ez=0.0;       /*campo eléctrico en la dirección e_z*/
	
return(q*(Ex+x5*B3/x7-x6*B2/x7));
}

double funcion5(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6,double x7){
	double B1=0.0;   /*campo magnetico en la dirección e_x*/
	double B2=0.0;   /*campo magnetico en la dirección e_y*/
	double B3=1.0;   /*campo magnetico en la direccion e_z*/
	
	double Ex=0.2;        /*campo eléctrico en dirección e_x */
	double Ey=0.0;        /*campo eléctrico en dirección e_y*/
	double Ez=0.0;        /*campo eléctrico en la dirección e_z*/
return(q*(Ey+x6*B1/x7-x4*B3/x7));
}

double funcion6(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6,double x7){
	double B1=0.0;   /*campo magnetico en la dirección e_x*/
	double B2=0.0;   /*campo magnetico en la dirección e_y*/
	double B3=1.0;   /*campo magnetico en la direccion e_z*/
	
	double Ex=0.20;       /*campo eléctrico en dirección e_x */
	double Ey=0.0;       /*campo eléctrico en dirección e_y*/
	double Ez=0.0;       /*campo eléctrico en la dirección e_z*/
	
return(q*(Ez+x4*B2/x7-x5*B1/x7));
}

double funcion7(double L, double x0, double x1, double x2,double x3,double x4,double x5,double x6,double x7){
	double Ex=0.20;       /*campo eléctrico en dirección e_x */
	double Ey=0.0;       /*campo eléctrico en dirección e_y*/
	double Ez=0.0;       /*campo eléctrico en la dirección e_z*/

return(q*((x4/x7)*Ex+(x5/x7)*Ey+(x6/x7)*Ez));

}
/*http://www.sc.ehu.es/sbweb/fisica3/numerico/diferencial/diferencial.html
https://www.youtube.com/watch?v=rxAUwddWfQk*/
/*https://www.youtube.com/watch?v=XNIWRHGqMlA*/
/*  gcc -O3 -march=native -ffast-math RK4.c -lm */
