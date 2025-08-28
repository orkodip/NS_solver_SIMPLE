//SOLUTION OF ENERGY EQUATION FROM EXISTING FLUID FLOW SOLUTION.
//THE FLOW SOLUTION IS READ FROM FILES GENERATED AFTER RUNNING "jet.cpp"
//Refer to Mookherjee et al:
//"Numerical investigation of a confined laminar jet impingement cooling of heat sources using nanofluids."
//Journal of Heat Transfer 142.8 (2020): 082301

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include "cfd_solvers.h"

using namespace std;

const int M=60,N=492;	//M->No of grids in vertical direction, N->No of grids in horizontal direction

class PVC
{
	float Re, Pr, dt, n;
	float ***U, ***V, ***TH;	//temperature
	float A_W, A_N, A_S, A_E, A_P, A_P0, S;	//coeffecients of momentum and p-dash equation
	float **U_A_P;	//coeffecients of momentum equation required in pressure correction and residual calculation
	float **U_A_N, **U_A_E, **U_A_W, **U_A_S, **U_S, **U_A_P0;	//coeffecients of momentum equation required in residual calculation
	float **U_prev,**V_prev;	//not required but has to be declared
	float *a,*b,*c,*d;	//required for TDMA
	int cnt;	//iteration counter	
	float zeta_f;
	
	void TH_bc_top();	//boundary conditions for TEMPERATURE
	void TH_bc_right();
	void TH_bc_bottom();
	void TH_bc_left();

	void CD_geo_U(int j, int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float F_e, float F_w, float F_n, float F_s);
	void CD_U(int j,int i,float F_e,float F_w,float F_n,float F_s,float D_e,float D_w,float D_n,float D_s);	//grid independent Central Difference scheme
	void CD_geo_V(int j, int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float F_e, float F_w, float F_n, float F_s);
	void CD_V(int j,int i,float F_e,float F_w,float F_n,float F_s,float D_e,float D_w,float D_n,float D_s);

	void QUICK_geo_U(int j,int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float *B_1w,float *B_2w,float *B_1e,float *B_2e,float *B_1n,float *B_2n,float *B_1s,float *B_2s,float *q_w,float *q_e,float *q_n,float *q_s,float F_e, float F_w, float F_n, float F_s);
	void QUICK_U(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s);	//grid independent QUICK scheme	
	void QUICK_geo_V(int j,int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float *B_1w,float *B_2w,float *B_1e,float *B_2e,float *B_1n,float *B_2n,float *B_1s,float *B_2s,float *q_w,float *q_e,float *q_n,float *q_s,float F_e, float F_w, float F_n, float F_s);
	void QUICK_V(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s);
	void CD_geo_TH(int j, int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float F_e, float F_w, float F_n, float F_s);
	void CD_TH(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s);
	void QUICK_geo_TH(int j,int i,float *G_1w,float *G_3w,float *G_1e,float *G_3e,float *G_1n,float *G_3n,float *G_1s,float *G_3s,float *K_1w,float *K_3w,float *K_1e,float *K_3e,float *K_1n,float *K_3n,float *K_1s,float *K_3s,float *r_w,float *r_e,float *r_n,float *r_s,float F_e, float F_w, float F_n, float F_s);
	void QUICK_TH(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s);

	void coeff_Energy(int j, int i);
	void Energy();	//Energy equation
	float Energy_res(FILE *);

	public:
		PVC(); ~PVC();	//memory allocation and de-allocation
		void getdata(float re, float pr, float n1, float DT, float zeta);		//initialisation of the variables
		void Energy_solver();
		void write();	//write pressure and velocity contour in a single file		
};

PVC :: PVC()
{
	U=(float ***) malloc(3*sizeof(float **));	//3 because -> 0: x co-ordinate, 1: y co-ordinate, 2: value of variable
	V=(float ***) malloc(3*sizeof(float **));	
	TH=(float ***) malloc(3*sizeof(float **));	
	
	for (int i=0;i<3;i++)
	{
		U[i]=(float **) malloc((N+1)*sizeof(float *));	//allocation of the rows
		V[i]=(float **) malloc((N+1)*sizeof(float *));	
		TH[i]=(float **) malloc((N+1)*sizeof(float *));	
	}
	
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<=N;j++)
		{
			U[i][j]=(float *) malloc((M+1)*sizeof(float));	//allocation of the columns
			V[i][j]=(float *) malloc((M+1)*sizeof(float));	
			TH[i][j]=(float *) malloc((M+1)*sizeof(float));
		}
	}
		
	 U_A_P = (float **) malloc ((N+1)*sizeof(float *));	//allocation of the rows
	 U_A_N = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_E = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_W = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_S = (float **) malloc ((N+1)*sizeof(float *));
	 U_S = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_P0 = (float **) malloc ((N+1)*sizeof(float *));
	 	 
	 for (int i=0;i<=N;i++)
	 {
		U_A_P[i]=(float *) malloc ((M+1)*sizeof(float));	//allocation of the columns
		U_A_N[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_E[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_W[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_S[i]=(float *) malloc ((M+1)*sizeof(float)); 
		U_S[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_P0[i]=(float *) malloc ((M+1)*sizeof(float));	 
	}	

	a=(float *)malloc(N*sizeof(float));	//lower diagonal of TDMA
	b=(float *)malloc(N*sizeof(float));	//main diagonal of TDMA
	c=(float *)malloc(N*sizeof(float));	//upper diagonal of TDMA
	d=(float *)malloc(N*sizeof(float));	//constant of TDMA

	cout<<endl<<"MEMORY ALLOCATED"<<endl;
}

PVC :: ~PVC()
{	
	free(U); free(V); free(TH); free(U_A_P); free(U_A_N); free(U_A_E); free(U_A_W); free(U_A_S); free(U_S); free(U_A_P0);
	free(a); free(b); free(c); free(d);
	cout<<endl<<"MEMORY RELEASED"<<endl;
}

void PVC :: getdata(float re, float pr, float n1, float DT, float zeta)
{
	Re=re; Pr=pr; dt=DT; zeta_f=zeta; n=n1;
	A_W=0.0; A_N=0.0; A_S=0.0; A_E=0.0; A_P=0.0; A_P0=0.0; S=0.0; cnt=0;
	
	cout<<endl<<"Reynolds Number = "<<Re<<", Prandlt Number = "<<Pr<<endl;
	cout<<endl<<"Time-Step = "<<dt<<", zeta_f = "<<zeta_f<<endl;

	char buffer[100];
	FILE *p = fopen("U.dat","r");

	for (int a=1;a<=3;a++)	//skip first three lines of the data file
		fgets(buffer,100,p);

	for (int HH=0;HH<=2;HH++)
		for (int j=1;j<=M;j++)
			for (int i=1;i<N;i++)
				fscanf (p,"%f",&U[HH][i][j]);
	fclose(p);

	p=fopen("V.dat","r");

	for (int a=1;a<=3;a++)	//skip first three lines of the data file
		fgets(buffer,100,p);

	for (int HH=0;HH<=2;HH++)
		for (int j=1;j<M;j++)
			for (int i=1;i<=N;i++)
				fscanf (p,"%f",&V[HH][i][j]);
	fclose(p);

	p=fopen("P.dat","r");

	for (int a=1;a<=3;a++)	//skip first three lines of the data file
		fgets(buffer,100,p);

	for (int HH=0;HH<2;HH++)
		for (int j=1;j<=M;j++)
			for (int i=1;i<=N;i++)
				fscanf (p,"%f",&TH[HH][i][j]);	//read the mesh only
	fclose(p);
	
	TH_bc_left();
	TH_bc_right();
	TH_bc_top();	
	TH_bc_bottom();
}

void PVC :: TH_bc_top()	//NEUMANN
{
	for (int i=1;i<=N;i++)
	{
		if (TH[0][i][M]<=0.5)
			TH[2][i][M]=0.0;	//DIRICHLET (jet inlet temperature)
		else
			TH[2][i][M]=TH[2][i][M-1];	//NEUMANN (adiabetic)
	}
}
void PVC :: TH_bc_bottom()	//NEUMANN
{
	for (int i=1;i<=N;i++)
	{
		if (TH[0][i][1]<=1.5)	//1st heater which is of half length
			TH[2][i][1]=TH[2][i][2]+n*(TH[1][i][2]-TH[1][i][1]);
		else if ((TH[0][i][1]>=3.5)&&(TH[0][i][1]<=6.5))	//2nd heater
			TH[2][i][1]=TH[2][i][2]+n*(TH[1][i][2]-TH[1][i][1]);
		else if ((TH[0][i][1]>=8.5)&&(TH[0][i][1]<=11.5))	//3rd heater
			TH[2][i][1]=TH[2][i][2]+n*(TH[1][i][2]-TH[1][i][1]);
		else if ((TH[0][i][1]>=13.5)&&(TH[0][i][1]<=16.5))	//4th heater
			TH[2][i][1]=TH[2][i][2]+n*(TH[1][i][2]-TH[1][i][1]);
		else if ((TH[0][i][1]>=18.5)&&(TH[0][i][1]<=21.5))	//5th heater
			TH[2][i][1]=TH[2][i][2]+n*(TH[1][i][2]-TH[1][i][1]);
		else
			TH[2][i][1]=TH[2][i][2];	//adiabetic
	}
}
void PVC :: TH_bc_left()	//NEUMANN (symetry)
{
	for (int j=1;j<=M;j++)
		TH[2][1][j]=TH[2][2][j];
}
void PVC :: TH_bc_right()	//NEUMANN (outflow)
{
	for (int j=1;j<=M;j++)
		TH[2][N][j]=TH[2][N-1][j];
}
#include "GI_schemes_SOU_QUICK_en.cpp"	//Grid Independent schemes written in GI_schemes.cpp file

void PVC :: coeff_Energy(int j, int i)
{
	float F_e=0.0, F_w=0.0, F_n=0.0, F_s=0.0;
	float D_e=0.0, D_w=0.0, D_n=0.0, D_s=0.0;	
	float q_n=0.0, q_e=0.0, q_w=0.0, q_s=0.0;	//flux terms of energy equation are 0 inside the domain

	A_W=0.0; A_N=0.0; A_S=0.0; A_E=0.0; A_P=0.0; S=0.0; A_P0=0.0;	//re-initialization of the variables

	A_P0=(U[0][i][j]-U[0][i-1][j])*(V[1][i][j]-V[1][i][j-1])/dt;	//transient coeffecient

	F_e = (V[1][i][j]-V[1][i][j-1])*U[2][i][j];	//flux strengths	
	F_w = (V[1][i][j]-V[1][i][j-1])*U[2][i-1][j];
	F_n = (U[0][i][j]-U[0][i-1][j])*V[2][i][j];
	F_s = (U[0][i][j]-U[0][i-1][j])*V[2][i][j-1];

	D_e = zeta_f*(V[1][i][j]-V[1][i][j-1])/((TH[0][i+1][j]-TH[0][i][j])*Re*Pr);	//diffusion coeffecients
	D_w = zeta_f*(V[1][i][j]-V[1][i][j-1])/((TH[0][i][j]-TH[0][i-1][j])*Re*Pr);
	D_n = zeta_f*(U[0][i][j]-U[0][i-1][j])/((TH[1][i][j+1]-TH[1][i][j])*Re*Pr);
	D_s = zeta_f*(U[0][i][j]-U[0][i-1][j])/((TH[1][i][j]-TH[1][i][j-1])*Re*Pr);

	if ((j==M-1)||(i==2)||(j==2)||(i==N-1))
		CD_TH(j,i,F_e,F_w,F_n,F_s,D_e,D_w,D_n,D_s);	//CENTRAL DIFFERENCE SCHEME FOR THE BORDER CV's
	else
		QUICK_TH(j,i,F_e,F_w,F_n,F_s,D_e,D_w,D_n,D_s);	//QUICK SCHEME

	A_P+=A_P0;	//addition of the psuedo-transient term
	
	U_A_P[i][j]=A_P;	//storing the coeffecients for use in residual calculation
	U_A_N[i][j]=A_N;	//VARIABLES ARE RE-USED HERE !!
	U_A_E[i][j]=A_E;	//MIGHT RESULT IN CONFLICT IN FUTURE !!
	U_A_W[i][j]=A_W;	//ERROR PRONE AREA !!
	U_A_S[i][j]=A_S;
	U_S[i][j]=S;
	U_A_P0[i][j]=A_P0;	
}

void PVC :: Energy()
{
	for (int z=M-1;z>=2;z--)	//the jacobian matrix changes for every row z=row
	{	
		//1st row of the Jacobian matrix
		coeff_Energy(z,2);
		a[0]=0.0;
		b[0]=A_P;
		c[0]=-A_E;		
		d[0]=S+A_W*TH[2][1][z]+A_N*TH[2][2][z+1]+A_S*TH[2][2][z-1]+A_P0*TH[2][2][z];
	
		for (int i=1;i<N-3;i++)	//generation of the coeff mat i=column of b matrix
		{
			coeff_Energy(z,i+2);
			c[i]=-A_E;
			b[i]=A_P;
			a[i]=-A_W;
			d[i]=S+A_N*TH[2][i+2][z+1]+A_S*TH[2][i+2][z-1]+A_P0*TH[2][i+2][z];
		}
	
		//last row of the Jacobian matrix
		coeff_Energy(z,N-1);
		a[N-3]=-A_W;
		b[N-3]=A_P;
		c[N-3]=0.0;
		d[N-3]=S+A_E*TH[2][N][z]+A_N*TH[2][N-1][z+1]+A_S*TH[2][N-1][z-1]+A_P0*TH[2][N-1][z];
		
		tdma_bs(N-2,a,b,c,d);	//TDMA with Backward substitution present in "cfd_solvers.h" header file
		//Theta value for the particular row is found out
		for (int i=2;i<=N-1;i++)
			TH[2][i][z]=a[i-2];	//updation of the temperature 
	}
}

void PVC :: Energy_solver()
{
	FILE *p;	//extracting temperature residuals
	p=fopen("temp_res.dat","w");	//generate temp_res.dat file
	fprintf (p,"TITLE = \"TEMPERATURE RESIDUAL PLOT\"\n");
	fprintf (p,"VARIABLES = \"Iteration\", \"Residual*10^7\"\n");
	cnt=0;	//re-initialization after flow field solution
	do
	{
		cnt++;		
		if ((cnt%1000)==0)
			cout<<endl<<"ITERATION "<<cnt<<endl;

		if ((cnt%10000)==0)	//export file for every specified iterations
		{
			write();
			cout<<endl<<"INTERMEDIATE FILE OUTPUT SUCCESSFUL"<<endl;
		}
		Energy();

		TH_bc_left();	//neumann boundary update
		TH_bc_right();
		TH_bc_top();
		TH_bc_bottom();
	}
	while (Energy_res(p) > 18e-7);

	fclose(p);
	cout<<endl<<"TEMPERATURE RESIDUAL FILE OUTPUT SUCCESSFUL"<<endl;
}

float PVC :: Energy_res(FILE *p)
{
	double cc=0.0;	//residual values
	
	for (int j=M-1;j>1;j--)
		for (int i=2;i<=N-1;i++)
			cc+=pow((-U_A_P[i][j]*TH[2][i][j]+U_A_E[i][j]*TH[2][i+1][j]+U_A_W[i][j]*TH[2][i-1][j]+U_A_N[i][j]*TH[2][i][j+1]+U_A_S[i][j]*TH[2][i][j-1]+U_A_P0[i][j]*TH[2][i][j]+U_S[i][j]),2.0);

	cc=sqrt(cc);

	fprintf(p,"%d\t%lf\n",cnt,cc*1e7);

	if ((cnt%30)==0)
		cout<<"ENERGY RESIDUAL="<<cc<<endl;
	
	return (cc);
}

void PVC :: write()
{
	FILE *p; float temp;
	p=fopen("temp.dat","w");	//generate vel_press.dat file
	fprintf (p,"TITLE = \"TEMPERATURE FIELD\"\n");	
	fprintf (p,"VARIABLES = \"X\", \"Y\", \"TH\"\n");	
	fprintf (p,"ZONE I=%d,J=%d,DATAPACKING=BLOCK\n",N,M);

	for (int j=1;j<=M;j++)	//print X co-ordinates
	{
		for (int i=1;i<=N;i++)
			fprintf (p," %.4f",TH[0][i][j]);
			
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");
	
	for (int j=1;j<=M;j++)	//print Y co-ordinates
	{
		for (int i=1;i<=N;i++)
			fprintf (p," %.4f",TH[1][i][j]);
			
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");

	for (int j=1;j<=M;j++)	//print temperature
	{
		for (int i=1;i<=N;i++)
			fprintf (p," %.4f",TH[2][i][j]);
			
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");

	fclose(p);
	cout<<endl<<"FILE OUTPUT SUCCESSFUL"<<endl;
}

int main()
{
	clock_t start = clock();
	
	PVC ms;

	ms.getdata (100.0,6.95599,1.0,0.01,1.06336);
	//getdata(float re, float pr, float n, float DT, float zeta)
	ms.Energy_solver();	//energy equation not coupled with momentum equation. (forced convection)
	ms.write();
		
	cout<<endl<<"END OF SIMULATION"<<endl;
	
	clock_t end = clock();
	printf("\nElapsed Time : %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
	return(0);
}
