//IMPRINGING JET WITH NANOFLUID IN FVM USING SIMPLE ALGORITHM AND A COMBINED SCHEME OF CENTRAL DIFFERENCE AND QUICK
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
const float LENGTH=50.5, HEIGHT=4.0;	//size of domain

class PVC
{
	float Re,dt,Gp;	//Reynolds number, time step, convergence criteria
	float ***U, ***V, ***P;	//x-velocity, y-velocity and pressure
	float **Pd;	//pressure dashed
	float **U_prev, **V_prev;
	float A_W, A_N, A_S, A_E, A_P, A_P0, S;	//coeffecients of momentum and p-dash equation
	float **U_A_P, **V_A_P;	//coeffecients of momentum equation required in pressure correction and residual calculation
	float **U_A_N, **U_A_E, **U_A_W, **U_A_S, **U_S, **U_A_P0;	//coeffecients of momentum equation required in residual calculation
	float **V_A_N, **V_A_E, **V_A_W, **V_A_S, **V_S, **V_A_P0;
	float *a,*b,*c,*d;	//required for TDMA
	float rex_P, theta;	//relaxation parameter and stone's cancellation parameter used in pressure correction equation.
	float V_in,U_out;	//average velocity at inlet and outlet
	int cnt;	//iteration counter	
	float r_f,eta_f;	//density and viscocity ratio

	void grid_gen();	//grid generator

	void U_bc_top();	//boundary conditions for X-velocity
	void U_bc_right();
	void U_bc_bottom();
	void U_bc_left();

	void V_bc_top();	//boundary conditions for Y-velocity
	void V_bc_right();
	void V_bc_bottom();
	void V_bc_left();

	void P_bc_top();	//boundary conditions for PRESSURE
	void P_bc_right();
	void P_bc_bottom();
	void P_bc_left();

	void CD_geo_U(int j, int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float F_e, float F_w, float F_n, float F_s);
	void CD_U(int j,int i,float F_e,float F_w,float F_n,float F_s,float D_e,float D_w,float D_n,float D_s);	//grid independent Central Difference scheme
	void CD_geo_V(int j, int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float F_e, float F_w, float F_n, float F_s);
	void CD_V(int j,int i,float F_e,float F_w,float F_n,float F_s,float D_e,float D_w,float D_n,float D_s);

	void QUICK_geo_U(int j,int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float *B_1w,float *B_2w,float *B_1e,float *B_2e,float *B_1n,float *B_2n,float *B_1s,float *B_2s,float *q_w,float *q_e,float *q_n,float *q_s,float F_e, float F_w, float F_n, float F_s);
	void QUICK_U(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s);	//grid independent QUICK scheme	
	void QUICK_geo_V(int j,int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float *B_1w,float *B_2w,float *B_1e,float *B_2e,float *B_1n,float *B_2n,float *B_1s,float *B_2s,float *q_w,float *q_e,float *q_n,float *q_s,float F_e, float F_w, float F_n, float F_s);
	void QUICK_V(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s);

	void coeff_U_mom(int j, int i);	//(row=j,column=i ) equation coeffecients calculation
	void coeff_V_mom(int j, int i);
	void coeff_P_dash(int j, int i);

	void U_mom();	//Navier-Stokes Momentum equations
	void V_mom();

	void P_dash();	//modified continuity equation for pressure correction
		
	float cont_source(FILE *);	//calculation of the continuity source term

	void P_corr();	//correction equations
	void U_corr();
	void V_corr();
		
	void U_mom_solver();
	void V_mom_solver();
	void P_dash_solver();	

	float U_mom_res();	//equation residuals
	float V_mom_res();
	float P_dash_res(int);
		
	public:		
		PVC();	//memory allocation via constructor
		~PVC();	//destructor used for freeing allocated memory
		
		void getdata(float re, float rexp, float r, float eta);		//initialisation of the variables
		void SIMPLE();	//SIMPLE algorithm
		void write();	//write pressure and velocity contour in a single file		
};

PVC :: PVC()
{
	U=(float ***) malloc(3*sizeof(float **));	//3 because -> 0: x co-ordinate, 1: y co-ordinate, 2: value of variable
	V=(float ***) malloc(3*sizeof(float **));	
	P=(float ***) malloc(3*sizeof(float **));	
	
	for (int i=0;i<3;i++)
	{
		U[i]=(float **) malloc((N+1)*sizeof(float *));	//allocation of the rows
		V[i]=(float **) malloc((N+1)*sizeof(float *));	
		P[i]=(float **) malloc((N+1)*sizeof(float *));	
	}
	
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<=N;j++)
		{
			U[i][j]=(float *) malloc((M+1)*sizeof(float));	//allocation of the columns
			V[i][j]=(float *) malloc((M+1)*sizeof(float));	
			P[i][j]=(float *) malloc((M+1)*sizeof(float));	
		}
	}
		
	 U_A_P = (float **) malloc ((N+1)*sizeof(float *));	//allocation of the rows
	 U_A_N = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_E = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_W = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_S = (float **) malloc ((N+1)*sizeof(float *));
	 U_S = (float **) malloc ((N+1)*sizeof(float *));
	 U_A_P0 = (float **) malloc ((N+1)*sizeof(float *));
	 V_A_P = (float **) malloc ((N+1)*sizeof(float *));
	 V_A_N = (float **) malloc ((N+1)*sizeof(float *));
	 V_A_E = (float **) malloc ((N+1)*sizeof(float *));
	 V_A_W = (float **) malloc ((N+1)*sizeof(float *));
	 V_A_S = (float **) malloc ((N+1)*sizeof(float *));
	 V_S = (float **) malloc ((N+1)*sizeof(float *));
	 V_A_P0 = (float **) malloc ((N+1)*sizeof(float *));
	 U_prev = (float **) malloc ((N+1)*sizeof(float *));
	 V_prev = (float **) malloc ((N+1)*sizeof(float *));
	 Pd = (float **) malloc ((N+1)*sizeof(float *));
	 
	 for (int i=0;i<=N;i++)
	 {
		U_A_P[i]=(float *) malloc ((M+1)*sizeof(float));	//allocation of the columns
		U_A_N[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_E[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_W[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_S[i]=(float *) malloc ((M+1)*sizeof(float)); 
		U_S[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_A_P0[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_A_P[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_A_N[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_A_E[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_A_W[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_A_S[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_S[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_A_P0[i]=(float *) malloc ((M+1)*sizeof(float));	 
		U_prev[i]=(float *) malloc ((M+1)*sizeof(float));	 
		V_prev[i]=(float *) malloc ((M+1)*sizeof(float));	 
		Pd[i]=(float *) malloc ((M+1)*sizeof(float));
	}

	a=(float *)malloc(N*sizeof(float));	//lower diagonal of TDMA
	b=(float *)malloc(N*sizeof(float));	//main diagonal of TDMA
	c=(float *)malloc(N*sizeof(float));	//upper diagonal of TDMA
	d=(float *)malloc(N*sizeof(float));	//constant of TDMA
	
	cout<<endl<<"MEMORY ALLOCATED"<<endl;
	grid_gen();	//generation of grid done here
}

PVC :: ~PVC()
{	
	free(U); free(V); free(P); free(Pd); free(U_prev); free(V_prev); free(U_A_P); free(U_A_N); free(U_A_E); free(U_A_W); free(U_A_S); free(U_S); free(U_A_P0);
	free(V_A_P); free(V_A_N); free(V_A_E); free(V_A_W); free(V_A_S); free(V_S); free(V_A_P0);
	cout<<endl<<"MEMORY RELEASED"<<endl;
}

void PVC :: getdata(float re, float rexp, float r, float eta)
{
	Re=re; rex_P=rexp; dt=0.005; Gp=0.2; theta=1.5;	//SOME PARAMETERS ARE HARD CODED HERE
	r_f=r; eta_f=eta;
	A_W=0.0; A_N=0.0; A_S=0.0; A_E=0.0; A_P=0.0; A_P0=0.0; S=0.0; cnt=0;
	
	cout<<endl<<"Reynolds Number = "<<Re<<", Pressure Relaxation = "<<rex_P<<endl<<endl;
	cout<<endl<<"R_f = "<<r_f<<", eta_f = "<<eta_f<<endl;
	
	U_bc_left();	//order of setting the boundary condition is important
	U_bc_right();
	U_bc_top();	
	U_bc_bottom();

	V_bc_left();
	V_bc_right();
	V_bc_top();	
	V_bc_bottom();
}

void PVC :: grid_gen()
{
	float dx[N+1],dy[M+1],y=0.0,factor_X1=1.001,factor_X2=1.006,factor_Y=1.025;

	for (int i=2;i<=12;i++)
		dx[i] = (0.5)/(12-1);	//12 grids in the jet region

	dx[13]=10.0*(factor_X1-1.0)/(pow(factor_X1,(200.0))-1.0);	//for increment in section 1
	for (int i=14;i<=212;i++)
		dx[i]=factor_X1*dx[i-1];

	dx[213]=40.0*(factor_X2-1.0)/(pow(factor_X2,(280.0))-1.0);	//for increment in section 2
	for (int i=214;i<=N;i++)
		dx[i]=factor_X2*dx[i-1];

	dy[2]=2.0*(factor_Y-1.0)/(pow(factor_Y,(M/2-1))-1.0);	//for increment
	for (int j=3;j<=M/2;j++)
		dy[j]=factor_Y*dy[j-1];

	factor_Y=1.0/factor_Y;	//for decrement

	dy[M/2+1]=2.0*(1.0-factor_Y)/(1.0-pow(factor_Y,(M/2)));

	for (int j=M/2+2;j<=M;j++)
		dy[j]=factor_Y*dy[j-1];

	for (int j=1;j<=M;j++)	//GENERATION OF SCALAR NODES
	{
		P[0][1][j]=0.0;	//left boundary
		P[1][1][j]=y;
		for (int i=2;i<=N;i++)
		{
			if (i==N)
				P[0][i][j]=LENGTH;	//right boundary
			else
				P[0][i][j]=P[0][i-1][j]+dx[i];

			P[1][i][j]=y;
		}
		if (j==M-1)
			y=HEIGHT;	//top boundary
		else
			y+=dy[j+1];
	}

	for (int j=1;j<=M;j++)	//GENERATION OF FORWARD STAGGERED X VELOCITY NODES
	{
		U[0][1][j]=0.0;	//left boundary
		U[1][1][j]=P[1][1][j];
		for (int i=2;i<N-1;i++)
		{
			U[0][i][j]=(P[0][i+1][j]+P[0][i][j])/2.0;	//PRACTICE A TYPE GRIDS
			U[1][i][j]=P[1][i][j];
		}
		U[0][N-1][j]=P[0][N][j];	//right boundary
		U[1][N-1][j]=P[1][N][j];
	}

	for (int j=1;j<M;j++)	//GENERATION OF FORWARD STAGGERED Y VELOCITY NODES
	{
		if (j==1)
			V[1][1][j]=P[1][1][j];	//bottom boundary
		else if (j==M-1)
			V[1][1][j]=P[1][1][M];	//top boundary
		else
			V[1][1][j]=(P[1][1][j]+P[1][1][j+1])/2.0;	//PRACTICE A TYPE GRIDS
		for (int i=1;i<=N;i++)
		{
			V[0][i][j]=P[0][i][1];
			V[1][i][j]=V[1][1][j];
		}
	}
	cout<<endl<<"GRID GENERATED"<<endl;		
}

void PVC :: U_bc_top()
{
	for (int i=1;i<N;i++)
		U[2][i][M]=0.0;	//DIRICHLET
}

void PVC :: U_bc_bottom()
{
	for (int i=1;i<N;i++)
		U[2][i][1]=0.0;	//DIRICHLET
}

void PVC :: U_bc_left()
{
	for (int j=1;j<=M;j++)
		U[2][1][j]=0.0;	//DIRICHLET
}

void PVC :: U_bc_right()
{
	for (int j=1;j<=M;j++)
		U[2][N-1][j]=U[2][N-2][j];	//NEUMANN

	float sum=0.0;
	for (int j=1;j<M;j++)
		sum+=(U[2][N-1][j+1]+U[2][N-1][j])*(U[1][N-1][j+1]-U[1][N-1][j]);
	U_out=0.5*sum/(U[1][N-1][M]-U[1][N-1][1]);	//average velocity at outlet
}

void PVC :: V_bc_top()
{
	for (int i=1;i<=N;i++)
	{
		if (V[0][i][M-1]<=0.5)
			V[2][i][M-1]=-1.0;	//DIRICHLET
		else
			V[2][i][M-1]=0.0;
	}
}

void PVC :: V_bc_bottom()
{
	for (int i=1;i<=N;i++)
		V[2][i][1]=0.0;	//DIRICHLET
}

void PVC :: V_bc_left()
{
	for (int j=1;j<M;j++)
		V[2][1][j]=V[2][2][j];	//NEUMANN
}

void PVC :: V_bc_right()
{
	for (int j=1;j<M;j++)
		V[2][N][j]=V[2][N-1][j];	//NEUMANN
}

void PVC :: P_bc_top()
{
	for (int i=1;i<=N;i++)
		P[2][i][M]=P[2][i][M-1];	//NEUMANN
}

void PVC :: P_bc_right()
{
	for (int j=1;j<=M;j++)
		P[2][N][j]=P[2][N-1][j];	//NEUMANN
}

void PVC :: P_bc_bottom()
{
	for (int i=1;i<=N;i++)
		P[2][i][1]=P[2][i][2];	//NEUMANN
}

void PVC :: P_bc_left()
{
	for (int j=1;j<=M;j++)
		P[2][1][j]=P[2][2][j];	//NEUMANN
}
#include "GI_schemes_SOU_QUICK_flow.cpp"	//Grid Independent schemes written in GI_schemes.cpp file

void PVC :: coeff_U_mom(int j, int i)
{
	float F_e=0.0, F_w=0.0, F_n=0.0, F_s=0.0;
	float D_e=0.0, D_w=0.0, D_n=0.0, D_s=0.0;	
	A_W=0.0; A_N=0.0; A_S=0.0; A_E=0.0; A_P=0.0; S=0.0; A_P0=0.0;	//re-initialization of the variables

	A_P0=(V[0][i+1][j]-V[0][i][j])*(V[1][i][j]-V[1][i][j-1])/dt;	//transient coeffecient

	F_e = 0.5*(V[1][i+1][j]-V[1][i+1][j-1])*(U[2][i+1][j]+U[2][i][j]);	//flux strengths	
	F_w = 0.5*(V[1][i][j]-V[1][i][j-1])*(U[2][i][j]+U[2][i-1][j]);
	F_n = 0.5*(V[0][i+1][j]-V[0][i][j])*(V[2][i+1][j]+V[2][i][j]);
	F_s = 0.5*(V[0][i+1][j-1]-V[0][i][j-1])*(V[2][i+1][j-1]+V[2][i][j-1]);
	
	D_e = (eta_f/r_f)*(V[1][i+1][j]-V[1][i+1][j-1])/((U[0][i+1][j]-U[0][i][j])*Re);	//diffusion coeffecients
	D_w = (eta_f/r_f)*(V[1][i][j]-V[1][i][j-1])/((U[0][i][j]-U[0][i-1][j])*Re);
	D_n = (eta_f/r_f)*(V[0][i+1][j]-V[0][i][j])/((U[1][i][j+1]-U[1][i][j])*Re);
	D_s = (eta_f/r_f)*(V[0][i+1][j-1]-V[0][i][j-1])/((U[1][i][j]-U[1][i][j-1])*Re);
	
	if ((j==M-1)||(i==2)||(j==2)||(i==N-2))
		CD_U(j,i,F_e,F_w,F_n,F_s,D_e,D_w,D_n,D_s);	//CENTRAL DIFFERENCE SCHEME FOR THE BORDER CV's
	else
		QUICK_U(j,i,F_e,F_w,F_n,F_s,D_e,D_w,D_n,D_s);	//QUICK SCHEME

	S+=(P[2][i][j]-P[2][i+1][j])*(V[1][i][j]-V[1][i][j-1])/r_f;	//addition of the pressure source term	
	A_P+=A_P0;	//addition of the psuedo-transient term
	
	U_A_P[i][j]=A_P;	//storing the A_P coeffecient for use in pressure correction
	U_A_N[i][j]=A_N;	//storing the other coeffecients for use in residual calculation
	U_A_E[i][j]=A_E;	
	U_A_W[i][j]=A_W;	
	U_A_S[i][j]=A_S;
	U_S[i][j]=S;
	U_A_P0[i][j]=A_P0;	
}

void PVC :: U_mom()
{
	for (int z=M-1;z>=2;z--)	//the jacobian matrix changes for every row z=row
	{	
		//1st row of the Jacobian matrix
		coeff_U_mom(z,2);
		a[0]=0.0;
		b[0]=A_P;
		c[0]=-A_E;		
		d[0]=S+A_W*U[2][1][z]+A_N*U[2][2][z+1]+A_S*U[2][2][z-1]+A_P0*U_prev[2][z];
	
		for (int i=1;i<N-4;i++)	//generation of the coeff mat i=column of b matrix
		{
			coeff_U_mom(z,i+2);
			c[i]=-A_E;
			b[i]=A_P;
			a[i]=-A_W;
			d[i]=S+A_N*U[2][i+2][z+1]+A_S*U[2][i+2][z-1]+A_P0*U_prev[i+2][z];	
		}
	
		//last row of the Jacobian matrix
		coeff_U_mom(z,N-2);
		a[N-4]=-A_W;
		b[N-4]=A_P;
		c[N-4]=0.0;
		d[N-4]=S+A_E*U[2][N-1][z]+A_N*U[2][N-2][z+1]+A_S*U[2][N-2][z-1]+A_P0*U_prev[N-2][z];
		
		tdma_bs(N-3,a,b,c,d);	//TDMA with Backward substitution present in "cfd_solvers.h" header file
		//Theta value for the particular row is found out
		for (int i=2;i<N-1;i++)
			U[2][i][z]=a[i-2];	//updation of the U velocity		
	}
}

void PVC :: U_mom_solver()
{
	for (int j=M;j>=1;j--)
		for (int i=1;i<=N-1;i++)
			U_prev[i][j]=U[2][i][j];	//storing the previous values of U velocity.
	do
	{
		U_mom();
		U_bc_right();	//neumann boundary update
	}
	while (U_mom_res() > 1e-5);
}

float PVC :: U_mom_res()
{
	double cc=0.0;	//residual values
	
	for (int j=M-1;j>=2;j--)
		for (int i=2;i<N-1;i++)
			cc+=pow((-U_A_P[i][j]*U[2][i][j]+U_A_E[i][j]*U[2][i+1][j]+U_A_W[i][j]*U[2][i-1][j]+U_A_N[i][j]*U[2][i][j+1]+U_A_S[i][j]*U[2][i][j-1]+U_A_P0[i][j]*U_prev[i][j]+U_S[i][j]),2.0);
	
	cc=sqrt(cc);
	//cout<<"U="<<cc<<endl;
	return (cc);
}

void PVC :: coeff_V_mom(int j, int i)
{
	float F_e=0.0, F_w=0.0, F_n=0.0, F_s=0.0;
	float D_e=0.0, D_w=0.0, D_n=0.0, D_s=0.0;	
	A_W=0.0; A_N=0.0; A_S=0.0; A_E=0.0; A_P=0.0; A_P0=0.0; S=0.0;	//re-initialization of the variables

	A_P0=(U[0][i][j]-U[0][i-1][j])*(U[1][i][j+1]-U[1][i][j])/dt;	//transient coeffecient

	F_e = 0.5*(U[2][i][j+1]+U[2][i][j])*(U[1][i][j+1]-U[1][i][j]);	//flux strengths
	F_w = 0.5*(U[2][i-1][j+1]+U[2][i-1][j])*(U[1][i-1][j+1]-U[1][i-1][j]);	
	F_n = 0.5*(V[2][i][j]+V[2][i][j+1])*(U[0][i][j+1]-U[0][i-1][j+1]);
	F_s = 0.5*(V[2][i][j]+V[2][i][j-1])*(U[0][i][j]-U[0][i-1][j]);
	
	D_e = (eta_f/r_f)*(U[1][i][j+1]-U[1][i][j])/(Re*(V[0][i+1][j]-V[0][i][j]));	//diffusion coeffecients
	D_w = (eta_f/r_f)*(U[1][i-1][j+1]-U[1][i-1][j])/(Re*(V[0][i][j]-V[0][i-1][j]));
	D_n = (eta_f/r_f)*(U[0][i][j+1]-U[0][i-1][j+1])/(Re*(V[1][i][j+1]-V[1][i][j]));
	D_s = (eta_f/r_f)*(U[0][i][j]-U[0][i-1][j])/(Re*(V[1][i][j]-V[1][i][j-1]));	
		
	if ((j==2)||(i==2)||(j==M-2)||(i==N-1))
		CD_V(j,i,F_e,F_w,F_n,F_s,D_e,D_w,D_n,D_s);	//CENTRAL DIFFERENCE SCHEME FOR THE BORDER CV's
	else
		QUICK_V(j,i,F_e,F_w,F_n,F_s,D_e,D_w,D_n,D_s);	//QUICK SCHEME

	S+=(P[2][i][j]-P[2][i][j+1])*(U[0][i][j]-U[0][i-1][j])/r_f;	//addition of the pressure source term	
	A_P+=A_P0;	//addition of psuedo-transient term

	V_A_P[i][j]=A_P;	//storing the A_P coeffecient for use in pressure correction
	V_A_N[i][j]=A_N;	//storing the other coeffecients for use in residual calculation
	V_A_E[i][j]=A_E;	
	V_A_W[i][j]=A_W;	
	V_A_S[i][j]=A_S;
	V_S[i][j]=S;
	V_A_P0[i][j]=A_P0;	
}

void PVC :: V_mom()	
{
	for (int z=M-2;z>=2;z--)	//the jacobian matrix changes for every row z=row 
	{	
		//1st row of the Jacobian matrix
		coeff_V_mom(z,2);	
		a[0]=0.0;
		b[0]=A_P;
		c[0]=-A_E;		
		d[0]=S+A_W*V[2][1][z]+A_N*V[2][2][z+1]+A_S*V[2][2][z-1]+A_P0*V_prev[2][z];
	
		for (int i=1;i<N-3;i++)	//generation of the coeff mat i=column of b matrix
		{
			coeff_V_mom(z,i+2);	
			c[i]=-A_E;
			b[i]=A_P;
			a[i]=-A_W;
			d[i]=S+A_N*V[2][i+2][z+1]+A_S*V[2][i+2][z-1]+A_P0*V_prev[i+2][z];	
		}
	
		//last row of the Jacobian matrix
		coeff_V_mom(z,N-1);
		a[N-3]=-A_W;
		b[N-3]=A_P;		
		c[N-3]=0.0;
		d[N-3]=S+A_E*V[2][N][z]+A_N*V[2][N-1][z+1]+A_S*V[2][N-1][z-1]+A_P0*V_prev[N-1][z];
	
		tdma_bs(N-2,a,b,c,d);	//TDMA with Backward substitution present in "cfd_solvers.h" header file
		//Theta value for the particular row is found out
		for (int i=2;i<=N-1;i++)
			V[2][i][z]=a[i-2];	//updation of the V velocity		
	}
}

void PVC :: V_mom_solver()
{
	for (int j=M-1;j>=1;j--)
		for (int i=1;i<=N;i++)
			V_prev[i][j]=V[2][i][j];	//storing the previous values of V velocity.
	do
	{
		V_mom();
		V_bc_left();
		V_bc_right();	//neumann boundary update
	}
	while (V_mom_res() > 1e-5);
}

float PVC :: V_mom_res()
{
	double cc=0.0;	//residual values
	
	for (int j=M-2;j>1;j--)
		for (int i=2;i<=N-1;i++)
			cc+=pow((-V_A_P[i][j]*V[2][i][j]+V_A_E[i][j]*V[2][i+1][j]+V_A_W[i][j]*V[2][i-1][j]+V_A_N[i][j]*V[2][i][j+1]+V_A_S[i][j]*V[2][i][j-1]+V_A_P0[i][j]*V_prev[i][j]+V_S[i][j]),2.0);

	cc=sqrt(cc);
	//cout<<"V="<<cc<<endl;
	return (cc);
}

void PVC :: coeff_P_dash(int j, int i)
{	
	A_W=0.0; A_N=0.0; A_S=0.0; A_E=0.0; A_P=0.0; S=0.0;	//re-initialization of the variables
	
	if (i!=2)
		A_W = pow((V[1][i][j]-V[1][i][j-1]),2)/U_A_P[i-1][j];		
	if (i!=N-1)
		A_E = pow((V[1][i][j]-V[1][i][j-1]),2)/U_A_P[i][j];		
	if (j!=2)
		A_S = pow((U[0][i][j]-U[0][i-1][j]),2)/V_A_P[i][j-1];		
	if (j!=M-1)
		A_N = pow((U[0][i][j]-U[0][i-1][j]),2)/V_A_P[i][j];
		
	A_P = A_E + A_W + A_N + A_S;
	S=(U[2][i][j]-U[2][i-1][j])*(V[1][i][j]-V[1][i][j-1])+(V[2][i][j]-V[2][i][j-1])*(U[0][i][j]-U[0][i-1][j]);	//continuity source term
}

void PVC :: P_dash()
{	
	for (int z=M-1;z>=2;z--)	//the jacobian matrix changes for every row z=row
	{	
		//1st row of the Jacobian matrix
		coeff_P_dash(z,2);	
		a[0]=0.0;
		b[0]=A_P-A_S*(theta-1.0);
		c[0]=-A_E;
		d[0]=A_W*Pd[1][z]+A_N*Pd[2][z+1]+A_S*(Pd[2][z-1]-(theta-1.0)*Pd[2][z])-S;
	
		for (int i=1;i<N-3;i++)	//generation of the coeff mat i=column
		{
			coeff_P_dash(z,i+2);	
			c[i]=-A_E;
			b[i]=(A_P-A_S*(theta-1.0));
			a[i]=-A_W;
			d[i]=A_N*Pd[i+2][z+1]+A_S*(Pd[i+2][z-1]-(theta-1.0)*Pd[i+2][z])-S;
		}
	
		//last row of the Jacobian matrix
		coeff_P_dash(z,N-1);
		a[N-3]=-A_W;
		b[N-3]=A_P-A_S*(theta-1.0);		
		c[N-3]=0.0;
		d[N-3]=A_E*Pd[N][z]+A_N*Pd[N-1][z+1]+A_S*(Pd[N-1][z-1]-(theta-1.0)*Pd[N-1][z])-S;
	
		tdma_bs(N-2,a,b,c,d);	//TDMA with Backward substitution present in "cfd_solvers.h" header file
		//Theta value for the particular row is found out
		for (int i=2;i<=N-1;i++)
			Pd[i][z]=a[i-2];	//updation of the P-dash
	}		
}

void PVC :: P_dash_solver()
{
	for (int j=1;j<=M;j++)
		for (int i=1;i<=N;i++)
			Pd[i][j]=0.0;		//P_dash is taken as 0 for every iteration

	int count=0;
	
	float cc_ini=P_dash_res(1);
	do
	{
		P_dash();
		count++;

		/*for (int i=1;i<=N;i++)	//neumann boundary condition for P-dash
		{
			Pd[i][M]=Pd[i][M-1];	//top boundary condition
			Pd[i][1]=Pd[i][2];	//bottom boundary condition
		}
		for (int j=1;j<=M;j++)
		{
			Pd[1][j]=Pd[2][j];	//left boundary condition
			Pd[N][j]=Pd[N-1][j];	//right boundary condition
		}*/		
	}
	while ((P_dash_res(count) > Gp*cc_ini)&&(count<=750));
}

float PVC :: P_dash_res(int count)
{
	double cc=0.0;	//residual values
	
	for (int j=M-1;j>1;j--)
	{
		for (int i=2;i<=N-1;i++)
		{
			coeff_P_dash(j,i);
			cc+=pow((A_E*Pd[i+1][j]+A_W*Pd[i-1][j]+A_N*Pd[i][j+1]+A_S*Pd[i][j-1]-S-A_P*Pd[i][j]),2.0);
		}
	}
	cc=sqrt(cc);
	//if (count%50==0)
	//	cout<<"P="<<cc<<endl;
	return (cc);
}

void PVC :: P_corr()
{
	for (int j=2;j<=M-1;j++)	
		for (int i=2;i<=N-1;i++)
			P[2][i][j]+=rex_P*Pd[i][j];		
}

void PVC :: U_corr()
{	
	for (int j=2;j<=M-1;j++)
		for (int i=2;i<N-1;i++)
			U[2][i][j]+=(1.0/U_A_P[i][j])*(Pd[i][j]-Pd[i+1][j])*(V[1][i][j]-V[1][i][j-1]);
}

void PVC :: V_corr()
{
	for (int j=2;j<M-1;j++)
		for (int i=2;i<=N-1;i++)			
			V[2][i][j]+=(1.0/V_A_P[i][j])*(Pd[i][j]-Pd[i][j+1])*(U[0][i][j]-U[0][i-1][j]);
}

void PVC :: SIMPLE()
{
	FILE *p;	//extracting mass residuals
	/*p=fopen("mass_res.dat","w");	//generate cont_res.dat file
	fprintf (p,"TITLE = \"MASS RESIDUAL PLOT for Re = %f\"\n",Re);	
	fprintf (p,"VARIABLES = \"Iteration\", \"Residual*10^7\"\n");*/
	do
	{
		cnt++;		
		if ((cnt%30)==0)
			cout<<endl<<"ITERATION "<<cnt<<endl;

		if ((cnt%100)==0)
		{
			write();
			cout<<endl<<"INTERMEDIATE FILE OUTPUT SUCCESSFUL"<<endl;
		}

		U_mom_solver(); V_mom_solver();
		P_dash_solver();
		P_corr(); U_corr(); V_corr();

		U_bc_right();	//neumann boundary update
		V_bc_left();
		V_bc_right();
	}
	while (cont_source(p)>=18e-7);
	
	P_bc_bottom();	//updation of the pressure boundary condition
	P_bc_left();
	P_bc_right();
	P_bc_top();
	
	//fclose(p);
	//cout<<endl<<"CONTINUITY RESIDUAL FILE OUTPUT SUCCESSFUL"<<endl;
}

float PVC :: cont_source(FILE *p)	//REQUIRED TO BE CHANGED !!!!
{
	double cont_s=0.0;
	
	for (int j=2;j<=M-1;j++)
		for (int i=2;i<=N-1;i++)
			cont_s+=fabs((U[2][i][j]-U[2][i-1][j])*(V[1][i][j]-V[1][i][j-1])+(V[2][i][j]-V[2][i][j-1])*(U[0][i][j]-U[0][i-1][j]));
	
	//fprintf(p,"%d\t%lf\n",cnt,cont_s*1e7);
		
	//if ((cnt%30)==0)
		cout<<"CONTINUITY="<<cont_s<<endl;
	
	return cont_s;
}

void PVC :: write()
{
	FILE *p;
	p=fopen("P.dat","w");

	fprintf (p,"TITLE = \"PRESSURE\"\n");	
	fprintf (p,"VARIABLES = \"X\", \"Y\", \"P\"\n");	
	fprintf (p,"ZONE I=%d,J=%d,DATAPACKING=BLOCK\n",N,M);

	for (int HH=0;HH<=2;HH++)
	{
		for (int j=1;j<=M;j++)
		{
			for (int i=1;i<=N;i++)
				fprintf (p," %.4f",P[HH][i][j]);
			fprintf (p,"\n");
		}
		fprintf (p,"\n\n");
	}
	fclose(p);

	p=fopen("U.dat","w");

	fprintf (p,"TITLE = \"X VELOCITY\"\n");	
	fprintf (p,"VARIABLES = \"X\", \"Y\", \"U\"\n");	
	fprintf (p,"ZONE I=%d,J=%d,DATAPACKING=BLOCK\n",N-1,M);

	for (int HH=0;HH<=2;HH++)
	{
		for (int j=1;j<=M;j++)
		{
			for (int i=1;i<N;i++)
				fprintf (p," %.4f",U[HH][i][j]);
			fprintf (p,"\n");
		}
		fprintf (p,"\n\n");
	}
	fclose(p);

	p=fopen("V.dat","w");

	fprintf (p,"TITLE = \"Y VELOCITY\"\n");	
	fprintf (p,"VARIABLES = \"X\", \"Y\", \"V\"\n");	
	fprintf (p,"ZONE I=%d,J=%d,DATAPACKING=BLOCK\n",N,M-1);

	for (int HH=0;HH<=2;HH++)
	{
		for (int j=1;j<M;j++)
		{
			for (int i=1;i<=N;i++)
				fprintf (p," %.4f",V[HH][i][j]);
			fprintf (p,"\n");
		}
		fprintf (p,"\n\n");
	}
	fclose(p);
	
	cout<<endl<<"FILE OUTPUT SUCCESSFUL"<<endl;
}

int main()
{
	clock_t start = clock();
	
	PVC ms;

	ms.getdata (100.0,0.001, 1.057739,1.162896);
	//getdata(float re, float rexp, float r, float eta);		//initialisation of the variables
	ms.SIMPLE();
	ms.write();
		
	cout<<endl<<"END OF SIMULATION"<<endl;
	
	clock_t end = clock();
	printf("\nElapsed Time : %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
	return(0);
}
