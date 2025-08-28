//THIS CODE IS CREATED TO READ THE INTERMEDIATE FILES AND CO-LOCATE THE VARIABLES FOR STREAMLINE PLOT IN TECPLOT
//THIS IS A STANDALONE CODE.

#include<iostream>
#include<cstdio>
#include<cstdlib>

using namespace std;

const int M=60,N=492;	//M->No of grids in vertical direction, N->No of grids in horizontal direction
const float Re = 100;	//Reynolds number

class cntr
{
	float ***P, ***U, ***V;

	public:
		cntr();
		~cntr();
		void fread();
		void test_write();
		void write();
};

cntr :: cntr()
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
}

cntr :: ~cntr()
{
	free(U); free(V); free(P);
}

void cntr :: fread()
{
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

	for (int HH=0;HH<=2;HH++)
		for (int j=1;j<=M;j++)
			for (int i=1;i<=N;i++)
				fscanf (p,"%f",&P[HH][i][j]);
	fclose(p);
}

void cntr :: test_write()
{
	FILE *p;
	p=fopen("U_test.dat","w");

	fprintf (p,"TITLE = \"TEST OF U\"\n");	
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

	p=fopen("V_test.dat","w");

	fprintf (p,"TITLE = \"TEST OF V\"\n");	
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

	p=fopen("P_test.dat","w");

	fprintf (p,"TITLE = \"TEST OF P\"\n");	
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

	cout<<"\nTEST FILE OUTPUT SUCCESSFULL"<<endl;
}

void cntr :: write()
{
	FILE *p; float temp;
	p=fopen("vel_press.dat","w");	//generate vel_press.dat file
	fprintf (p,"TITLE = \"CO-LOCATED PLOTS FOR Re = %f\"\n",Re);	
	fprintf (p,"VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"P\"\n");	
	fprintf (p,"ZONE I=%d,J=%d,DATAPACKING=BLOCK\n",N,M);

	for (int j=1;j<=M;j++)	//print X co-ordinates
	{
		for (int i=1;i<=N;i++)
			fprintf (p," %.4f",P[0][i][j]);
			
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");
	
	for (int j=1;j<=M;j++)	//print Y co-ordinates
	{
		for (int i=1;i<=N;i++)
			fprintf (p," %.4f",P[1][i][j]);
			
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");

	for (int j=1;j<=M;j++)	//print X velocity
	{
		for (int i=1;i<=N;i++)
		{
			if (i==1)	//left boundary
				fprintf (p," %.4f",U[2][i][j]);
			else if (i==N)	//right boundary
				fprintf (p," %.4f",U[2][i-1][j]);
			else	//inner domain
			{
				temp=((U[0][i][j]-P[0][i][j])*U[2][i-1][j]+(P[0][i][j]-U[0][i-1][j])*U[2][i][j])/(U[0][i][j]-U[0][i-1][j]);	//lever rule
				fprintf (p," %.4f",temp);
			}
		}
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");

	for (int j=1;j<=M;j++)	//print Y velocity
	{
		for (int i=1;i<=N;i++)
		{
			if (j==1)	//bottom boundary
				fprintf (p," %.4f",V[2][i][j]);
			else if (j==M)	//top boundary
				fprintf (p," %.4f",V[2][i][j-1]);
			else	//inner domain
			{
				temp=((V[1][i][j]-P[1][i][j])*V[2][i][j-1]+(P[1][i][j]-V[1][i][j-1])*V[2][i][j])/(V[1][i][j]-V[1][i][j-1]);	//lever rule
				fprintf (p," %.4f",temp);
			}
		}
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");

	for (int j=1;j<=M;j++)	//print pressure
	{
		for (int i=1;i<=N;i++)
			fprintf (p," %.4f",P[2][i][j]);
			
		fprintf (p,"\n");
	}
	fprintf (p,"\n\n");

	fclose(p);
	cout<<endl<<"FILE OUTPUT SUCCESSFUL"<<endl;
}

int main()
{
	cntr ms;

	ms.fread();
	ms.write();

	return (0);
}
