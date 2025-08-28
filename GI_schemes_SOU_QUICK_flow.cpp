//CENTRAL DIFFERENCE SCHEME AND COMBINED SOU WITH DEFERRED QUICK SCHEME FOR A FORWARD STAGGERED NON UNIFORM GRID.

void PVC :: CD_geo_U(int j, int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float F_e, float F_w, float F_n, float F_s)
{
	float D_1w=0.0,D_2w=0.0, D_1e=0.0,D_2e=0.0, D_1n=0.0,D_2n=0.0, D_1s=0.0,D_2s=0.0;

	D_1w = (V[0][i][j]-U[0][i-1][j])*CHK(F_w)+(U[0][i][j]-V[0][i][j])*CHK(-F_w);

	D_2w = (U[0][i][j]-V[0][i][j])*CHK(F_w)+(V[0][i][j]-U[0][i-1][j])*CHK(-F_w);

	*A_1w = D_1w/(D_1w+D_2w); *A_2w = D_2w/(D_1w+D_2w);
//------------------------------------------------------------------------------------------------------------------------------
	D_1e = (V[0][i+1][j]-U[0][i][j])*CHK(F_e)+(U[0][i+1][j]-V[0][i+1][j])*CHK(-F_e);

	D_2e = (U[0][i+1][j]-V[0][i+1][j])*CHK(F_e)+(V[0][i+1][j]-U[0][i][j])*CHK(-F_e);

	*A_1e = D_1e/(D_1e+D_2e); *A_2e = D_2e/(D_1e+D_2e);
//---------------------------------------------------------------------------------------------------------------------------
	D_1n = (V[1][i][j]-U[1][i][j])*CHK(F_n)+(U[1][i][j+1]-V[1][i][j])*CHK(-F_n);

	D_2n = (U[1][i][j+1]-V[1][i][j])*CHK(F_n)+(V[1][i][j]-U[1][i][j])*CHK(-F_n);

	*A_1n = D_1n/(D_1n+D_2n); *A_2n = D_2n/(D_1n+D_2n);
//---------------------------------------------------------------------------------------------------------------------------
	D_1s = (V[1][i][j-1]-U[1][i][j-1])*CHK(F_s)+(U[1][i][j]-V[1][i][j-1])*CHK(-F_s);

	D_2s = (U[1][i][j]-V[1][i][j-1])*CHK(F_s)+(V[1][i][j-1]-U[1][i][j-1])*CHK(-F_s);

	*A_1s = D_1s/(D_1s+D_2s); *A_2s = D_2s/(D_1s+D_2s);
}

void PVC :: CD_geo_V(int j, int i,float *A_1w,float *A_2w,float *A_1e,float *A_2e,float *A_1n,float *A_2n,float *A_1s,float *A_2s,float F_e, float F_w, float F_n, float F_s)
{
	float D_1w=0.0,D_2w=0.0, D_1e=0.0,D_2e=0.0, D_1n=0.0,D_2n=0.0, D_1s=0.0,D_2s=0.0;

	D_1w = (U[0][i-1][j]-V[0][i-1][j])*CHK(F_w)+(V[0][i][j]-U[0][i-1][j])*CHK(-F_w);

	D_2w = (V[0][i][j]-U[0][i-1][j])*CHK(F_w)+(U[0][i-1][j]-V[0][i-1][j])*CHK(-F_w);

	*A_1w = D_1w/(D_1w+D_2w); *A_2w = D_2w/(D_1w+D_2w);
//---------------------------------------------------------------------------------------------------------------------------
	D_1e = (U[0][i][j]-V[0][i][j])*CHK(F_e)+(V[0][i+1][j]-U[0][i][j])*CHK(-F_e);

	D_2e = (V[0][i+1][j]-U[0][i][j])*CHK(F_e)+(U[0][i][j]-V[0][i][j])*CHK(-F_e);

	*A_1e = D_1e/(D_1e+D_2e); *A_2e = D_2e/(D_1e+D_2e);
//---------------------------------------------------------------------------------------------------------------------------
	D_1n = (U[1][i][j+1]-V[1][i][j])*CHK(F_n)+(V[1][i][j+1]-U[1][i][j+1])*CHK(-F_n);

	D_2n = (V[1][i][j+1]-U[1][i][j+1])*CHK(F_n)+(U[1][i][j+1]-V[1][i][j])*CHK(-F_n);

	*A_1n = D_1n/(D_1n+D_2n); *A_2n = D_2n/(D_1n+D_2n);
//---------------------------------------------------------------------------------------------------------------------------
	D_1s = (U[1][i][j]-V[1][i][j-1])*CHK(F_s)+(V[1][i][j]-U[1][i][j])*CHK(-F_s);

	D_2s = (V[1][i][j]-U[1][i][j])*CHK(F_s)+(U[1][i][j]-V[1][i][j-1])*CHK(-F_s);

	*A_1s = D_1s/(D_1s+D_2s); *A_2s = D_2s/(D_1s+D_2s);
}

void PVC :: CD_U(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s)
{
	float A_1w=0.0,A_2w=0.0, A_1e=0.0,A_2e=0.0, A_1n=0.0,A_2n=0.0, A_1s=0.0,A_2s=0.0;
	float B_E=0.0, B_W=0.0, B_N=0.0, B_S=0.0, B_P=0.0;

	CD_geo_U(j,i,&A_1w,&A_2w,&A_1e,&A_2e,&A_1n,&A_2n,&A_1s,&A_2s,F_e,F_w,F_n,F_s);

	A_E = D_e+MAX(0.0,-F_e);

	A_W = D_w+MAX(F_w,0.0);

	A_N = D_n+MAX(0.0,-F_n);

	A_S = D_s+MAX(F_s,0.0);

	A_P = A_E+A_W+A_N+A_S+(F_e-F_w)+(F_n-F_s);

	B_E = -A_1e*(MAX(F_e,0.0)+MAX(0.0,-F_e));

	B_W = -A_1w*(MAX(F_w,0.0)+MAX(0.0,-F_w));

	B_N = -A_1n*(MAX(F_n,0.0)+MAX(0.0,-F_n));

	B_S = -A_1s*(MAX(F_s,0.0)+MAX(0.0,-F_s));

	B_P = -(B_E+B_W+B_N+B_S);

	S = B_E*U[2][i+1][j]+B_W*U[2][i-1][j]+B_N*U[2][i][j+1]+B_S*U[2][i][j-1]+B_P*U[2][i][j];
}

void PVC :: CD_V(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s)
{
	float A_1w=0.0,A_2w=0.0, A_1e=0.0,A_2e=0.0, A_1n=0.0,A_2n=0.0, A_1s=0.0,A_2s=0.0;
	float B_E=0.0, B_W=0.0, B_N=0.0, B_S=0.0, B_P=0.0;

	CD_geo_V(j,i,&A_1w,&A_2w,&A_1e,&A_2e,&A_1n,&A_2n,&A_1s,&A_2s,F_e,F_w,F_n,F_s);

	A_E = D_e+MAX(0.0,-F_e);

	A_W = D_w+MAX(F_w,0.0);

	A_N = D_n+MAX(0.0,-F_n);

	A_S = D_s+MAX(F_s,0.0);

	A_P = A_E+A_W+A_N+A_S+(F_e-F_w)+(F_n-F_s);

	B_E = -A_1e*(MAX(F_e,0.0)+MAX(0.0,-F_e));

	B_W = -A_1w*(MAX(F_w,0.0)+MAX(0.0,-F_w));

	B_N = -A_1n*(MAX(F_n,0.0)+MAX(0.0,-F_n));

	B_S = -A_1s*(MAX(F_s,0.0)+MAX(0.0,-F_s));

	B_P = -(B_E+B_W+B_N+B_S);

	S = B_E*V[2][i+1][j]+B_W*V[2][i-1][j]+B_N*V[2][i][j+1]+B_S*V[2][i][j-1]+B_P*V[2][i][j];
}

void PVC :: QUICK_geo_U(int j,int i,float *G_1w,float *G_3w,float *G_1e,float *G_3e,float *G_1n,float *G_3n,float *G_1s,float *G_3s,float *K_1w,float *K_3w,float *K_1e,float *K_3e,float *K_1n,float *K_3n,float *K_1s,float *K_3s,float *r_w,float *r_e,float *r_n,float *r_s,float F_e, float F_w, float F_n, float F_s)
{
	float D_1w=0.0,D_2w=0.0,D_3w=0.0, D_1e=0.0,D_2e=0.0,D_3e=0.0, D_1n=0.0,D_2n=0.0,D_3n=0.0, D_1s=0.0,D_2s=0.0,D_3s=0.0;

	D_1w = (V[0][i][j]-U[0][i-1][j])*CHK(F_w)+(U[0][i][j]-V[0][i][j])*CHK(-F_w);

	D_2w = (U[0][i][j]-V[0][i][j])*CHK(F_w)+(V[0][i][j]-U[0][i-1][j])*CHK(-F_w);

	D_3w = (V[0][i][j]-U[0][i-2][j])*CHK(F_w)+(U[0][i+1][j]-V[0][i][j])*CHK(-F_w);

	*G_1w = D_1w/(D_3w-D_1w); *G_3w = D_3w/(D_3w-D_1w);

	*K_1w = (D_1w+D_2w)/(D_3w-D_1w); *K_3w = (D_2w+D_3w)/(D_3w-D_1w);

	*r_w = D_1w*D_3w/((D_1w+D_2w)*(D_2w+D_3w));
//---------------------------------------------------------------------------------------------------------------------------
	D_1e = (V[0][i+1][j]-U[0][i][j])*CHK(F_e)+(U[0][i+1][j]-V[0][i+1][j])*CHK(-F_e);

	D_2e = (U[0][i+1][j]-V[0][i+1][j])*CHK(F_e)+(V[0][i+1][j]-U[0][i][j])*CHK(-F_e);

	D_3e = (V[0][i+1][j]-U[0][i-1][j])*CHK(F_e)+(U[0][i+2][j]-V[0][i+1][j])*CHK(-F_e);
	
	*G_1e = D_1e/(D_3e-D_1e); *G_3e = D_3e/(D_3e-D_1e);

	*K_1e = (D_1e+D_2e)/(D_3e-D_1e); *K_3e = (D_2e+D_3e)/(D_3e-D_1e);

	*r_e = D_1e*D_3e/((D_1e+D_2e)*(D_2e+D_3e));
//---------------------------------------------------------------------------------------------------------------------------
	D_1n = (V[1][i][j]-U[1][i][j])*CHK(F_n)+(U[1][i][j+1]-V[1][i][j])*CHK(-F_n);

	D_2n = (U[1][i][j+1]-V[1][i][j])*CHK(F_n)+(V[1][i][j]-U[1][i][j])*CHK(-F_n);

	D_3n = (V[1][i][j]-U[1][i][j-1])*CHK(F_n)+(U[1][i][j+2]-V[1][i][j])*CHK(-F_n);
	
	*G_1n = D_1n/(D_3n-D_1n); *G_3n = D_3n/(D_3n-D_1n);

	*K_1n = (D_1n+D_2n)/(D_3n-D_1n); *K_3n = (D_2n+D_3n)/(D_3n-D_1n);

	*r_n = D_1n*D_3n/((D_1n+D_2n)*(D_2n+D_3n));
//---------------------------------------------------------------------------------------------------------------------------
	D_1s = (V[1][i][j-1]-U[1][i][j-1])*CHK(F_s)+(U[1][i][j]-V[1][i][j-1])*CHK(-F_s);

	D_2s = (U[1][i][j]-V[1][i][j-1])*CHK(F_s)+(V[1][i][j-1]-U[1][i][j-1])*CHK(-F_s);

	D_3s = (V[1][i][j-1]-U[1][i][j-2])*CHK(F_s)+(U[1][i][j+1]-V[1][i][j-1])*CHK(-F_s);
	
	*G_1s = D_1s/(D_3s-D_1s); *G_3s = D_3s/(D_3s-D_1s);

	*K_1s = (D_1s+D_2s)/(D_3s-D_1s); *K_3s = (D_2s+D_3s)/(D_3s-D_1s);

	*r_s = D_1s*D_3s/((D_1s+D_2s)*(D_2s+D_3s));
}

void PVC :: QUICK_geo_V(int j,int i,float *G_1w,float *G_3w,float *G_1e,float *G_3e,float *G_1n,float *G_3n,float *G_1s,float *G_3s,float *K_1w,float *K_3w,float *K_1e,float *K_3e,float *K_1n,float *K_3n,float *K_1s,float *K_3s,float *r_w,float *r_e,float *r_n,float *r_s,float F_e, float F_w, float F_n, float F_s)
{
	float D_1w=0.0,D_2w=0.0,D_3w=0.0, D_1e=0.0,D_2e=0.0,D_3e=0.0, D_1n=0.0,D_2n=0.0,D_3n=0.0, D_1s=0.0,D_2s=0.0,D_3s=0.0;

	D_1w = (U[0][i-1][j]-V[0][i-1][j])*CHK(F_w)+(V[0][i][j]-U[0][i-1][j])*CHK(-F_w);

	D_2w = (V[0][i][j]-U[0][i-1][j])*CHK(F_w)+(U[0][i-1][j]-V[0][i-1][j])*CHK(-F_w);

	D_3w = (U[0][i-1][j]-V[0][i-2][j])*CHK(F_w)+(V[0][i+1][j]-U[0][i-1][j])*CHK(-F_w);

	*G_1w = D_1w/(D_3w-D_1w); *G_3w = D_3w/(D_3w-D_1w);

	*K_1w = (D_1w+D_2w)/(D_3w-D_1w); *K_3w = (D_2w+D_3w)/(D_3w-D_1w);

	*r_w = D_1w*D_3w/((D_1w+D_2w)*(D_2w+D_3w));
//---------------------------------------------------------------------------------------------------------------------------
	D_1e = (U[0][i][j]-V[0][i][j])*CHK(F_e)+(V[0][i+1][j]-U[0][i][j])*CHK(-F_e);

	D_2e = (V[0][i+1][j]-U[0][i][j])*CHK(F_e)+(U[0][i][j]-V[0][i][j])*CHK(-F_e);

	D_3e = (U[0][i][j]-V[0][i-1][j])*CHK(F_e)+(V[0][i+2][j]-U[0][i][j])*CHK(-F_e);
	
	*G_1e = D_1e/(D_3e-D_1e); *G_3e = D_3e/(D_3e-D_1e);

	*K_1e = (D_1e+D_2e)/(D_3e-D_1e); *K_3e = (D_2e+D_3e)/(D_3e-D_1e);

	*r_e = D_1e*D_3e/((D_1e+D_2e)*(D_2e+D_3e));
//---------------------------------------------------------------------------------------------------------------------------
	D_1n = (U[1][i][j+1]-V[1][i][j])*CHK(F_n)+(V[1][i][j+1]-U[1][i][j+1])*CHK(-F_n);

	D_2n = (V[1][i][j+1]-U[1][i][j+1])*CHK(F_n)+(U[1][i][j+1]-V[1][i][j])*CHK(-F_n);

	D_3n = (U[1][i][j+1]-V[1][i][j-1])*CHK(F_n)+(V[1][i][j+2]-U[1][i][j+1])*CHK(-F_n);
	
	*G_1n = D_1n/(D_3n-D_1n); *G_3n = D_3n/(D_3n-D_1n);

	*K_1n = (D_1n+D_2n)/(D_3n-D_1n); *K_3n = (D_2n+D_3n)/(D_3n-D_1n);

	*r_n = D_1n*D_3n/((D_1n+D_2n)*(D_2n+D_3n));
//---------------------------------------------------------------------------------------------------------------------------
	D_1s = (U[1][i][j]-V[1][i][j-1])*CHK(F_s)+(V[1][i][j]-U[1][i][j])*CHK(-F_s);

	D_2s = (V[1][i][j]-U[1][i][j])*CHK(F_s)+(U[1][i][j]-V[1][i][j-1])*CHK(-F_s);

	D_3s = (U[1][i][j]-V[1][i][j-2])*CHK(F_s)+(V[1][i][j+1]-U[1][i][j])*CHK(-F_s);
	
	*G_1s = D_1s/(D_3s-D_1s); *G_3s = D_3s/(D_3s-D_1s);

	*K_1s = (D_1s+D_2s)/(D_3s-D_1s); *K_3s = (D_2s+D_3s)/(D_3s-D_1s);

	*r_s = D_1s*D_3s/((D_1s+D_2s)*(D_2s+D_3s));
}

void PVC :: QUICK_U(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s)
{
	float G_1w=0.0,G_3w=0.0, G_1e=0.0,G_3e=0.0, G_1n=0.0,G_3n=0.0, G_1s=0.0,G_3s=0.0;
	float K_1w=0.0,K_3w=0.0, K_1e=0.0,K_3e=0.0, K_1n=0.0,K_3n=0.0, K_1s=0.0,K_3s=0.0;
	float r_w=0.0,r_e=0.0,r_n=0.0,r_s=0.0;
	float B_EE=0.0,B_WW=0.0,B_NN=0.0,B_SS=0.0,B_P=0.0;

	QUICK_geo_U(j,i,&G_1w,&G_3w,&G_1e,&G_3e,&G_1n,&G_3n,&G_1s,&G_3s,&K_1w,&K_3w,&K_1e,&K_3e,&K_1n,&K_3n,&K_1s,&K_3s,&r_w,&r_e,&r_n,&r_s,F_e,F_w,F_n,F_s);

	A_E = D_e+G_3e*MAX(-F_e,0.0)+G_1w*MAX(-F_w,0.0);

	A_W = D_w+G_1e*MAX(F_e,0.0)+G_3w*MAX(F_w,0.0);

	A_N = D_n+G_3n*MAX(-F_n,0.0)+G_1s*MAX(-F_s,0.0);

	A_S = D_s+G_1n*MAX(F_n,0.0)+G_3s*MAX(F_s,0.0);

	A_P = A_E+A_W+A_N+A_S+F_e*G_3e-F_w*G_3w-F_s*G_3s+F_n*G_3n-G_1e*MAX(F_e,0.0)-G_1n*MAX(F_n,0.0)-G_1w*MAX(0.0,-F_w)-G_1s*MAX(0.0,-F_s);

	S = -r_e*(U_prev[i+1][j]-K_3e*U_prev[i][j]+K_1e*U_prev[i-1][j])*MAX(F_e,0.0)+r_w*(U_prev[i][j]-K_3w*U_prev[i-1][j]
		+K_1w*U_prev[i-2][j])*MAX(F_w,0.0)-r_n*(U_prev[i][j+1]-K_3n*U_prev[i][j]+K_1n*U_prev[i][j-1])*MAX(F_n,0.0)
		+r_s*(U_prev[i][j]-K_3s*U_prev[i][j-1]+K_1s*U_prev[i][j-2])*MAX(F_s,0.0)+r_e*(U_prev[i][j]-K_3e*U_prev[i+1][j]
		+K_1e*U_prev[i+2][j])*MAX(-F_e,0.0)-r_w*(U_prev[i-1][j]-K_3w*U_prev[i][j]+K_1w*U_prev[i+1][j])*MAX(-F_w,0.0)
		+r_n*(U_prev[i][j]-K_3n*U_prev[i][j+1]+K_1n*U_prev[i][j+2])*MAX(-F_n,0.0)-r_s*(U_prev[i][j-1]-K_3s*U_prev[i][j]
		+K_1s*U_prev[i][j+1])*MAX(-F_s,0.0)-G_1w*U[2][i-2][j]*MAX(F_w,0.0)-G_1s*U[2][i][j-2]*MAX(F_s,0.0)-G_1e*U[2][i+2][j]*MAX(-F_e,0.0)
		-G_1n*U[2][i][j+2]*MAX(-F_n,0.0);
}

void PVC :: QUICK_V(int j, int i, float F_e, float F_w, float F_n, float F_s, float D_e, float D_w, float D_n, float D_s)
{
	float G_1w=0.0,G_3w=0.0, G_1e=0.0,G_3e=0.0, G_1n=0.0,G_3n=0.0, G_1s=0.0,G_3s=0.0;
	float K_1w=0.0,K_3w=0.0, K_1e=0.0,K_3e=0.0, K_1n=0.0,K_3n=0.0, K_1s=0.0,K_3s=0.0;
	float r_w=0.0,r_e=0.0,r_n=0.0,r_s=0.0;
	float B_EE=0.0,B_WW=0.0,B_NN=0.0,B_SS=0.0,B_P=0.0;

	QUICK_geo_V(j,i,&G_1w,&G_3w,&G_1e,&G_3e,&G_1n,&G_3n,&G_1s,&G_3s,&K_1w,&K_3w,&K_1e,&K_3e,&K_1n,&K_3n,&K_1s,&K_3s,&r_w,&r_e,&r_n,&r_s,F_e,F_w,F_n,F_s);

	A_E = D_e+G_3e*MAX(-F_e,0.0)+G_1w*MAX(-F_w,0.0);

	A_W = D_w+G_1e*MAX(F_e,0.0)+G_3w*MAX(F_w,0.0);

	A_N = D_n+G_3n*MAX(-F_n,0.0)+G_1s*MAX(-F_s,0.0);

	A_S = D_s+G_1n*MAX(F_n,0.0)+G_3s*MAX(F_s,0.0);

	A_P = A_E+A_W+A_N+A_S+F_e*G_3e-F_w*G_3w-F_s*G_3s+F_n*G_3n-G_1e*MAX(F_e,0.0)-G_1n*MAX(F_n,0.0)-G_1w*MAX(0.0,-F_w)-G_1s*MAX(0.0,-F_s);

	S = -r_e*(V_prev[i+1][j]-K_3e*V_prev[i][j]+K_1e*V_prev[i-1][j])*MAX(F_e,0.0)+r_w*(V_prev[i][j]-K_3w*V_prev[i-1][j]
		+K_1w*V_prev[i-2][j])*MAX(F_w,0.0)-r_n*(V_prev[i][j+1]-K_3n*V_prev[i][j]+K_1n*V_prev[i][j-1])*MAX(F_n,0.0)
		+r_s*(V_prev[i][j]-K_3s*V_prev[i][j-1]+K_1s*V_prev[i][j-2])*MAX(F_s,0.0)+r_e*(V_prev[i][j]-K_3e*V_prev[i+1][j]
		+K_1e*V_prev[i+2][j])*MAX(-F_e,0.0)-r_w*(V_prev[i-1][j]-K_3w*V_prev[i][j]+K_1w*V_prev[i+1][j])*MAX(-F_w,0.0)
		+r_n*(V_prev[i][j]-K_3n*V_prev[i][j+1]+K_1n*V_prev[i][j+2])*MAX(-F_n,0.0)-r_s*(V_prev[i][j-1]-K_3s*V_prev[i][j]
		+K_1s*V_prev[i][j+1])*MAX(-F_s,0.0)-G_1w*V[2][i-2][j]*MAX(F_w,0.0)-G_1s*V[2][i][j-2]*MAX(F_s,0.0)-G_1e*V[2][i+2][j]*MAX(-F_e,0.0)
		-G_1n*V[2][i][j+2]*MAX(-F_n,0.0);
}
