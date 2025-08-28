/*UNIVERSAL HEADER FILE FOR DIFFERENT CFD APPLICATIONS

CREATED BY - ORKODIP MOOKHERJEE (08/01/2018)

CONTENTS

1) THOMAS' TDMA ALGORITHM WITH BACKWARD SUBSTITUTION
2) MAXIMUM OF TWO NUMBERS
3) CHECK WHETHER A NUMBER IS GREATER THAN ZERO OR NOT

*/
void tdma_bs(int n, float *a, float *b, float *c, float *d)
{
	//TDMA ALGORITHM STARTS
	for (int i=0;i<n;i++)
	{
		if (i==0)
		{
			c[0]/=b[0];
			d[0]/=b[0];
		}
		else
		{
			c[i]/=(b[i]-a[i]*c[i-1]);
			d[i]=(d[i]-a[i]*d[i-1])/(b[i]-a[i]*c[i-1]);
		}
	}
	//TDMA ALGORITHM ENDS

	//BACKWARD SUBSTITUTION BEGINS

	a[n-1]=d[n-1];
	for (int i=(n-2);i>=0;i--)
		a[i]=(d[i]-c[i]*a[i+1]);

	//BACKWARD SUBSTITUTION ENDS
}

float MAX (float a, float b)
{
	if (a>=b)
		return a;
	else
		return b;
}

float CHK(float a)
{
	if (a>=0.0)
		return 1.0;
	else
		return 0.0;
}
