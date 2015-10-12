#include "mex.h"
#include <math.h>

//從已有的一維陣列創造 M[W][H]的二維陣列
//用完後直接free(M)即可
double** TransformM2D(double* in, int h, int w)
{
	double **out = (double**)malloc(sizeof(double*)*w);
	out[0] = in;
	for (int i = 1; i < w; i++)
		out[i] = out[i - 1] + h;
	return out;
}

// 建立新的二維陣列 M[W][H]
// 用完後需呼叫FreeM2D(M)
double** CreateM2D(int h, int w)
{
	double **out = (double**)malloc(sizeof(double*)*w);
	out[0] = (double*)malloc(sizeof(double)*w*h);
	for (int i = 1; i < w; i++)
		out[i] = out[i - 1] + h;
	return out;
}

void FreeM2D(double** m)
{
	free(m[0]);
	free(m);
}

double wij(double** g, int h, int w, int yi, int xi, int yj, int xj, int r, double sigmaI, double sigmaX)
{
	double dist,diff;
	if (xi < 0 || xi >= w || yi <0 || yi >= h || xj < 0 || xj >= w || yj <0 || yj >= h)
		return 0;
	dist = sqrt((double)(xi - xj)*(xi - xj) + (yi - yj)*(yi - yj));
	if (dist>r)
		return 0;
	diff = g[xi][yi] - g[xj][yj];
	if (diff < 0)
		diff = -diff;
	return exp(-diff / sigmaI)*exp(-dist / sigmaX);
}

// 輸入：二維灰階圖、Initial Guess(-1->outside 1->inside)、參數vector(r,sigmaI,sigmaX,itermax)
// 輸出：segmentation result(-1->ouside 1->inside)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int r,itermax,H,W,i,j,m,n;
	double sigmaI, sigmaX;
	double* pd;
	double td,tsum;
	double **G, **InitGuess, **IG, **R, **D;
	pd = mxGetPr(prhs[2]);
	r = (int)(*pd);
	sigmaI = *(pd + 1);
	sigmaX = *(pd + 2);
	itermax = (int)(*(pd + 3));
	H = mxGetM(prhs[0]);
	W = mxGetN(prhs[0]);
	G = TransformM2D(mxGetPr(prhs[0]), H, W);
	InitGuess = TransformM2D(mxGetPr(prhs[1]), H, W);
	IG = CreateM2D(H, W);
	memcpy(IG[0], InitGuess[0], sizeof(double)*H*W);
	plhs[0] = mxCreateDoubleMatrix(H, W, mxREAL);
	R = TransformM2D(mxGetPr(plhs[0]), H, W);
	D = CreateM2D(H, W);
	double PWP = 0, PWN = 0, NWN = 0;
	double MWP, MWN, MWM, PWP2, PWN2, NWN2;
	double b, b2, E, E2;
	double DP = 0, DN = 0;
	int iter = 0;
	bool flag = true, flag2;

	for (i = 0; i < H; i++)
	for (j = 0; j < W; j++)
	{
		tsum = 0;
		for (m = -r; m <=r;m++)
		for (n = -r; n <= r; n++)
		{
			if (i + m < 0 || i + m >= H || j + n < 0 || j + n >= W)
				continue;
			td = wij(G, H, W, i, j, i + m, j + n, r, sigmaI, sigmaX);
			tsum += td;
			if (IG[j][i] == 1)
			{
				if (IG[j + n][i + m] == 1)
					PWP += td;
				else
					PWN -= td;
			}
			else
			{
				if (IG[j + n][i + m] == -1)
					NWN += td;
			}
		}
		D[j][i] = tsum;
	}
	
	while (iter < itermax&&flag == true)
	{
		iter++;
		flag = false;
		//在開始更新前計算目前IG下的DP、DN、b、E
		PWP2 = PWP;
		PWN2 = PWN;
		NWN2 = NWN;
		DP = 0;
		DN = 0;
		for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
		{
			if (IG[j][i] == 1)
				DP += D[j][i];
			else
				DN += D[j][i];
		}
		b = DP / DN;
		E = (PWP + b * 2 * PWN + b*b*NWN) / (DP + b*b*DN);
		//在IG中找出區域的邊界做為growing或shrinking的candidate，並進行更新運算
		//將更新的結果放入R中，並同步更新PWP2、PWN2、NWN2

		
		for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
		{
			R[j][i] = IG[j][i];
			
			if (IG[j][i] == -1)
			{
				flag2 = false;
				for (m = -1; m <= 1; m++)
				for (n = -1; n <= 1; n++)
				{
					if (i + m < 0 || i + m >= H || j + n < 0 || j + n >= W)
						continue;
					if (IG[j + n][i + m] == 1) //growing candidate
						flag2 = true;
				}
				if (flag2 == true)
				{
					MWP = 0;
					MWM = wij(G, H, W, i, j, i, j, r, sigmaI, sigmaX);
					MWN = 0;
					for (m = -r; m <= r; m++)
					for (n = -r; n <= r; n++)
					{
						if (i + m < 0 || i + m >= H || j+n < 0 || j + n >= W)
							continue;
						td = wij(G, H, W, i, j, i + m, j + n, r, sigmaI, sigmaX);
						if (IG[j + n][i + m] == 1)
							MWP += td;
						else
							MWN -= td;
					}
					
					b2 = (DP + D[j][i]) / (DN - D[j][i]);
					E2 = (PWP + b2 * 2 * PWN + b2*b2*NWN + (b2 + 1) * 2 * MWP + (b2*b2 + b2 * 2 + 1)*MWM + (b2 + 1)*b2 * 2 * MWN) / ((DP + D[j][i]) + b2*b2*(DN - D[j][i]));
					//printf("%d %d %lf %lf %lf %lf %lf %lf %lf\n", j, i, E, E2,MWP,MWM,MWN,b,b2);
					if (E2>E)
					{
						printf("更新%d %d(growing)\n", j, i);
						R[j][i] = 1;
						flag = true;
						PWP2 = PWP2 + MWP * 2 + MWM;
						PWN2 = PWN2 + MWP + MWN + MWM;
						NWN2 = NWN2 + MWN * 2 + MWM;
					}
					
				}
			}
			else
			{
				flag2 = false;
				for (m = -1; m <= 1; m++)
				for (n = -1; n <= 1; n++)
				{
					if (i + m < 0 || i + m >= H || j+n < 0 || j + n >= W)
						continue;
					if (IG[j + n][i + m] == -1) //shrinking candidate
						flag2 = true;
				}
				if (flag2 == true)
				{
					MWP = 0;
					MWM = wij(G, H, W, i, j, i, j, r, sigmaI, sigmaX);
					MWN = 0;
					for (m = -r; m <= r; m++)
					for (n = -r; n <= r; n++)
					{
						if (i + m < 0 || i + m >= H || j+n < 0 || j + n >= W)
							continue;
						td = wij(G, H, W, i, j, i + m, j + n, r, sigmaI, sigmaX);
						if (IG[j + n][i + m] == 1)
							MWP -= td;
						else
							MWN += td;
					}
					b2 = (DP - D[j][i]) / (DN + D[j][i]);
					E2 = (PWP + b2 * 2 * PWN + b2*b2*NWN + (b2 + 1) * 2 * MWP + (b2*b2 + b2 * 2 + 1)*MWM + (b2 + 1)*b2 * 2 * MWN) / ((DP - D[j][i]) + b2*b2*(DN + D[j][i]));
					if (E2>E)
					{
						R[j][i] = -1;
						flag = true;
						printf("更新%d %d(shrinking)\n", j, i);
						PWP2 = PWP2 + MWP * 2 + MWM;
						PWN2 = PWN2 + MWP + MWN + MWM;
						NWN2 = NWN2 + MWN * 2 + MWM;
					}
				}
			}
		}
		//把結果R copy給IG，以進行下一次運算
		memcpy(IG[0], R[0], sizeof(double)*H*W);
		PWP = PWP2;
		PWN = PWN2;
		NWN = NWN2;
	}
	FreeM2D(D);
	FreeM2D(IG);
	free(G);
	free(InitGuess);
	free(R);
}