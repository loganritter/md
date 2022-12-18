#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

void initializeSimulation();
void initializeParameters();
void singleStep();
void initialVelocity();
void initialPosition();
void initialAcceleration();
void initialDiffusion();
void leapfrogStep(int);
void applyBoundary();
void calculateForce();
void evaluateProperty();
void collectProperty(int);
void allocateMemory();
void thermostat();
void buildNeighborList();
void printResults(FILE*);
void printRDF();
void radialDist();
void printDiffusion();
void evaluateDiffusion();
void resetDiffusion();
void collectDiffusion();

typedef struct
{
	double x, y, z;
} VecR;

typedef struct
{
	int x, y, z;
} VecI;

typedef struct
{
	VecR r, rv, ra;
	double mass;
	char type;
} Mol;

typedef struct
{
	double val, sum, sum2;
} Prop;

typedef struct
{
	VecR *orgR, *rTrue;
	double *rrDiffuse;
	int count;
} TBuf;

Prop pressure, KE, E_total, Temp, PE;
VecR region, v_sum, vSum;
VecR rSum;
VecI cell, unit_cell;
Mol *mol;
TBuf *tBuf, *tBufAA,*tBufBB, *tBufAB;
FILE *fp;

#define AllocMem(a, n, t)\
	a = (t*)malloc((n) * sizeof(t))

#define V_add(v1, v2, v3)\
	(v1).x = (v2).x + (v3).x,\
	(v1).y = (v2).y + (v3).y,\
	(v1).z = (v2).z + (v3).z

#define V_sub(v1, v2, v3)\
	(v1).x = (v2).x - (v3).x,\
	(v1).y = (v2).y - (v3).y,\
	(v1).z = (v2).z - (v3).z

#define VV_sub(v1, v2)\
	V_sub(v1, v1, v2)

#define V_dot(v1, v2)\
	((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)

#define Vs_add(v1, v2, s3, v3)\
	(v1).x = (v2).x + (s3) * (v3).x,\
	(v1).y = (v2).y + (s3) * (v3).y,\
	(v1).z = (v2).z + (s3) * (v3).z

# define VCSum(v)\
	((v).x + (v).y + (v).z)

#define V_set(v,sx,sy,sz)\
	(v).x = sx,\
	(v).y = sy,\
	(v).z = sz

#define V_setAll(v, s)\
	V_set(v, s, s, s)

#define V_zero(v)\
	V_setAll(v, 0)

#define V_VS_add(v1, s, v2)\
	Vs_add(v1, v1, s, v2)

#define V_length_sq(v)\
	V_dot(v, v)

#define V_cell_wrap(t)\
	if (m2v.t >= cell.t)\
	{\
		m2v.t = 0;\
		shift.t = region.t;\
	}\
	else if (m2v.t < 0)\
	{\
		m2v.t = cell.t -1;\
		shift.t = -region.t;\
	}

#define V_cell_wrapAll()\
{\
	V_cell_wrap(x);\
	V_cell_wrap(y);\
	V_cell_wrap(z);\
}

#define V_wrap(v, t)\
	if (v.t >= 0.5 * region.t)\
		v.t-= region.t;\
	else if (v.t < -0.5 * region.t)\
		v.t += region.t

#define V_wrapAll(v)\
{\
	V_wrap (v, x);\
	V_wrap(v, y);\
	V_wrap(v, z);\
}

#define sqr(x) (x * x)
#define cube(x) (x * x * x)

#define DO_MOL\
	for(n=0; n<nMol; n++)
	
#define v_Mul(v1, v2, v3)\
	(v1).x = (v2).x * (v3).x,\
	(v1).y = (v2).y * (v3).y,\
	(v1).z = (v2).z * (v3).z

#define V_scale(v, s)\
	(v).x *= s,\
	(v).y *= s,\
	(v).z *= s

#define VS_copy(v2, s1, v1)\
	(v2).x = (s1) * (v1).x,\
	(v2).y = (s1) * (v1).y,\
	(v2).z = (s1) * (v1).z

#define V_prod(v)\
	((v).x * (v).y * (v).z)

#define V_div(v1, v2, v3)\
	(v1).x = (v2).x / (v3).x,\
	(v1).y = (v2).y / (v3).y,\
	(v1).z = (v2).z / (v3).z

#define VV_add(v1, v2)\
	(V_add(v1, v1, v2))

#define prop_zero(v)\
	v.sum = 0.,\
	v.sum2 = 0.

#define prop_accum(v)\
	v.sum += v.val,\
	v.sum2 += sqr(v.val)

#define Max(x1, x2)\
	(((x1)>(x2)) ? (x1):(x2))

#define prop_avg(v, n)\
	v.sum /= n,\
	v.sum2 = sqrt(Max(v.sum2 / n - sqr(v.sum), 0.))

#define prop_est(v)\
	v.sum, v.sum2

#define V_linear(c, r)\
	((c).z * (r).y + (c).y) * (r).x + (c).x

# define N_OFFSET 14
# define OFFSET_VALS \
{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{-1,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1},{-1,1,1},{-1,0,1},{-1,-1,1},\
{0,-1,1},{1,-1,1}}

# define Vprod(c)\
	(c).x * (c).y * (c).z

# define DO_CELL(j, m)\
	for(j = cell_list[m]; j >= 0;j = cell_list[j])
	
#define Nint(x) \
	(((x) < 0.) ? (- (int) (0.5 - (x))): ((int) (0.5 + (x))))

int step_limit;
int step_avg = 100;
double dt = 0.001;
double temperature = 0.6;
double density = 1.2;
int step_count = 0;
double time_now = 0;
int nMol = 864;
double r_cutoff, uSum, virSum, v_magnitude, dist_hist, nebr_shell = 0.4;
int *cell_list;
int adjust_temp = 50, step_equil = 50000;
int *nebr_tab, nebr_tab_len, nebr_tab_max, nebr_now = 1, nebr_tab_fac = 100;
double mass_ratio = 1.0;
double mass2 = 1.0;
double mass1;
double eps_aa = 1.0, eps_ab = 1.5, eps_bb = 0.5;
double sigma_aa = 1.0, sigma_ab = 0.8, sigma_bb = 0.88;
double *rr_diff_avg_aa, *rr_diff_avg_bb, *rr_diff_avg_ab;
int tally_avg_diffuse, diffuse_limit = 1000, nbuff_diffuse = 250, nval_diffuse = 2800, step_diffuse = 100;
int nMolA = 0, nMolB = 0;
double delta_R, latticCorr, dr_region = 4.0, *hisRdfAB, *hisRdfBB, *hisRdfAA;
int size_hist_rdf = 200, count_number = 0, count_limit = 100, stepRdf = 50;
double *MSDAA, *MSDBB, *MSDAB;
double Q;

int main()
{
	clock_t start, finish;
	int simulate, progress;
	double duration;
	
	start = clock();

	initializeSimulation();
	initializeParameters();

	simulate = 1;
	fp = fopen("output.dat","w");
	fprintf(fp, "Steps E_total Pressure Temp U\n");

	progress = floor(step_limit / 10);

	while(simulate)
	{
        if(step_count == progress)
            printf("[ 10 |");
        else if(step_count == 2*progress)
            printf(" 20 |");
        else if(step_count == 3*progress)
            printf(" 30 |");
        else if(step_count == 4*progress)
            printf(" 40 |");
        else if(step_count == 5*progress)
            printf(" 50 |");
        else if(step_count == 6*progress)
            printf(" 60 |");
        else if(step_count == 7*progress)
            printf(" 70 |");
        else if(step_count == 8*progress)
            printf(" 80 |");
        else if(step_count == 9*progress)
            printf(" 90 |");
        else if(step_count == 10*progress - 1)
            printf(" 100 ]\n");
        fflush(stdout);

		singleStep();
		if (step_count >= step_limit)
			simulate = 0;
	}
	fclose(fp);
	
	finish = clock();
	duration = (double)(finish - start)/CLOCKS_PER_SEC;
	printf("Simulation time: %2.3f seconds\n", duration);

	return(0);
}

void initializeSimulation()
{
	int M = 1;

	step_limit = (nval_diffuse + (diffuse_limit - 1) * nval_diffuse / nbuff_diffuse) * step_diffuse + step_equil;
	mass1 = mass2 * mass_ratio;

	while(4*M*M*M < nMol)
		M++;

	unit_cell.x = M;
	unit_cell.y = M;
	unit_cell.z = M;

	r_cutoff = 2.5 * sigma_ab;
	VS_copy(region, 1.0/pow(density/4.0, 1.0/3.0), unit_cell);
	VS_copy(cell, 1.0/(r_cutoff + nebr_shell), region);
	nebr_tab_max = nebr_tab_fac * nMol;

	allocateMemory();
}

void allocateMemory()
{
	int nb;
	
	AllocMem(mol, nMol, Mol);
	AllocMem(cell_list, Vprod(cell) + nMol, int);
	AllocMem(nebr_tab, 2 * nebr_tab_max, int);

	AllocMem(hisRdfAA, size_hist_rdf, double);
	AllocMem(hisRdfAB, size_hist_rdf, double);
	AllocMem(hisRdfBB, size_hist_rdf, double);

	AllocMem(rr_diff_avg_aa, nval_diffuse, double);
	AllocMem(rr_diff_avg_bb, nval_diffuse, double);
	AllocMem(rr_diff_avg_ab, nval_diffuse, double);

	AllocMem(MSDAA, nval_diffuse, double);
	AllocMem(MSDBB, nval_diffuse, double);
	AllocMem(MSDAB, nval_diffuse, double);

	AllocMem(tBufAA, nbuff_diffuse, TBuf);
	AllocMem(tBufBB, nbuff_diffuse, TBuf);
	AllocMem(tBufAB, nbuff_diffuse, TBuf);
	AllocMem(tBuf, nbuff_diffuse, TBuf);

	for(nb = 0; nb < nbuff_diffuse; nb++)
	{
		AllocMem(tBuf[nb].orgR, nMol, VecR);
		AllocMem(tBuf[nb].rTrue, nMol, VecR);
		AllocMem(tBuf[nb].rrDiffuse, nval_diffuse, double);

		AllocMem(tBufAA[nb].orgR, nMol, VecR);
		AllocMem(tBufAA[nb].rTrue, nMol, VecR);
		AllocMem(tBufAA[nb].rrDiffuse, nval_diffuse, double);

		AllocMem(tBufBB[nb].orgR, nMol, VecR);
		AllocMem(tBufBB[nb].rTrue, nMol, VecR);
		AllocMem(tBufBB[nb].rrDiffuse, nval_diffuse, double);

		AllocMem(tBufAB[nb].orgR, nMol, VecR);
		AllocMem(tBufAB[nb].rTrue, nMol, VecR);
		AllocMem(tBufAB[nb].rrDiffuse, nval_diffuse, double);
	}
}

double gaussian()
{
	static int available = 0;
	double fac, r_sq, v1, v2;
	static double gset;
	
	if(!available)
	{
		do
		{
			v1 = 2 * rand()/(double)RAND_MAX -1;
			v2 = 2 * rand()/(double)RAND_MAX -1;
			r_sq = v1*v1+v2*v2;
		} while(r_sq >= 1 || r_sq == 0);

		fac = sqrt(-2.0*log(r_sq)/r_sq);
		gset = v1*fac;
		available = 1;

		return v2*fac;
	}
	else
	{
		available = 0;
		return gset;
	}
}

void initializeParameters()
{
	initialPosition();
	initialVelocity();
	initialAcceleration();
	initialDiffusion();
}

void singleStep()
{
	step_count++;
	time_now = step_count * dt;

	leapfrogStep(1);
	applyBoundary();

	if(nebr_now)
	{
		nebr_now = 0;
		dist_hist = 0;
		buildNeighborList();
	}

	calculateForce();
	leapfrogStep(2);
	evaluateProperty();

	if((step_count < step_equil) && (step_count % adjust_temp == 0))
		thermostat();

	collectProperty(1);

	if((step_count >= step_equil) && ((step_count - step_equil) % stepRdf == 0))
		radialDist();
	if((step_count >= step_equil) && ((step_count - step_equil) % step_diffuse == 0))
		evaluateDiffusion();
	if(step_count % step_avg == 0)
	{
		collectProperty(2);
		printResults(fp);
		collectProperty(0);
	}
}

void initialPosition()
{
	VecR gap, c;
	int nx, ny, nz, n, m, t, d;
	double mratio;

	n = 0;
	V_div(gap, region, unit_cell);

	for(nz = 0; nz < unit_cell.z; nz++)
		for(ny = 0; ny < unit_cell.y; ny++)
			for(nx = 0; nx < unit_cell.x; nx++)
			{
				V_set(c, nx + 0.25, ny + 0.25, nz + 0.25);
				v_Mul(c, c, gap);
				V_VS_add(c, -0.5, region);
			
				for(m = 0; m < 4; m++)
				{
					mol[n].r = c;
					if(m != 3)
					{
						if(m != 0)
							(mol[n].r).x += 0.5 * gap.x;
						if(m != 1)
							(mol[n].r).y += 0.5 * gap.y;
						if(m != 2)
							(mol[n].r).z += 0.5 * gap.z;
					}
					n++;
				}
			}

	for(t = 0; t < nMol; t++)
	{
		if(t % 5 == 0)
		{
			mol[t].type = 'B';
			mol[t].mass = mass2;
		}
		else
		{
			mol[t].type = 'A';
			mol[t].mass = mass1;
		}
	}
	
	for (d = 0; d < nMol; d++)
	{
		if(mol[d].type == 'A')
			nMolA = nMolA+1;
		else if(mol[d].type =='B')
			nMolB = nMolB+1;
	}
	
	mratio = (double)mass1/((double)mass2);
	Q = (double)(sqr((mratio*nMolA + nMolB))/(nMolA*nMolB));
}

void initialVelocity()
{
	int i, n;

	DO_MOL
	{
		mol[n].rv.x = gaussian();
		mol[n].rv.y = gaussian();
		mol[n].rv.z = gaussian();
	}

	double vCM[3] = {0,0,0};
	double m = nMolA*mass1 + nMolB*mass2;

	DO_MOL
	{
		vCM[0] += mol[n].rv.x * mol[n].mass;
		vCM[1] += mol[n].rv.y * mol[n].mass;
		vCM[2] += mol[n].rv.z * mol[n].mass;
	}

	for(i = 0; i < 3; i++)
		vCM[i] /= m;

	DO_MOL
	{
		mol[n].rv.x -= vCM[0];
		mol[n].rv.y -= vCM[1];
		mol[n].rv.z -= vCM[2];
	}

	v_magnitude = sqrt(3*(nMol - 1) * temperature/(nMolA*mass1 + nMolB*mass2));
}

void initialAcceleration()
{
	int n;

	DO_MOL
		V_zero(mol[n].ra);
}

void initialDiffusion()
{
	int nb;

	for(nb = 0; nb < nbuff_diffuse; nb++)
	   	tBuf[nb].count = -nb * nval_diffuse / nbuff_diffuse;

	resetDiffusion();
}

void resetDiffusion()
{
	int j;
	tally_avg_diffuse = 0;

	for(j = 0; j < nval_diffuse; j ++)
	{
		rr_diff_avg_aa[j] = 0.0;
		rr_diff_avg_ab[j] = 0.0;
		rr_diff_avg_bb[j] = 0.0;
	}

	for(j = 0; j < nval_diffuse; j++)
	{
		rr_diff_avg_aa[j] = 0.0;
		rr_diff_avg_ab[j] = 0.0;
		rr_diff_avg_bb[j] = 0.0;
	}
}

void calculateForce()
{
	VecR dr;
	int n, j1, j2;
	double rr, r6, rri, rr_cutoff, fAA, fAB, fBA, fBB, uVal;

	rr_cutoff = sqr(r_cutoff);
	
	DO_MOL
		V_zero(mol[n].ra);
	
	uSum = 0.0;
	virSum = 0.0;
	
	for(n = 0; n < nebr_tab_len; n++)
	{
		j1 = nebr_tab[2*n];
		j2 = nebr_tab[2*n+1];
		V_sub(dr,mol[j1].r,mol[j2].r);
		V_wrapAll(dr);
		rr= V_length_sq(dr);
		
		if(rr < rr_cutoff)
		{
			if(mol[j1].type == 'A' && mol[j2].type == 'A')
			{
				rri = sigma_aa * sigma_aa / rr;
				r6 = cube(rri);

				uVal = 4.0 * eps_aa * r6 * (r6 - 1);
				fAA = 48.0 * eps_aa * (1.0 / (sigma_aa * sigma_aa)) * r6 * rri * (r6 - 0.5);

				V_VS_add(mol[j1].ra, fAA/mass1, dr);
				V_VS_add(mol[j2].ra, -fAA/mass1, dr);

				virSum = virSum + fAA * rr;
				uSum += uVal;
			}

			else if(mol[j1].type == 'A' && mol[j2].type == 'B')
			{
				rri = sigma_ab * sigma_ab / rr;
				r6 = cube(rri);

				uVal = 4.0 * eps_ab * r6 * (r6 - 1);
				fAB = 48.0 * eps_ab * (1. / (sigma_ab * sigma_ab)) * r6 * rri * (r6 - 0.5);

				V_VS_add(mol[j1].ra, fAB/mass1, dr);
				V_VS_add(mol[j2].ra, -fAB/mass2, dr);

				virSum = virSum + fAB * rr;
				uSum += uVal;
			}
			
			else if(mol[j1].type == 'B' && mol[j2].type == 'A')
			{
				rri = sigma_ab * sigma_ab / rr;
				r6 = cube(rri);

				uVal = 4.0 * eps_ab * r6 * (r6 - 1);
				fBA = 48.0 * eps_ab * (1.0 / (sigma_ab * sigma_ab)) * r6 * rri * (r6 - 0.5);

				V_VS_add(mol[j1].ra, fBA/mass2, dr);
				V_VS_add(mol[j2].ra, -fBA/mass1, dr);

				virSum = virSum + fBA * rr;
				uSum += uVal;
			}
			
			else if(mol[j1].type == 'B' && mol[j2].type == 'B')
			{
				rri = sigma_bb*sigma_bb/rr;
				r6 = cube(rri);

				uVal = 4.0 * eps_bb * r6 * (r6 - 1);
				fBB = 48.0 * eps_bb * (1.0 / (sigma_bb * sigma_bb)) * r6 * rri * (r6 - 0.5);

				V_VS_add (mol[j1].ra, fBB/mass2, dr);
				V_VS_add(mol[j2].ra, -fBB/mass2, dr);

				virSum = virSum + fBB * rr;
				uSum += uVal;
			}
		}
	}
}

void radialDist()
{
	VecR dr;
	int n, j1, j2;
	double rr, volume, normalize_AA, normalize_AB, normalize_BB, n_sqr;
	
	if(count_number == 0)
		for(n = 0; n < size_hist_rdf;n++)
		{
			hisRdfAA[n]=0.;
			hisRdfAB[n]=0.;
			hisRdfBB[n]=0.;
		}
	
	delta_R = dr_region / size_hist_rdf;
	
	for(j1 = 0; j1 < nMol - 1; j1++)
		for(j2 = j1 + 1; j2 < nMol; j2++)
		{
			V_sub(dr, mol[j1].r, mol[j2].r);
			V_wrapAll(dr)
				rr = V_length_sq(dr);
			if(rr < sqr(dr_region))
			{
				n = sqrt(rr) / delta_R;
				if(mol[j1].type == 'A' && mol[j2].type =='A')
					hisRdfAA[n]++;
				if(mol[j1].type == 'B' && mol[j2].type =='B')
					hisRdfBB[n]++;
				if((mol[j1].type == 'A' && mol[j2].type == 'B')||(mol[j1].type == 'B' && mol[j2].type == 'A'))
					hisRdfAB[n]++;
			}
		}\

	count_number++;

	if(count_number == count_limit)
	{
		volume = V_prod(region);

		normalize_AA = volume/(2.0 * M_PI * sqr(nMolA) * cube(delta_R) * count_number);
		normalize_AB = volume/(4.0 * M_PI * nMolA* nMolB * cube(delta_R) * count_number);
		normalize_BB = volume/(2.0 * M_PI * sqr(nMolB) * cube(delta_R) * count_number);

		for(n = 0; n < size_hist_rdf; n++)
		{
			n_sqr = (n + 0.5) * (n + 0.5);
			hisRdfAA[n] *= normalize_AA / n_sqr;
			hisRdfAB[n] *= normalize_AB / n_sqr;
			hisRdfBB[n] *= normalize_BB / n_sqr;
		}

		printRDF();
		count_number = 0;
	}
}

void buildNeighborList()
{
	VecR dr, inv_wid, rs, shift;
	VecI c, m1v, m2v, v_off[] = OFFSET_VALS;
	
	int n, m1x, m1y, m1z, m1, m2, offset;
	int cc = 0;
	int j1 = 0;
	int j2 = 0;
	double rrNbr, rr;

	rrNbr = sqr((r_cutoff + nebr_shell));
	V_div(inv_wid, cell, region);

	for(n = nMol; n < nMol + V_prod(cell); n++)
		cell_list[n]= -1;

	DO_MOL
	{
		Vs_add(rs, mol[n].r, 0.5, region);
		v_Mul(c, rs, inv_wid);
		cc = V_linear(c, cell) + nMol;
		cell_list[n]=cell_list[cc];
		cell_list[cc] = n;
	}
	
	nebr_tab_len = 0;
	for(m1z = 0; m1z < cell.z; m1z++)
	{
		for(m1y = 0; m1y<cell.y; m1y++)
		{
			for(m1x = 0; m1x < cell.x; m1x++)
			{
				V_set(m1v, m1x, m1y, m1z);
				m1 = V_linear(m1v, cell)+nMol;

				for(offset = 0; offset < N_OFFSET; offset++)
				{
					V_add(m2v,m1v,v_off[offset]);
					V_zero(shift);
					V_cell_wrapAll();
					m2 = V_linear(m2v, cell) + nMol;
					DO_CELL(j1,m1)
						DO_CELL(j2,m2)
						{
							if(m1 != m2 || j2 < j1)
							{
								V_sub(dr, mol[j1].r, mol[j2].r);
								VV_sub(dr, shift);
								rr= V_length_sq(dr);
								if(rr < rrNbr)
								{
									nebr_tab[2*nebr_tab_len] = j1;
									nebr_tab[2*nebr_tab_len+1] = j2;
									nebr_tab_len++;
								}
							}
						}
				}
			}
		}
	}
}

void leapfrogStep(int iteration)
{
	int n;
	
	if(iteration == 1)
		DO_MOL
		{
			V_VS_add(mol[n].rv, 0.5*dt, mol[n].ra);
			V_VS_add(mol[n].r, dt, mol[n].rv);
		}
	
	else
		DO_MOL
			V_VS_add(mol[n].rv, 0.5*dt, mol[n].ra);
}

void applyBoundary()
{
	int n;

	DO_MOL
		V_wrapAll(mol[n].r);
}

void collectProperty(int checkpoint)
{
	if(checkpoint == 0)
	{
		prop_zero(E_total);
		prop_zero(KE);
		prop_zero(pressure);
		prop_zero(Temp);
		prop_zero(PE);
	}

	else if(checkpoint == 1)
	{
		prop_accum(E_total);
		prop_accum(KE);
		prop_accum(pressure);
		prop_accum(Temp);
		prop_accum(PE);
	}
	
	else if(checkpoint == 2)
	{
		prop_avg(E_total,step_avg);
		prop_avg(pressure, step_avg);
		prop_avg(KE, step_avg);
		prop_avg(Temp,step_avg);
		prop_avg(PE, step_avg);
	}
}

void thermostat()
{
	int n;
	double vScale, vv_sum;
	
	vv_sum = 0.0;

	DO_MOL
		vv_sum += mol[n].mass * V_length_sq(mol[n].rv);
	
	vScale = v_magnitude / sqrt(vv_sum/(nMolA * mass1 + nMolB * mass2));
	
	DO_MOL
		V_scale(mol[n].rv, vScale);
}

void evaluateProperty()
{
	int n;
	double v_sqr;
	double vv_sum = 0.0;
	double vvMax = 0.0;
	
	V_zero(v_sum);

	DO_MOL
	{
		VV_add(v_sum, mol[n].rv);
		v_sqr = V_length_sq(mol[n].rv);
		vv_sum += v_sqr * mol[n].mass;
		vvMax = Max(vvMax,v_sqr);
	}
	
	KE.val = 0.5 * vv_sum / nMol;
	PE.val = uSum/nMol;
	E_total.val = KE.val + PE.val;

	pressure.val = density * (vv_sum + virSum) / (3 * nMol);

	Temp.val = vv_sum / (3 * nMol);

	dist_hist += sqrt(vvMax) * dt;

	if(dist_hist > 0.5 * nebr_shell)
		nebr_now = 1;
}

void printRDF()
{
	int n;
	double rb;
	FILE *fpAA, *fpBB, *fpAB;
	
	fpAA = fopen("AA_RDF.dat","a");
	fpAB = fopen("AB_RDF.dat","a");
	fpBB = fopen("BB_RDF.dat","a");
	
	fprintf(fpAA, "r g(r)\n");
	fprintf(fpAB, "r g(r)\n");
	fprintf(fpBB, "r g(r)\n");
	
	for(n = 0; n < size_hist_rdf; n++)
	{
		rb = (n + 0.5) * dr_region / size_hist_rdf;
		
		fprintf(fpAA, "%8.4f %8.4f\n", rb, hisRdfAA[n]);
		fprintf(fpAB, "%8.4f %8.4f\n", rb, hisRdfAB[n]);
		fprintf(fpBB, "%8.4f %8.4f\n", rb, hisRdfBB[n]);
	}
	
	fclose(fpAA);
	fclose(fpAB);
	fclose(fpBB);
}

void evaluateDiffusion()
{
	VecR dr;
	int n, nb, ni;
	
	for(nb = 0; nb < nbuff_diffuse; nb ++)
	{
		if(tBuf[nb].count == 0)
		{
			DO_MOL
			{
				if(mol[n].type == 'A')
				{
					tBufAA[nb].orgR[n] = mol[n].r;
					tBufAA[nb].rTrue[n] = mol[n].r;
				}
				if(mol[n].type =='B')
				{
					tBufBB[nb].orgR[n] = mol[n].r;
					tBufBB[nb].rTrue[n] = mol[n].r;
				}
			}
		}

		if(tBuf[nb].count >= 0)
		{
			ni = tBuf[nb].count;
			tBufAA[nb].rrDiffuse[ni] = 0.;
			tBufBB[nb].rrDiffuse[ni] = 0.;
			tBufAB[nb].rrDiffuse[ni] = 0.;
			V_setAll(rSum, 0.);

			DO_MOL{
				if (mol[n].type == 'A')
				{
					V_sub(dr, tBufAA[nb].rTrue[n], mol[n].r);
					V_div(dr, dr, region);

					dr.x = Nint(dr.x);
					dr.y = Nint(dr.y);
					dr.z = Nint(dr.z);

					v_Mul(dr, dr, region);
					V_add(tBufAA[nb].rTrue[n], mol[n].r, dr);
					V_sub(dr, tBufAA[nb].rTrue[n], tBufAA[nb].orgR[n]);

					VV_add(rSum, dr);
					tBufAA[nb].rrDiffuse[ni] += V_length_sq(dr);
				}
				
				if (mol[n].type == 'B')
				{
					V_sub(dr, tBufBB[nb].rTrue[n], mol[n].r);
					V_div(dr, dr, region);

					dr.x = Nint(dr.x);
					dr.y = Nint(dr.y);
					dr.z = Nint(dr.z);

					v_Mul(dr, dr, region);
					V_add(tBufBB[nb].rTrue[n], mol[n].r, dr);
					V_sub(dr, tBufBB[nb].rTrue[n], tBufBB[nb].orgR[n]);

					tBufBB[nb].rrDiffuse[ni] += V_length_sq(dr);
				}
			}

			tBufAB[nb].rrDiffuse[ni] = V_length_sq(rSum);
		}

		tBuf[nb].count++;
	}

	collectDiffusion();
}

void collectDiffusion()
{
	int nb, j;
	double facAA, facBB, facAB;
	
	for(nb = 0; nb < nbuff_diffuse; nb++)
	{
		if(tBuf[nb].count == nval_diffuse)
		{
			for(j = 0; j < nval_diffuse; j++)
			{
				rr_diff_avg_aa[j] += tBufAA[nb].rrDiffuse[j];
				rr_diff_avg_bb[j] += tBufBB[nb].rrDiffuse[j];
				rr_diff_avg_ab[j] += tBufAB[nb].rrDiffuse[j];
			}
			
			tBuf[nb].count = 0;
			tally_avg_diffuse++;
			
			if(tally_avg_diffuse == diffuse_limit)
			{
				for(j = 1; j < nval_diffuse; j++)
				{
					MSDAA[j] = rr_diff_avg_aa[j] / (diffuse_limit * nMolA);
					MSDBB[j] = rr_diff_avg_bb[j] / (diffuse_limit * nMolB);
					MSDAB[j] = rr_diff_avg_ab[j] / (diffuse_limit * nMol);
				}
				facAA = 1.0 / (6 * nMolA * step_diffuse * dt * diffuse_limit);
				facBB = 1.0 / (6 * nMolB * step_diffuse * dt * diffuse_limit);
				facAB = (Q/ (6 * nMol * step_diffuse * dt * diffuse_limit));

				for(j = 1; j < nval_diffuse; j++)
				{
					rr_diff_avg_aa[j] *= facAA / j;
					rr_diff_avg_bb[j] *= facBB / j;
					rr_diff_avg_ab[j] *= facAB / j;
				}

				printDiffusion();
				resetDiffusion();
			}
		}
	}
}

void printResults(FILE *fp)
{
	fprintf(fp,"%5d %7.4f %7.4f %7.4f %7.4f\n", step_count, E_total.val, pressure.sum, Temp.sum, PE.sum);
}

void printDiffusion()
{
	int j;
	double tVal;
	FILE *fpDaa, *fpDbb, *fpDab;
	
	fpDaa = fopen("AA_diffusion.dat","w");
	fpDbb = fopen("BB_diffusion.dat","w");
	fpDab = fopen("AB_diffusion.dat","w");
	
	fprintf(fpDaa, "time AA_Diffusion MSD_AA\n");
	fprintf(fpDbb, "time BB_Diffusion MSD_BB\n");
	fprintf(fpDab, "time AB_Diffusion MSD_AB\n");
	
	for (j = 0; j < nval_diffuse; j++)
	{
		tVal = j * step_diffuse * dt;
		fprintf(fpDaa,"%8.4f %8.4f %8.4f\n", tVal, rr_diff_avg_aa[j], MSDAA[j]);
		fprintf(fpDbb,"%8.4f %8.4f %8.4f\n", tVal, rr_diff_avg_bb[j], MSDBB[j]);
		fprintf(fpDab,"%8.4f %8.4f %8.4f\n", tVal, rr_diff_avg_ab[j], MSDAB[j]);
	}
	
	fclose(fpDaa);
	fclose(fpDab);
	fclose(fpDbb);
}
