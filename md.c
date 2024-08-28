#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define kB 1.38064852e-23						//Bolztmann's constant (J/K)
#define NA 6.02214e23							//Avogadro's constant

#define N_STEPS 10000							//Number of simulation steps
#define BIN 100								//Binning number for radial distribution histogram
#define DIMENSIONS 3							//Number of coordinates (x,y,z)
#define N 256								//Number of particles
#define s 3.405								//Sigma of gas (in Angstroms)
#define p_star 1.1							//Density of gas in (dimentionless)
#define T_star 2.7							//Temperature of gas (dimensionless)
#define dt_star 0.001							//Time step (dimensionless)
#define MASS 39.948							//Mass of gas (amu)
#define epsilon 1.65401693e-21						//Epsilon of gas (J)
#define THERMO (N_STEPS / 2.0)						//Thermostat for specified portion of N_STEPS (i.e. first quarter, half, etc.)
#define APPLY_T 5							//Apply thermostat every X number of steps
static char atom[] = "Ar";						//Atom type

/* THE ABOVE PARAMETERS CAN BE ADJUSTED TO CHANGE THE PROPERTIES OF THE GAS IN THE MACROS BELOW							*
 * Lennard-Jones parameters (epsilon and sigma) taken from: https://www.sciencedirect.com/science/article/pii/002199917590042X?via%3Dihub 	*/

#define p (p_star / pow(s,3.0))						//Density (A^-3)
#define L pow((N/p),1.0/3.0)						//Box length (A)
#define e (epsilon / kB)						//Energy of gas (K)
#define m (MASS * 10.0 / NA / kB)					//Conversion of amu to K*ps^2/A^2
#define T (T_star * e)							//Conversion of temperature to K
#define dt (dt_star * sqrt((m*pow(s,2.0)) / e))				//Conversion of time step to ps

static double r[N][DIMENSIONS];
static double v[N][DIMENSIONS];
static double a[N][DIMENSIONS];

//New array declaration for multicomponent system for sigma and epsilon parameters
//static double sigma[N], epsilon[N];

void crystalLattice();
void wrapToBox();
void initializeVelocities();
double calculateAcceleration();
void velocityVerlet();
double potentialEnergy();
double kineticEnergy();
double generateGaussian();
double meanSquaredVelocity();
void thermostat();
void radialDist(FILE *fp);
void MSD(FILE *fp);
//void VACF(FILE *fp);

int main()
{
	int i,j,k,n;
	int progress;
	double Temp,Press,Pavg,Tavg,V,PE,KE,sqPE,sqKE,PEavg,KEavg,mv,ulrc,plrc,cvPE,cvKE,vSum,v2;

	clock_t start,end;
	double sim_time;

        FILE *ftraj, *fvel, *fener, *fmv, *ftemp, *fpress;
        ftraj = fopen("traj.xyz","w");
	fvel = fopen("velocities.xyz","w");
        fener = fopen("energy.dat","w");
        fmv = fopen("momentum.dat","w");
        ftemp = fopen("temp.dat","w");
        fpress = fopen("pressure.dat","w");

	start = clock();

	crystalLattice();
	calculateAcceleration();
	v2 = meanSquaredVelocity();

	Pavg = 0;
	Tavg = 0;
	PEavg = 0;
	KEavg = 0;
	sqPE = 0;
	sqKE = 0;
	vSum = 0.0;

	progress = floor(N_STEPS / 10);

	for(n=0;n<=N_STEPS;n++)
	{
		wrapToBox();

		if(n == progress)
			printf("[ 10 |");
		else if(n == 2*progress)
			printf(" 20 |");
		else if(n == 3*progress)
			printf(" 30 |");
		else if(n == 4*progress)
			printf(" 40 |");
		else if(n == 5*progress)
			printf(" 50 |");
		else if(n == 6*progress)
			printf(" 60 |");
		else if(n == 7*progress)
			printf(" 70 |");
		else if(n == 8*progress)
			printf(" 80 |");
		else if(n == 9*progress)
			printf(" 90 |");
		else if(n == 10*progress)
			printf(" 100 ]\n");
		fflush(stdout);

		V = calculateAcceleration();

		//Apply thermostat
		if(n != 0 && n % APPLY_T == 0 && n < THERMO)
			thermostat();

		//Momentum at each step
		for(i=0;i<N;i++)
			for(j=0;j<3;j++)
				vSum += v[i][j] / N;
		mv = m * vSum;
	       	fprintf(fmv,"%d\t %lf\n",n,mv);

		//Write atom position to trajectory file
	        for(i=0;i<1;i++)
			{
					fprintf(ftraj,"%d\n\n",N);
					for(j=0;j<N;j++)
					{
			fprintf(ftraj,"%s\t",atom);

							for(k=0;k<3;k++)
									fprintf(ftraj,"%lf\t",r[j][k]);
							fprintf(ftraj,"\n");
					}
			}
		
		velocityVerlet();

		//Write atom velocities to file
                for(i=0;i<1;i++)
                {
                        fprintf(fvel,"%d\n\n",N);
                        for(j=0;j<N;j++)
                        {
				fprintf(fvel,"%s\t",atom);

                                for(k=0;k<3;k++)
                                        fprintf(fvel,"%lf\t",v[j][k]);
                                fprintf(fvel,"\n");
                        }
                }

		//Collect simulation parameters after thermostat switches off
		if(n >= THERMO)
		{
			PE = potentialEnergy();
			KE = kineticEnergy();
			fprintf(fener,"%lf\t %lf\t %lf\n",PE,KE,PE+KE);
			sqPE += PE*PE;
			sqKE += KE*KE;
			PEavg += PE;
                        KEavg += KE;

			//Temperature from kinetic theory of gases (kB = 1 in reduced units)
			Temp = (2.0/3.0) * KE;
			fprintf(ftemp,"%d\t %lf\n",n,Temp);
			Tavg += Temp;

			//Pressure from virial theorem (kB = 1 in reduced units)
			Press = (p*s*s*s)/N/3.0 * (v2 + V/e);
			fprintf(fpress,"%d\t %lf\n",n,Press);
			Pavg += Press;
		}
	}
	fclose(ftraj);
	fclose(fvel);

        printf("*****************************************************************************\n");
        printf("Calculating Radial Distribution...\n");
        radialDist(ftraj);
        printf("Done!\n");

        printf("*****************************************************************************\n");
        printf("Calculating Mean Squared Displacement...\n");
        MSD(ftraj);
        printf("Done!\n");
/*
	printf("*****************************************************************************\n");
        printf("Calculating Velocity Autocorrelation...\n");
        VACF(fvel);
        printf("Done!\n");
*/
	sqPE /= (N_STEPS - THERMO);
	sqKE /= (N_STEPS - THERMO);
	PEavg /= (N_STEPS - THERMO);
	KEavg /= (N_STEPS - THERMO);
	Tavg /= (N_STEPS - THERMO);
	Pavg /= (N_STEPS - THERMO);

	//Specific heat from fluctuation in average potential/kinetic energy and average temp
	cvPE = (3.0/2.0) * pow(1.0 - 2*N*(sqPE-pow(PEavg,2.0))/3.0/pow(Tavg,2.0),-1.0);
	cvKE = (3.0/2.0) * pow(1.0 - 2*N*(sqKE-pow(KEavg,2.0))/3.0/pow(Tavg,2.0),-1.0);

	//Long-range correction to pressure (from Comp. Sim. of Liquids - Allen/Tildesly, pg. 65)
	plrc = (16.0/3.0) * M_PI * pow(p_star,2.0) * ((2.0/3.0) * pow(s/L/2.0,9.0) - pow(s/L/2.0,3.0));

	//Long-range correction to energy (from Comp. Sim. of Liquids - Allen/Tildesly, pg. 65)
	ulrc = (8.0/3.0) * M_PI * N * p_star * e * ((1.0/3.0) * pow(s/L/2.0,9.0) - pow(s/L/2.0,3.0));

	printf("*****************************************************************************\n");
	printf("Momentum: %lf\n",mv);
	printf("*****************************************************************************\n");
	printf("Starting Temperature (K): %lf\n",T);
        printf("Average Temperature (K): %lf\n",Tavg);
	printf("Percent Difference: %.2lf%%\n",(Tavg-T)/T*100.0);
	printf("*****************************************************************************\n");
	printf("--- PRESSURE ---\n");
	printf("Average Reduced Pressure: %lf\n",Pavg);
	printf("Pressure from Long-Range Correction: %lf\n",plrc);
	printf("Average Reduced Pressure with Long-Range Correction: %lf\n",(Pavg+plrc));
	printf("*****************************************************************************\n");
	printf("--- ENERGY ---\n");
	printf("Average Reduced Potential Energy: %lf\n",PEavg/e);
	printf("Energy from Long-Range Correction per Particle: %lf\n",ulrc/e/N);
	printf("Average Reduced Potential Energy with Long-Range Correction: %lf\n",(PEavg+ulrc/N)/e);
	printf("*****************************************************************************\n");
	printf("---SPECIFIC HEAT AT CONSTANT VOLUME---\n");
	printf("Specific Heat from Fluctuations in Potenial Energy: %lf\n",cvPE);
	printf("Specific Heat from Fluctuations in Kinetic Energy: %lf\n",cvKE);
	printf("*****************************************************************************\n");

	end = clock();

	sim_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("SIMULATION RUNTIME: %.2lf seconds\n",sim_time);
	printf("*****************************************************************************\n");

	return(0);
}
/*
// *** UNCOMMENT FOR SIMPLE CUBIC *** //
void crystalLattice()
{
	int i,j,k,n;

	//Number of atoms in each x,y,z direction rounded up to the nearest whole number
	n = ceil(pow(N,1.0/3.0));

	//Atom spacing in given x,y,z direction
	double dr = L / n;

	//Index for total number of atoms, N
	int index = 0;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			for(k=0;k<n;k++)
			{
				r[index][0] = i * dr;
				r[index][1] = j * dr;
				r[index][2] = k * dr;

				index++;
			}

	initializeVelocities();

	//Print initial positions and velocities
//	for(i=0;i<N;i++)
//		printf("%lf\t %lf\t %lf\n",r[i][0],r[i][1],r[i][2]);
//	for(i=0;i<N;i++)
//		printf("%lf\t %lf\t %lf\n",v[i][0],v[i][1],v[i][2]);

}
*/
void crystalLattice()
{
	int i,j,k;
	int index = 0;
	double n,dr,dr2;

	n = pow(N/4.0,1.0/3.0);
	dr = L / n;
	dr2 = dr / 2.0;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			for(k=0;k<n;k++)
			{
				if(index < N)
				{
					r[index][0] = i * dr;
					r[index][1] = j * dr;
					r[index][2] = k * dr;
				}
				index++;

				if(index < N)
				{
					r[index][0] = i * dr + dr2;
					r[index][1] = j * dr + dr2;
					r[index][2] = k * dr;
				}
				index++;

				if(index < N)
				{
					r[index][0] = i * dr;
					r[index][1] = j * dr + dr2;
					r[index][2] = k * dr + dr2;
				}
				index++;
	
				if(index < N)
				{
					r[index][0] = i * dr + dr2;
					r[index][1] = j * dr;
					r[index][2] = k * dr + dr2;
				}
				index++;
			}

	initializeVelocities();

//        for(i=0;i<N;i++)
//        	printf("%lf\t %lf\t %lf\n",r[i][0],r[i][1],r[i][2]);
}

void wrapToBox()
{
	int i,j;

	for (i=0;i<N;i++)
		for (j=0;j<3;j++)
			r[i][j] += -L * floor(r[i][j]/L);
}

void initializeVelocities()
{
	int i,j;
	double vScale;
	double v2 = 0;
	double instant_temp = 0;
	double vCM[3] = {0,0,0};

	//Assign random initial velocities
	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			v[i][j] = generateGaussian();

	//Calculating center of mass velocity
        for(i=0;i<N;i++)
                for(j=0;j<3;j++)
                        vCM[j] += v[i][j];

        for(i=0;i<3;i++)
                vCM[i] /= N;

        //Subtract off center of mass velocity
        for(i=0;i<N;i++)
                for(j=0;j<3;j++)
                        v[i][j] -= vCM[j];

	//Scale initial velocity with scaled temperature (i.e. initial temperature from T* against temperature "right now")
	//instant_temp = m * v^2 / ((3N - 3) * kB) -- kB = 1 in reduced units
        for(i=0;i<N;i++)
                for(j=0;j<3;j++)
                        v2 += v[i][j] * v[i][j] / N;

        for(i=0;i<N;i++)
                instant_temp += m * v2;

        instant_temp /= (3.0*N - 3.0);

        vScale = sqrt(T / instant_temp);
	
	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			v[i][j] *= vScale;
}

double calculateAcceleration()
{
	int i,j,k;
	double F,r2,r6,r12,sor;
	double V = 0;

	//Position of i relative to j
	double rij[3];

	//Initialize acceleration to 0
	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			a[i][j] = 0;

	//Loop over all distinct pairs i,j
	for(i=0;i<N-1;i++)
		for(j=i+1;j<N;j++)
		{
			r2 = 0.0;
			for(k=0;k<3;k++)
			{
				//Component-by-componenent position of i relative to j
				rij[k] = r[i][k] - r[j][k];

				//Periodic boundary conditions
				rij[k] += -L + trunc(rij[k]/L);
	                        while(rij[k] >= 0.5*L)
					rij[k] -= L;
				while(rij[k] < -0.5*L)
	            	        	rij[k] += L;

				//Dot product of the position component
				r2 += rij[k] * rij[k];
			}

			if (r2 < 0.25*L*L)
			{
				sor = (s*s) / r2;
				r6 = sor * sor * sor;
				r12 = r6 * r6;

				F = 48.0 * e/r2 * (r12 - 0.5*r6);

				//Virial sum for pressure calculation
				V += F * r2;

				for(k=0;k<3;k++)
			        {
			                a[i][k] += rij[k] * F/m;
			                a[j][k] -= rij[k] * F/m;
			        }
			}
		}

	return V;
}

void velocityVerlet()
{
	int i,j;

	//Calculate acceleration from forces at current position
	calculateAcceleration();

	//Update positions and velocity with current velocity and acceleration
	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
		{
			r[i][j] += v[i][j]*dt + 0.5*a[i][j]*pow(dt,2.0);
			v[i][j] += 0.5*a[i][j]*dt;
		}
		
	//Update acceleration from current position
	calculateAcceleration();

	//Update velocity with updated acceleration
	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			v[i][j] += 0.5*a[i][j]*dt;
}

double potentialEnergy()
{
	int i,j,k;
	double sor,r2,r6,r12;
	double PE = 0.0;
	double rij[3];

	for(i=0;i<N-1;i++)
		for(j=i+1;j<N;j++)
		{
			r2 = 0.0;
	                for(k=0;k<3;k++)
			{
				rij[k] = r[i][k] - r[j][k];
	
				rij[k] += -L + trunc(rij[k]/L);
	                        while(rij[k] >= L/2.0)
					rij[k] -= L;
				while(rij[k] < -L/2.0)
	            	        	rij[k] += L;

	                	r2 += rij[k] * rij[k];
			}

			if(r2 < 0.25*L*L)
			{	
				sor = (s*s) / r2;
				r6 = sor * sor * sor;
				r12 = r6 * r6;

				PE += 4.0 * e * (r12 - r6) / N;
			}
		}

	return PE;
}

double kineticEnergy()
{
	int i,j;
	double v2 = 0.0;
	double KE = 0.0;

	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			v2 += v[i][j] * v[i][j];

	KE += 0.5 * m * v2 / N;

	return KE;
}

//Generation of random number sampling via the Marsaglia polar method (shamelessly stolen from Wikipedia: https://en.wikipedia.org/wiki/Marsaglia_polar_method#Implementation)
double generateGaussian()
{
	static bool available = false;
	static double gset;
	double fac,rsq,v1,v2;

	if(!available)
	{
		do
		{
			v1 = 2.0 * rand() / (double)RAND_MAX - 1.0;
			v2 = 2.0 * rand() / (double)RAND_MAX - 1.0;
			rsq = v1*v1 + v2*v2;
		} while(rsq >= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0 * log(rsq)/rsq);
		gset = v1 * fac;
		available = true;

		return v2*fac;
	}

	else
	{
		available = false;
		return gset;
	}
}

double meanSquaredVelocity()
{
	int i;
	double v2;
	double vx2 = 0.0;
	double vy2 = 0.0;
	double vz2 = 0.0;

	for(i=0;i<N;i++)
	{
		vx2 += v[i][0]*v[i][0];
		vy2 += v[i][1]*v[i][1];
		vz2 += v[i][2]*v[i][2];
	}

	v2 = (vx2+vy2+vz2) / N;

	return v2;
}

void thermostat()
{
	int i,j;
        double vScale;
	double v2 = 0.0;
        double instant_temp = 0.0;

	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			v2 += v[i][j] * v[i][j] / N;

        for(i=0;i<N;i++)
                instant_temp += m * v2;

        instant_temp /= (3.0*N - 3.0);

        vScale = sqrt(T / instant_temp);

        for(i=0;i<N;i++)
                for(j=0;j<3;j++)
                        v[i][j] *= vScale;
}

void radialDist(FILE *fp)
{
        int i,j,n;
        int row;
        double dx,dy,dz,r,dr,rij;
        double g[BIN];
        double x[N],y[N],z[N];

        FILE *fgr, *fxyz;

        dr = L / 2.0 / BIN;

        for(i=0;i<BIN;i++)
                g[i] = 0.0;

        fxyz = fopen("traj.xyz","r");
        for(n=0;n<N_STEPS;n++)
        {
                //Assign the first row in trajectory file as -1
                for(row=-1;row<N;row++)
                {
                        if(row == -1)
                                fscanf(fxyz,"%*d\n\n"); //Ignore the atom type and empty space beneath
                        else
                                fscanf(fxyz,"%*s %lf %lf %lf\n",&x[row],&y[row],&z[row]); //Each array contains the coordinates for N atoms at each x,y,z position. For unknown reasons, this method ignores the final block of coordinates, but this is deemed insignificant for calculating g(r)
                }

                if(n >= THERMO)
                {
                        for(i=0;i<N-1;i++)
                                for(j=i+1;j<N;j++)
                                {
					//Apply boundary conditions after thermostating
                                        dx = x[i] - x[j];
                                        dx += -L * round(dx/L);
					//while(dx >= L/2.0)
					//	dx -= L;
					//while(dx < -L/2.0)
					//	dx += L;

                                        dy = y[i] - y[j];
                                        dy += -L * round(dy/L);
					//while(dy >= L/2.0)
                                        //        dy -= L;
                                        //while(dy < -L/2.0)
                                        //        dy += L;

                                        dz = z[i] - z[j];
                                        dz += -L * round(dz/L);
					//while(dz >= L/2.0)
                                        //        dz -= L;
                                        //while(dz < -L/2.0)
                                        //        dz += L;

                                        rij = sqrt(dx*dx + dy*dy + dz*dz);
					
					if(rij < L/2.0)
					{
	                                        g[(int)(rij/dr)] += 2.0;
					}
                                }
                }
        }
	fclose(fxyz);

        fgr = fopen("radial_dist.dat","w");
        fprintf(fgr,"r\t\t g(r)\n");
        for(i=0;i<BIN;i++)
        {
                r = (i+0.5) * dr;
                g[i] /= THERMO;
                g[i] /= 4.0*M_PI/3.0 * (pow(i+1.0,3.0) - pow(i,3.0)) * pow(dr,3.0) * p;
                fprintf(fgr,"%lf\t%lf\n",r,g[i]/N);
        }
}

void MSD(FILE *fp)
{
    int i,j,n;
    double dx0, dy0, dz0;
    double dx, dy, dz;
    double r2_0, r2_t, msd;
    double x0[N], y0[N], z0[N];
    double x[N], y[N], z[N];

    FILE *fxyz, *fmsd;

    fmsd = fopen("msd.dat", "w");
    fprintf(fmsd, "Time (ps)\t MSD (A^2)\n");

    fxyz = fopen("traj.xyz", "r");
    for (n = 0; n < N_STEPS; n++)
    {
        if (n == THERMO)
        {
            for (int row = -1; row < N; row++)
            {
                if (row == -1)
                    fscanf(fxyz, "%*d\n\n");
                else
                    fscanf(fxyz, "%*s %lf %lf %lf\n", &x0[row], &y0[row], &z0[row]);
            }
        }

        for (int row = -1; row < N; row++)
        {
            if (row == -1)
                fscanf(fxyz, "%*d\n\n");
            else
                fscanf(fxyz, "%*s %lf %lf %lf\n", &x[row], &y[row], &z[row]);
        }

        if (n >= THERMO)
        {
            double sum_r2 = 0.0;
            int count = 0;
            for (i = 0; i < N - 1; i++)
            {
                for (j = i + 1; j < N; j++)
                {
                    dx0 = x0[i] - x0[j];
                    dy0 = y0[i] - y0[j];
                    dz0 = z0[i] - z0[j];
                    dx0 -= L * round(dx0 / L);
                    dy0 -= L * round(dy0 / L);
                    dz0 -= L * round(dz0 / L);
                    r2_0 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;

                    dx = x[i] - x[j];
                    dy = y[i] - y[j];
                    dz = z[i] - z[j];
                    dx -= L * round(dx / L);
                    dy -= L * round(dy / L);#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define kB 1.38064852e-23						//Bolztmann's constant (J/K)
#define NA 6.02214e23							//Avogadro's constant

#define N_STEPS 10000							//Number of simulation steps
#define BIN 100								//Binning number for radial distribution histogram
#define DIMENSIONS 3							//Number of coordinates (x,y,z)
#define N 256								//Number of particles
#define s 3.405								//Sigma of gas (in Angstroms)
#define p_star 1.1							//Density of gas in (dimentionless)
#define T_star 2.7							//Temperature of gas (dimensionless)
#define dt_star 0.001							//Time step (dimensionless)
#define MASS 39.948							//Mass of gas (amu)
#define epsilon 1.65401693e-21						//Epsilon of gas (J)
#define THERMO (N_STEPS / 2.0)						//Thermostat for specified portion of N_STEPS (i.e. first quarter, half, etc.)
#define APPLY_T 5							//Apply thermostat every X number of steps
static char atom[] = "Ar";						//Atom type

/* THE ABOVE PARAMETERS CAN BE ADJUSTED TO CHANGE THE PROPERTIES OF THE GAS IN THE MACROS BELOW							*
 * Lennard-Jones parameters (epsilon and sigma) taken from: https://www.sciencedirect.com/science/article/pii/002199917590042X?via%3Dihub 	*/

#define p (p_star / pow(s,3.0))						//Density (A^-3)
#define L pow((N/p),1.0/3.0)						//Box length (A)
#define e (epsilon / kB)						//Energy of gas (K)
#define m (MASS * 10.0 / NA / kB)					//Conversion of amu to K*ps^2/A^2
#define T (T_star * e)							//Conversion of temperature to K
#define dt (dt_star * sqrt((m*pow(s,2.0)) / e))				//Conversion of time step to ps

static double r[N][DIMENSIONS];
static double v[N][DIMENSIONS];
static double a[N][DIMENSIONS];

//New array declaration for multicomponent system for sigma and epsilon parameters
//static double sigma[N], epsilon[N];

void crystalLattice();
void wrapToBox();
void initializeVelocities();
double calculateAcceleration();
void velocityVerlet();
double potentialEnergy();
double kineticEnergy();
double generateGaussian();
double meanSquaredVelocity();
void thermostat();
void radialDist(FILE *fp);
void MSD(FILE *fp);
//void VACF(FILE *fp);

int main()
{
	int i,j,k,n;
	int progress;
	double Temp,Press,Pavg,Tavg,V,PE,KE,sqPE,sqKE,PEavg,KEavg,mv,ulrc,plrc,cvPE,cvKE,vSum,v2;

	clock_t start,end;
	double sim_time;

	FILE *ftraj, *fvel, *fener, *fmv, *ftemp, *fpress;
	ftraj = fopen("traj.xyz","w");
	fvel = fopen("velocities.xyz","w");
	fener = fopen("energy.dat","w");
	fmv = fopen("momentum.dat","w");
	ftemp = fopen("temp.dat","w");
	fpress = fopen("pressure.dat","w");

	start = clock();

	crystalLattice();
	calculateAcceleration();
	v2 = meanSquaredVelocity();

	Pavg = 0;
	Tavg = 0;
	PEavg = 0;
	KEavg = 0;
	sqPE = 0;
	sqKE = 0;
	vSum = 0.0;

	progress = floor(N_STEPS / 10);

	for(n=0; n<=N_STEPS; n++)
	{
		wrapToBox();

		if(n == progress)
			printf("[ 10 |");
		else if(n == 2*progress)
			printf(" 20 |");
		else if(n == 3*progress)
			printf(" 30 |");
		else if(n == 4*progress)
			printf(" 40 |");
		else if(n == 5*progress)
			printf(" 50 |");
		else if(n == 6*progress)
			printf(" 60 |");
		else if(n == 7*progress)
			printf(" 70 |");
		else if(n == 8*progress)
			printf(" 80 |");
		else if(n == 9*progress)
			printf(" 90 |");
		else if(n == 10*progress)
			printf(" 100 ]\n");
		fflush(stdout);

		V = calculateAcceleration();

		//Apply thermostat
		if(n != 0 && n % APPLY_T == 0 && n < THERMO)
			thermostat();

		//Momentum at each step
		for(i=0; i<N; i++)
			for(j=0; j<3; j++)
				vSum += v[i][j] / N;
		mv = m * vSum;
		fprintf(fmv,"%d\t %lf\n",n,mv);

		//Write atom position to trajectory file
		for(i=0; i<1; i++)
		{
			fprintf(ftraj,"%d\n\n",N);
			for(j=0; j<N; j++)
			{
				fprintf(ftraj,"%s\t",atom);

				for(k=0; k<3; k++)
					fprintf(ftraj,"%lf\t",r[j][k]);
				fprintf(ftraj,"\n");
			}
		}

		velocityVerlet();

		//Write atom velocities to file
		for(i=0; i<1; i++)
		{
			fprintf(fvel,"%d\n\n",N);
			for(j=0; j<N; j++)
			{
				fprintf(fvel,"%s\t",atom);

				for(k=0; k<3; k++)
					fprintf(fvel,"%lf\t",v[j][k]);
				fprintf(fvel,"\n");
			}
		}

		//Collect simulation parameters after thermostat switches off
		if(n >= THERMO)
		{
			PE = potentialEnergy();
			KE = kineticEnergy();
			fprintf(fener,"%lf\t %lf\t %lf\n",PE,KE,PE+KE);
			sqPE += PE*PE;
			sqKE += KE*KE;
			PEavg += PE;
			KEavg += KE;

			//Temperature from kinetic theory of gases (kB = 1 in reduced units)
			Temp = (2.0/3.0) * KE;
			fprintf(ftemp,"%d\t %lf\n",n,Temp);
			Tavg += Temp;

			//Pressure from virial theorem (kB = 1 in reduced units)
			Press = (p*s*s*s)/N/3.0 * (v2 + V/e);
			fprintf(fpress,"%d\t %lf\n",n,Press);
			Pavg += Press;
		}
	}
	fclose(ftraj);
	fclose(fvel);

	printf("*****************************************************************************\n");
	printf("Calculating Radial Distribution...\n");
	radialDist(ftraj);
	printf("Done!\n");

	printf("*****************************************************************************\n");
	printf("Calculating Mean Squared Displacement...\n");
	MSD(ftraj);
	printf("Done!\n");
	/*
		printf("*****************************************************************************\n");
	        printf("Calculating Velocity Autocorrelation...\n");
	        VACF(fvel);
	        printf("Done!\n");
	*/
	sqPE /= (N_STEPS - THERMO);
	sqKE /= (N_STEPS - THERMO);
	PEavg /= (N_STEPS - THERMO);
	KEavg /= (N_STEPS - THERMO);
	Tavg /= (N_STEPS - THERMO);
	Pavg /= (N_STEPS - THERMO);

	//Specific heat from fluctuation in average potential/kinetic energy and average temp
	cvPE = (3.0/2.0) * pow(1.0 - 2*N*(sqPE-pow(PEavg,2.0))/3.0/pow(Tavg,2.0),-1.0);
	cvKE = (3.0/2.0) * pow(1.0 - 2*N*(sqKE-pow(KEavg,2.0))/3.0/pow(Tavg,2.0),-1.0);

	//Long-range correction to pressure (from Comp. Sim. of Liquids - Allen/Tildesly, pg. 65)
	plrc = (16.0/3.0) * M_PI * pow(p_star,2.0) * ((2.0/3.0) * pow(s/L/2.0,9.0) - pow(s/L/2.0,3.0));

	//Long-range correction to energy (from Comp. Sim. of Liquids - Allen/Tildesly, pg. 65)
	ulrc = (8.0/3.0) * M_PI * N * p_star * e * ((1.0/3.0) * pow(s/L/2.0,9.0) - pow(s/L/2.0,3.0));

	printf("*****************************************************************************\n");
	printf("Momentum: %lf\n",mv);
	printf("*****************************************************************************\n");
	printf("Starting Temperature (K): %lf\n",T);
	printf("Average Temperature (K): %lf\n",Tavg);
	printf("Percent Difference: %.2lf%%\n",(Tavg-T)/T*100.0);
	printf("*****************************************************************************\n");
	printf("--- PRESSURE ---\n");
	printf("Average Reduced Pressure: %lf\n",Pavg);
	printf("Pressure from Long-Range Correction: %lf\n",plrc);
	printf("Average Reduced Pressure with Long-Range Correction: %lf\n",(Pavg+plrc));
	printf("*****************************************************************************\n");
	printf("--- ENERGY ---\n");
	printf("Average Reduced Potential Energy: %lf\n",PEavg/e);
	printf("Energy from Long-Range Correction per Particle: %lf\n",ulrc/e/N);
	printf("Average Reduced Potential Energy with Long-Range Correction: %lf\n",(PEavg+ulrc/N)/e);
	printf("*****************************************************************************\n");
	printf("---SPECIFIC HEAT AT CONSTANT VOLUME---\n");
	printf("Specific Heat from Fluctuations in Potenial Energy: %lf\n",cvPE);
	printf("Specific Heat from Fluctuations in Kinetic Energy: %lf\n",cvKE);
	printf("*****************************************************************************\n");

	end = clock();

	sim_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("SIMULATION RUNTIME: %.2lf seconds\n",sim_time);
	printf("*****************************************************************************\n");

	return(0);
}
/*
// *** UNCOMMENT FOR SIMPLE CUBIC *** //
void crystalLattice()
{
	int i,j,k,n;

	//Number of atoms in each x,y,z direction rounded up to the nearest whole number
	n = ceil(pow(N,1.0/3.0));

	//Atom spacing in given x,y,z direction
	double dr = L / n;

	//Index for total number of atoms, N
	int index = 0;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			for(k=0;k<n;k++)
			{
				r[index][0] = i * dr;
				r[index][1] = j * dr;
				r[index][2] = k * dr;

				index++;
			}

	initializeVelocities();

	//Print initial positions and velocities
//	for(i=0;i<N;i++)
//		printf("%lf\t %lf\t %lf\n",r[i][0],r[i][1],r[i][2]);
//	for(i=0;i<N;i++)
//		printf("%lf\t %lf\t %lf\n",v[i][0],v[i][1],v[i][2]);

}
*/
void crystalLattice()
{
	int i,j,k;
	int index = 0;
	double n,dr,dr2;

	n = pow(N/4.0,1.0/3.0);
	dr = L / n;
	dr2 = dr / 2.0;

	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			for(k=0; k<n; k++)
			{
				if(index < N)
				{
					r[index][0] = i * dr;
					r[index][1] = j * dr;
					r[index][2] = k * dr;
				}
				index++;

				if(index < N)
				{
					r[index][0] = i * dr + dr2;
					r[index][1] = j * dr + dr2;
					r[index][2] = k * dr;
				}
				index++;

				if(index < N)
				{
					r[index][0] = i * dr;
					r[index][1] = j * dr + dr2;
					r[index][2] = k * dr + dr2;
				}
				index++;

				if(index < N)
				{
					r[index][0] = i * dr + dr2;
					r[index][1] = j * dr;
					r[index][2] = k * dr + dr2;
				}
				index++;
			}

	initializeVelocities();

//        for(i=0;i<N;i++)
//        	printf("%lf\t %lf\t %lf\n",r[i][0],r[i][1],r[i][2]);
}

void wrapToBox()
{
	int i,j;

	for (i=0; i<N; i++)
		for (j=0; j<3; j++)
			r[i][j] += -L * floor(r[i][j]/L);
}

void initializeVelocities()
{
	int i,j;
	double vScale;
	double v2 = 0;
	double instant_temp = 0;
	double vCM[3] = {0,0,0};

	//Assign random initial velocities
	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v[i][j] = generateGaussian();

	//Calculating center of mass velocity
	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			vCM[j] += v[i][j];

	for(i=0; i<3; i++)
		vCM[i] /= N;

	//Subtract off center of mass velocity
	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v[i][j] -= vCM[j];

	//Scale initial velocity with scaled temperature (i.e. initial temperature from T* against temperature "right now")
	//instant_temp = m * v^2 / ((3N - 3) * kB) -- kB = 1 in reduced units
	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v2 += v[i][j] * v[i][j] / N;

	for(i=0; i<N; i++)
		instant_temp += m * v2;

	instant_temp /= (3.0*N - 3.0);

	vScale = sqrt(T / instant_temp);

	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v[i][j] *= vScale;
}

double calculateAcceleration()
{
	int i,j,k;
	double F,r2,r6,r12,sor;
	double V = 0;

	//Position of i relative to j
	double rij[3];

	//Initialize acceleration to 0
	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			a[i][j] = 0;

	//Loop over all distinct pairs i,j
	for(i=0; i<N-1; i++)
		for(j=i+1; j<N; j++)
		{
			r2 = 0.0;
			for(k=0; k<3; k++)
			{
				//Component-by-componenent position of i relative to j
				rij[k] = r[i][k] - r[j][k];

				//Periodic boundary conditions
				rij[k] += -L + trunc(rij[k]/L);
				while(rij[k] >= 0.5*L)
					rij[k] -= L;
				while(rij[k] < -0.5*L)
					rij[k] += L;

				//Dot product of the position component
				r2 += rij[k] * rij[k];
			}

			if (r2 < 0.25*L*L)
			{
				sor = (s*s) / r2;
				r6 = sor * sor * sor;
				r12 = r6 * r6;

				F = 48.0 * e/r2 * (r12 - 0.5*r6);

				//Virial sum for pressure calculation
				V += F * r2;

				for(k=0; k<3; k++)
				{
					a[i][k] += rij[k] * F/m;
					a[j][k] -= rij[k] * F/m;
				}
			}
		}

	return V;
}

void velocityVerlet()
{
	int i,j;

	//Calculate acceleration from forces at current position
	calculateAcceleration();

	//Update positions and velocity with current velocity and acceleration
	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
		{
			r[i][j] += v[i][j]*dt + 0.5*a[i][j]*pow(dt,2.0);
			v[i][j] += 0.5*a[i][j]*dt;
		}

	//Update acceleration from current position
	calculateAcceleration();

	//Update velocity with updated acceleration
	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v[i][j] += 0.5*a[i][j]*dt;
}

double potentialEnergy()
{
	int i,j,k;
	double sor,r2,r6,r12;
	double PE = 0.0;
	double rij[3];

	for(i=0; i<N-1; i++)
		for(j=i+1; j<N; j++)
		{
			r2 = 0.0;
			for(k=0; k<3; k++)
			{
				rij[k] = r[i][k] - r[j][k];

				rij[k] += -L + trunc(rij[k]/L);
				while(rij[k] >= L/2.0)
					rij[k] -= L;
				while(rij[k] < -L/2.0)
					rij[k] += L;

				r2 += rij[k] * rij[k];
			}

			if(r2 < 0.25*L*L)
			{
				sor = (s*s) / r2;
				r6 = sor * sor * sor;
				r12 = r6 * r6;

				PE += 4.0 * e * (r12 - r6) / N;
			}
		}

	return PE;
}

double kineticEnergy()
{
	int i,j;
	double v2 = 0.0;
	double KE = 0.0;

	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v2 += v[i][j] * v[i][j];

	KE += 0.5 * m * v2 / N;

	return KE;
}

//Generation of random number sampling via the Marsaglia polar method (shamelessly stolen from Wikipedia: https://en.wikipedia.org/wiki/Marsaglia_polar_method#Implementation)
double generateGaussian()
{
	static bool available = false;
	static double gset;
	double fac,rsq,v1,v2;

	if(!available)
	{
		do
		{
			v1 = 2.0 * rand() / (double)RAND_MAX - 1.0;
			v2 = 2.0 * rand() / (double)RAND_MAX - 1.0;
			rsq = v1*v1 + v2*v2;
		} while(rsq >= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0 * log(rsq)/rsq);
		gset = v1 * fac;
		available = true;

		return v2*fac;
	}

	else
	{
		available = false;
		return gset;
	}
}

double meanSquaredVelocity()
{
	int i;
	double v2;
	double vx2 = 0.0;
	double vy2 = 0.0;
	double vz2 = 0.0;

	for(i=0; i<N; i++)
	{
		vx2 += v[i][0]*v[i][0];
		vy2 += v[i][1]*v[i][1];
		vz2 += v[i][2]*v[i][2];
	}

	v2 = (vx2+vy2+vz2) / N;

	return v2;
}

void thermostat()
{
	int i,j;
	double vScale;
	double v2 = 0.0;
	double instant_temp = 0.0;

	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v2 += v[i][j] * v[i][j] / N;

	for(i=0; i<N; i++)
		instant_temp += m * v2;

	instant_temp /= (3.0*N - 3.0);

	vScale = sqrt(T / instant_temp);

	for(i=0; i<N; i++)
		for(j=0; j<3; j++)
			v[i][j] *= vScale;
}

void radialDist(FILE *fp)
{
	int i,j,n;
	int row;
	double dx,dy,dz,r,dr,rij;
	double g[BIN];
	double x[N],y[N],z[N];

	FILE *fgr, *fxyz;

	dr = L / 2.0 / BIN;

	for(i=0; i<BIN; i++)
		g[i] = 0.0;

	fxyz = fopen("traj.xyz","r");
	for(n=0; n<N_STEPS; n++)
	{
		//Assign the first row in trajectory file as -1
		for(row=-1; row<N; row++)
		{
			if(row == -1)
				fscanf(fxyz,"%*d\n\n"); //Ignore the atom type and empty space beneath
			else
				fscanf(fxyz,"%*s %lf %lf %lf\n",&x[row],&y[row],&z[row]); //Each array contains the coordinates for N atoms at each x,y,z position. For unknown reasons, this method ignores the final block of coordinates, but this is deemed insignificant for calculating g(r)
		}

		if(n >= THERMO)
		{
			for(i=0; i<N-1; i++)
				for(j=i+1; j<N; j++)
				{
					//Apply boundary conditions after thermostating
					dx = x[i] - x[j];
					dx += -L * round(dx/L);
					//while(dx >= L/2.0)
					//	dx -= L;
					//while(dx < -L/2.0)
					//	dx += L;

					dy = y[i] - y[j];
					dy += -L * round(dy/L);
					//while(dy >= L/2.0)
					//        dy -= L;
					//while(dy < -L/2.0)
					//        dy += L;

					dz = z[i] - z[j];
					dz += -L * round(dz/L);
					//while(dz >= L/2.0)
					//        dz -= L;
					//while(dz < -L/2.0)
					//        dz += L;

					rij = sqrt(dx*dx + dy*dy + dz*dz);

					if(rij < L/2.0)
					{
						g[(int)(rij/dr)] += 2.0;
					}
				}
		}
	}
	fclose(fxyz);

	fgr = fopen("radial_dist.dat","w");
	fprintf(fgr,"r\t\t g(r)\n");
	for(i=0; i<BIN; i++)
	{
		r = (i+0.5) * dr;
		g[i] /= THERMO;
		g[i] /= 4.0*M_PI/3.0 * (pow(i+1.0,3.0) - pow(i,3.0)) * pow(dr,3.0) * p;
		fprintf(fgr,"%lf\t%lf\n",r,g[i]/N);
	}
}

void MSD(FILE *fp)
{
	int i,j,n;
	double dx0, dy0, dz0;
	double dx, dy, dz;
	double r2_0, r2_t, msd;
	double x0[N], y0[N], z0[N];
	double x[N], y[N], z[N];

	FILE *fxyz, *fmsd;

	fmsd = fopen("msd.dat", "w");
	fprintf(fmsd, "Time (ps)\t MSD (A^2)\n");

	fxyz = fopen("traj.xyz", "r");
	for (n = 0; n < N_STEPS; n++)
	{
		if (n == THERMO)
		{
			for (int row = -1; row < N; row++)
			{
				if (row == -1)
					fscanf(fxyz, "%*d\n\n");
				else
					fscanf(fxyz, "%*s %lf %lf %lf\n", &x0[row], &y0[row], &z0[row]);
			}
		}

		for (int row = -1; row < N; row++)
		{
			if (row == -1)
				fscanf(fxyz, "%*d\n\n");
			else
				fscanf(fxyz, "%*s %lf %lf %lf\n", &x[row], &y[row], &z[row]);
		}

		if (n >= THERMO)
		{
			double sum_r2 = 0.0;
			int count = 0;
			for (i = 0; i < N - 1; i++)
			{
				for (j = i + 1; j < N; j++)
				{
					dx0 = x0[i] - x0[j];
					dy0 = y0[i] - y0[j];
					dz0 = z0[i] - z0[j];
					dx0 -= L * round(dx0 / L);
					dy0 -= L * round(dy0 / L);
					dz0 -= L * round(dz0 / L);
					r2_0 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;

					dx = x[i] - x[j];
					dy = y[i] - y[j];
					dz = z[i] - z[j];
					dx -= L * round(dx / L);
					dy -= L * round(dy / L);
					dz -= L * round(dz / L);
					r2_t = dx * dx + dy * dy + dz * dz;

					sum_r2 += (r2_t - r2_0) * (r2_t - r2_0);
					count++;
				}
			}

			msd = sum_r2 / count;
			fprintf(fmsd, "%d\t %lf\n", n, msd);
		}
	}

	fclose(fxyz);
	fclose(fmsd);
}
                    dz -= L * round(dz / L);
                    r2_t = dx * dx + dy * dy + dz * dz;

                    sum_r2 += (r2_t - r2_0) * (r2_t - r2_0);
                    count++;
                }
            }

            msd = sum_r2 / count;
            fprintf(fmsd, "%d\t %lf\n", n, msd);
        }
    }

    fclose(fxyz);
    fclose(fmsd);
}
