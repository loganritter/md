#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define kB 1.38064852e-23									//Bolztmann's constant (J/K)
#define NA 6.02214e23										//Avogadro's constant
#define N_STEPS 20000										//Number of simulation steps

#define DIMENSIONS 3										//Number of coordinates (x,y,z)
#define N 27												//Number of particles
#define s 3.405												//Sigma of gas (in Angstroms)
#define p_star 0.6											//Density of gas in (dimentionless)
#define T_star 1.24											//Temperature of gas (dimensionless)
#define dt_star 0.001										//Time step (dimensionless)
#define MASS 39.948											//Mass of gas (amu)
#define epsilon 1.65401693e-21								//Epsilon of gas (J)
static char atom[] = "Ar";									//Atom type

/* THE ABOVE PARAMETERS CAN BE ADJUSTED TO CHANGE THE PROPERTIES OF THE GAS IN THE MACROS BELOW															*
 * Lennard-Jones parameters (epsilon and sigma for Argon) taken from: https://www.sciencedirect.com/science/article/pii/002199917590042X?via%3Dihub 	*/

#define NPB (N / (p_star/pow(s,3.0)))						//Density of gas per box (NPB)
#define L pow(NPB,1.0/3.0)									//Box size (Angstroms)
#define p (p_star / pow(s,3.0))								//Density of gas in A^(-3.0)
#define e (epsilon / kB)									//Energy of gas (K)
#define m (MASS * 10.0 / NA / kB)							//Conversion of amu to K*Ps^2/A^2
#define T (T_star * e)										//Conversion of temperature to K
#define dt (dt_star * sqrt((m*pow(s,2.0)) / e))				//Conversion of time step to s

static double r[N][DIMENSIONS] = {{0.0}};
static double v[N][DIMENSIONS] = {{0.0}};
static double a[N][DIMENSIONS] = {{0.0}};

void crystalLattice();
void wrapToBox();
void initializeVelocities();
double calculateAcceleration();
void velocityVerlet();
double potentialEnergy();
double kineticEnergy();
double generateGaussian();
double meanVelocitySquared();
void thermostat();

int main()
{
	int i,j,k,n;
	int progress;
	double Temp,Press,Pavg,Tavg,V,PE,KE,v2,mv;

	FILE *fpos, *fener, *fopen();
	fpos = fopen("traj.xyz","w");
	fener = fopen("energy.dat","w");

	crystalLattice();
	V = calculateAcceleration();
	v2 = meanVelocitySquared();

	Pavg = 0;
	Tavg = 0;
	mv = 0;

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

		//Apply thermostat for first half of simulation
		if(n != 0)
			if(n % 5 == 0)
				if(n < N_STEPS/2.0)
					thermostat();

		//Write atom position to trajectory file
	    for(i=0;i<1;i++)
        {
			fprintf(fpos,"%d\n\n",N);
            for(j=0;j<N;j++)
            {
				fprintf(fpos,"%s\t",atom);
                for(k=0;k<3;k++)
					fprintf(fpos,"%lf\t",r[j][k]);
                fprintf(fpos,"\n");
            }
        }
		
		velocityVerlet();

		//Collect simulation parameters after thermostat switches off
		if(n >= N_STEPS/2.0)
		{		
			PE = potentialEnergy();
			KE = kineticEnergy();

			fprintf(fener,"%lf\t %lf\t %lf\n",PE,KE,PE+KE);

			//Temperature from kinetic theory of gas (kB = 1 in reduced units) - N is included in KE function
			Temp = (2.0/3.0) * KE;
	
			//Pressure from virial theorem (kB = 1 in reduced units)
			Press = (p/m) * (N*Temp - (1.0/(3.0*Temp)) * V);

			Tavg += Temp;
			Pavg += Press;
		}
	}
	
	Tavg /= N_STEPS / 2.0;
	Pavg /= N_STEPS / 2.0;

	//Momentum per step
	mv += sqrt(v2) / n;

	//Pavg *= e / pow(s,3.0);
	
	printf("Average Momentum: %lf\n",mv);
    printf("Average Temp (K): %lf\n",Tavg);
	printf("Average Pressure (reduced): %lf\n",Pavg);

	return(0);
}

//TO DO: make fcc
void crystalLattice()
{
	int i,j,k,n;

	//Number of atoms in each x,y,z direction rounded up to the nearest whole number
	n = ceil(pow(N,1.0/3.0));

	//Atom spacing in given x,y,z direction
	double dr = L / n;

	//Index for total number of atoms, N
	int atom_tracker = 0;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			for(k=0;k<n;k++)
	       	{
				//Ensures the total number of atoms on the lattice does not exceed N
				if(n < N)
				{
					r[atom_tracker][0] = i * dr;
					r[atom_tracker][1] = j * dr;
					r[atom_tracker][2] = k * dr;
				}
				atom_tracker++;
			}

	initializeVelocities();

//	for(i=0;i<N;i++)
//		printf("%lf\t %lf\t %lf\n",r[i][0],r[i][1],r[i][2]);
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
			v[i][j] = generateGaussian() - 0.5;

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
        for(i=0;i<N;i++)
                for(j=0;j<3;j++)
                        v2 += v[i][j] * v[i][j];

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
	double F,r2,sor;
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
				rij[k] += -L * trunc(rij[k]/L);

				//Cutoff
				while(rij[k] >= L/2.0)
					rij[k] -= L;

				while(rij[k] < -L/2.0)
                                       rij[k] += L;

				//Square of the position component
				r2 += rij[k] * rij[k];
			}
			sor = s / sqrt(r2);
			F = 48.0 * e/sqrt(r2) * (pow(sor,12.0) - 0.5*pow(sor,6.0));

			//Virial sum for pressure calculation
			V += F * sqrt(r2);

			for(k=0;k<3;k++)
		    {
				a[i][k] += rij[k] * F/m;
		        a[j][k] -= rij[k] * F/m;
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
				rij[k] += -L * trunc(rij[k]/L);
	
	            while(rij[k] >= L/2.0)
					rij[k] -= L;
	
				while(rij[k] < -L/2.0)
				   	rij[k] += L;

	            r2 += rij[k] * rij[k];
			}
	
			sor = s / sqrt(r2);
			r6 = pow(sor,6.0);
			r12 = pow(sor,12.0);

			PE += 4.0 * e * (r12 - r6) / N;
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

//Generation of random number sampling via the Marsaglia polar method (https://en.wikipedia.org/wiki/Marsaglia_polar_method)
double generateGaussian()
{
	static bool evaluate = false;
	static double sDeviate;
	double u1,u2,usq,ufac;

 	if(!evaluate)
    {
		do
        {
			u1 = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            u2 = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            usq = u1*u1 + u2*u2;
        } while(usq >= 1.0 || usq == 0.0);

        ufac = sqrt((-2.0 * log(usq))/usq);
        sDeviate = u1 * ufac;
        evaluate = true;

        return u2*ufac;
    }

	else
    {
		evaluate = false;
		return sDeviate;
    }
}

double meanVelocitySquared()
{
	int i;
	double vx2 = 0.0;
	double vy2 = 0.0;
	double vz2 = 0.0;
	double v2;

	for(i=0;i<N;i++)
	{
		vx2 += v[i][0] * v[i][0];
		vy2 += v[i][1] * v[i][1];
		vz2 += v[i][2] * v[i][2];
	}
	v2 = (vx2 + vy2 + vz2) / N;

	return v2;
}

void thermostat(double v2)
{
	int i,j;
        double vScale;
        double instant_temp = 0.0;
        v2 = meanVelocitySquared();

        for(i=0;i<N;i++)
                instant_temp += m * v2;

        instant_temp /= (3.0*N - 3.0);

        vScale = sqrt(T / instant_temp);

        for(i=0;i<N;i++)
                for(j=0;j<3;j++)
                        v[i][j] *= vScale;
}
