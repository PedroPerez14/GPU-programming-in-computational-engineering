//   Shallow Water Equations 2D Simulator (SWE2D)
//  ===========================================================================================
//   Square grid
//  ===========================================================================================
//   Version 1.0 - May 2026
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza   
//  =========================================================================================== 

#include "define.h"
#include "lib/shallow_water.h"
#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
	#include <cblas.h>
	#include <cublas_v2.h>
#endif

// Inputs and storage
char filename[1024];
FILE *fp;
FILE *logFile;
FILE *qFile;

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


// MAIN CODE FUNTION
int main(){


// OPENING DATA STORAGE

	char tempfile[1024];

	//i=system ("rm outputFiles/*.out");	// Remove old files (linux only)

	logFile=fopen("outputFiles/log.out","w");
	qFile=fopen("outputFiles/discharge.out","w");
	
	#if DISPLAY
	Gnuplot gp;			// Graphical output pipe (linux only)
	#endif



// SIMULATOR HEADERS
	printf("\n   SHALLOW WATER EQUATIONS  <<<>>> 2D SIMULATOR");
	printf("\n   Version 1.0 - May 2026");
	printf("\n   Computational Hydraulics Group - University of Zaragoza");
	printf("\n   --------------------------------------------------");

	printf("\n\n>> Case:	");
	printf("\n");

	fprintf(logFile,"\n   SHALLOW WATER EQUATIONS <<<>>>2D SIMULATOR");
	fprintf(logFile,"\n   Version 1.0 - May 2026");
	fprintf(logFile,"\n   Computational Hydraulics Group - University of Zaragoza");
	fprintf(logFile,"\n   --------------------------------------------------");

	fprintf(logFile,"\n\n>> Case:	");
	fprintf(logFile,"\n");

	#ifdef __CUDACC__
		queryAndSetDevice(0);

		cublasHandle_t handle;
		cublasCreate(&handle);
	#endif

// GEOMETRY DATA
	// Load geometry data	
	int nX, nY; // Number of cells in X and Y directions
	int nCells;
	double dx; // Grid spacing
	double Xmin, Ymin; // Domain limits
	double Length, Width; // Domain size

	read_raster_header("input/elevation.input", &nX, &nY, &Xmin, &Ymin, &dx);

	nCells = nX*nY;
	Length = dx*nX; 
	Width = dx*nY;

	printf("\n\n>> Geometry: ");
	printf("\n     nCells %d ", nCells);
	printf("\n     DeltaX %lf - (Xmin,Ymin) = (%lf,%lf) - Width %lf - Length %lf", dx, Xmin, Ymin, Width, Length);			
	printf("\n");

	fprintf(logFile,"\n\n>> Geometry: ");
	fprintf(logFile,"\n     nCells %d ", nCells);
	fprintf(logFile,"\n     DeltaX %lf - (Xmin,Ymin) = (%lf,%lf) - Width %lf - Length %lf", dx, Xmin, Ymin, Width, Length);			
	fprintf(logFile,"\n");


	// Change array declarations to 2D
	double *x, *y; // Cell center coordinates
	double area = dx*dx;

	x = (double*) malloc( (nX+1)*(nY+1)*sizeof(double) );
	y = (double*) malloc( (nX+1)*(nY+1)*sizeof(double) );

	// Initialize grid coordinates
	for(int j = 0; j <= nY; j++) { //rows
		for(int i = 0; i <= nX; i++) { //cols
			x[IDX(i,j,nX+1)] = Xmin + i*dx; // i for x-direction
			y[IDX(i,j,nX+1)] = Ymin + (nY-j)*dx; // j for y-direction
		}
	}

// MAIN VARIABLE DEFINITION

	// Declare variables
	clock_t CPUtime;

	int iter;										// Iteraction index
	int nout;										// Output index
		
	double *h, *qx, *qy;									// Primitive variables
	double *ux, *uy;									// Wave speed
	double *n;										// 2D averaged Manning Coefficient	
	double *zb;										// Bed surface level	
	
	double *DU1;									// Variation conserved variable h in Dt
	double *DU2;									// Variation conserved variable qx in Dt
	double *DU3;									// Variation Conserved variable qy in Dt

	double t, dt;
	double dtSW;

	// Mass error monitor
	double massWt0;									// Initial mass t = t0
	double massWtn;									// Mass at the time t = tn
	double Qin;										// Discharge at intlet for the time t = t0
	double Qout;									// Discharge at outlet for the time t = t0
	double massWerror;
	double Qbalance;

	// Allocate memory 

	h = (double*) malloc( nCells* sizeof(double) );
	qx = (double*) malloc( nCells* sizeof(double) );
	qy = (double*) malloc( nCells* sizeof(double) );
	ux = (double*) malloc( nCells* sizeof(double) );
	uy = (double*) malloc( nCells* sizeof(double) );

	n = (double*) malloc( nCells* sizeof(double) );
	zb = (double*) malloc( nCells* sizeof(double) );

	DU1 = (double*) malloc( nCells* sizeof(double) );
	DU2 = (double*) malloc( nCells* sizeof(double) );
	DU3 = (double*) malloc( nCells* sizeof(double) );


	#ifdef __CUDACC__
	double *d_h, *d_qx, *d_qy;									// Primitive variables
	double *d_ux, *d_uy;									// Wave speed
	double *d_n;										// 2D averaged Manning Coefficient	
	double *d_zb;										// Bed surface level	
	
	double *d_DU1;									// Variation conserved variable h in Dt
	double *d_DU2;									// Variation conserved variable qx in Dt
	double *d_DU3;									// Variation Conserved variable qy in Dt

	cudaMalloc((void**) &d_h, nCells*sizeof(double));
	cudaMalloc((void**) &d_qx, nCells*sizeof(double));
	cudaMalloc((void**) &d_qy, nCells*sizeof(double));
	cudaMalloc((void**) &d_ux, nCells*sizeof(double));
	cudaMalloc((void**) &d_uy, nCells*sizeof(double));

	cudaMalloc((void**) &d_n, nCells*sizeof(double));
	cudaMalloc((void**) &d_zb, nCells*sizeof(double));

	cudaMalloc((void**) &d_DU1, nCells*sizeof(double));
	cudaMalloc((void**) &d_DU2, nCells*sizeof(double));
	cudaMalloc((void**) &d_DU3, nCells*sizeof(double));
	#endif

	//read raster elevation
	read_raster("input/elevation.input", nX, nY, zb);
	
	printf("\n\n>> Geometry loaded");
	printf("\n");
	fprintf(logFile,"\n\n>> Geometry loaded");
	fprintf(logFile,"\n");



// SIMULATION SETUP
	double simTime;
	double Toutput;
	double CFL;
	double n1;

	sprintf(tempfile,"input/simulation.input");
	read_simulation_setup(tempfile, &(simTime), 
		&(Toutput), &(CFL), &(n1));


	double QIN, HIN;
	double HOUT, ZSOUT;
	QIN = -1.0;
	HIN = -1.0;
	HOUT = -1.0;
	ZSOUT = -1.0;

	sprintf(tempfile,"input/simulation.input");
	read_boundary_configuration(tempfile, 
		&(QIN), &(HIN), 
		&(HOUT), &(ZSOUT));	


	printf("\n   --------------------------------------------------");
	printf("\n>> Simulation setup loaded:");
	printf("\n     Simulation time: %.1lf s",simTime);
	printf("\n     CFL: %.2lf",CFL);
	printf("\n     Friction term activated: %d",friction);
	printf("\n     Friction model: %d",Fmodel);
	printf("\n     Number of cells: %d %d Total: %d",nX, nY, nCells);
	printf("\n     dx: %.2lf m", dx);
	printf("\n     Data saved each %d iterations",Niter);
	printf("\n");
	printf("\n     Fluid density: %04.0lf kg/m3",rhow);
	printf("\n     Dynamic viscosity: %.6lf Pa.s",mu);
	printf("\n     Manning coefficient: %.3lf",n1);



	fprintf(logFile,"\n\n   --------------------------------------------------");
	fprintf(logFile,"\n>> Simulation setup loaded:");
	fprintf(logFile,"\n     Simulation time: %.1lf s",simTime);
	fprintf(logFile,"\n     CFL: %.2lf",CFL);
	fprintf(logFile,"\n     Friction term activated: %d",friction);
	fprintf(logFile,"\n     Friction model: %d",Fmodel);
	fprintf(logFile,"\n     Number of cells: %d %d Total: %d",nX, nY, nCells);
	fprintf(logFile,"\n     dx: %.2lf m",dx);
	fprintf(logFile,"\n     Data saved each %d iterations",Niter);
	fprintf(logFile,"\n");
	fprintf(logFile,"\n     Fluid density: %04.0lf kg/m3",rhow);
	fprintf(logFile,"\n     Dynamic viscosity: %.6lf Pa.s",mu);
	fprintf(logFile,"\n     Manning coefficient: %.3lf",n1);
 

	printf("\n\n>> Simulation setup loaded");
	printf("\n");
	fprintf(logFile,"\n\n>> Simulation setup loaded");
	fprintf(logFile,"\n");


// VARIABLE INITIALIZATION
	iter = 1;					// Iteration number
	t = 0.0;					// Initial time

	initialize_variables( nCells, h,  qx,  qy,  ux,  uy, 
		n, DU1,  DU2, DU3);

	printf("\n\n>> Flow initialization completed");
	fprintf(logFile,"\n\n>> Flow initialization completed");



/////////////////////////////////////////////////////////////////////////
// INITIAL CONDITIONS 
/////////////////////////////////////////////////////////////////////////
	
	//read hini. assume qx,qy=0
	read_raster("input/hini.input", nX, nY, h);

	// Flow variables calculation
	compute_initial_flow_variables( nCells,
		h,  qx,  qy,  ux,  uy,  n,  n1);
		
	#ifdef __CUDACC__

	cudaMemcpy(d_h, h, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_qx, qx, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_qy, qy, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ux, ux, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_uy, uy, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_zb, zb, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_n, n, nCells*sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(d_DU1, DU1, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_DU2, DU2, nCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_DU3, DU3, nCells*sizeof(double), cudaMemcpyHostToDevice);

	#endif

	printf("\n\n>> Initial conditions loaded");
	fprintf(logFile,"\n\n>> Initial conditions loaded");

	// Initial mass balance
	#ifdef __CUDACC__
	cublasDasum(handle, nCells, d_h, 1, &massWtn);
	massWtn *= area;
	#else
	compute_water_mass( nCells, h,  area,  &(massWtn));
	#endif
	printf("mass:%lf\n",massWtn);
	Qin = 0.0;	
	Qout = 0.0;
	Qbalance = 0.0;

	fprintf(logFile,"\n>> INITIAL FLOW MASS %.6lf m3",massWtn);	
	fprintf(logFile,"\n");

	// Store initial conditions
	//Cell data
	nout = 0;
	sprintf(filename, "outputFiles/celldata%d.vtk",nout);
	write_vtk_cells(filename, nX, nY, x, y, 
			zb, h, ux, uy);
	nout += 1;


	// Discharge data
	fprintf(qFile,"%.3lf\t %.6lf\t %.6lf\t %.6lf\n",
		t,
		Qin,
		Qout,
		Qbalance);

	printf("\n\nPress INTRO key to start ...");
	printf("\n\n");
	getchar();


CPUtime = clock();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////     TIME LOOP     /////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
while(t <= simTime) {	

// TIME STEP COMPUTATION
	compute_flow_time_step_2D(nX, nY,
			h, ux, uy, dx, 
			&dtSW);
	dt = CFL * dtSW;



// DISPLAY TIME AND TIME STEP
	if(iter%Niter == 0){
		printf("\n--------------------------------------------------------------------------------------------------------");
		printf("\n>> Time: %.3lf seg  -  Iteration %d  -  Time step %.6lf",t,iter,dt);
	}

// FLUXES CALCULATION
	// X-direction fluxes
	#ifdef __CUDACC__
		compute_x_fluxes <<< (nX-1)*nY/nThreads+1 , nThreads >>>
		(nX, nY,
    		d_h, d_qx, d_qy, d_ux, d_uy, d_zb, d_n, 
		d_DU1, d_DU2, d_DU3, 
		dx);
	#else
		compute_x_fluxes(nX, nY,
    		h, qx, qy, ux, uy, zb, n, DU1, DU2, DU3, dx);
	#endif

	// Y-direction fluxes
	#ifdef __CUDACC__
		compute_y_fluxes<<<(nX*(nY-1))/nThreads+1,nThreads>>>(
		nX, nY,
		d_h, d_qx, d_qy, d_ux, d_uy, d_zb, d_n, 
		d_DU1, d_DU2, d_DU3, 
		dx);
	#else
		compute_y_fluxes(nX, nY,
		h, qx, qy, ux, uy, zb, n, DU1, DU2, DU3, dx);
	#endif



// NEGATIVE WATER DEPTH CHECKING
	check_depth_positivity( nCells,
		h,  DU1,  dx,  &(dt));

	  
// SHALLOW WATER CELL UPDATE
#ifdef __CUDACC__
	update_cells_2D <<< nCells/nThreads+1 , nThreads>>>
		(nCells,
		d_h, d_qx, d_qy, d_ux, d_uy,
		d_DU1, d_DU2, d_DU3,
		dx, dt);
#else
	update_cells_2D(nCells,
			h, qx, qy, ux, uy,
			DU1, DU2, DU3,
			dx, dt);
#endif
	
	// Wet/dry conditions
	#ifdef __CUDACC__
		wet_dry_x<<<(nX*(nY-1))/nThreads+1,nThreads>>>(nX, nY,
			d_h, d_qx, d_ux, d_zb);
		wet_dry_y<<<(nX*(nY-1))/nThreads+1,nThreads>>>(nX, nY,
			d_h, d_qy, d_uy, d_zb);
	#else
		wet_dry_x(nX, nY,
		h, qx, ux, zb);

		wet_dry_y(nX, nY,
		h, qy, uy, zb);
	#endif

	// Boundary conditions
	#ifdef __CUDACC__
		set_west_boundary <<< nY/nThreads+1, nThreads>>>
		(nX, nY, d_h, d_qx, d_qy, d_ux, d_uy, QIN, HIN);
		set_east_boundary <<< nY/nThreads+1, nThreads>>>
		(nX, nY, d_h, d_qx, d_qy, d_ux, d_uy, d_zb, QIN, HIN);
		set_north_boundary <<< nX/nThreads+1, nThreads>>>
		(nX, nY, d_h, d_qx, d_qy, d_ux, d_uy);
		set_south_boundary <<< nX/nThreads+1, nThreads>>>
		(nX, nY, d_h, d_qx, d_qy, d_ux, d_uy);
	#else
  		set_west_boundary(nX, nY, h, qx, qy, ux, uy, QIN, HIN);
   		set_east_boundary(nX, nY, h, qx, qy, ux, uy, zb, QIN, HIN);
   		set_north_boundary(nX, nY, h, qx, qy, ux, uy);
   		set_south_boundary(nX, nY, h, qx, qy, ux, uy);
	#endif

// SIMULATION MONITORS
	massWt0 = massWtn;
	#ifdef __CUDACC__
	cublasDasum(handle, nCells, d_h, 1, &massWtn);
	massWtn *= area;
	#else
	compute_water_mass( nCells, h,  area,  &(massWtn));
	#endif

	//printf("mass:%lf\n",massWtn);


	//Compute mass error
	if(massWt0 != 0.0) { 
		massWerror = (massWtn - (massWt0-(Qout-Qin)*dt)) / massWt0; 
		if(fabs(massWerror) < 1e-16) massWerror = 1e-16;
	} else { 
		massWerror = 0.0; 
	}

	Qbalance += (Qin-Qout)*dt;

	
	if(iter%Niter == 0){
		printf("\n\tMass error %.3e\t",massWerror);
		printf("\n\tWater discharge IN %.6lf OUT %.6lf  [m3/s]\n",Qin,Qout);			
	}


	Qin = 0.0; // West inflow
    Qout = 0.0; // East outflow
		

// UPDATE TIME
	iter++;
	t = t+dt;

	if(dt <= 1e-9) { break;}
	
	
	//write files!!!
	
// DATA OUTPUT
	if(t >= nout*Toutput){
		//Cell data
		//output vtk
		sprintf(filename, "outputFiles/celldata%d.vtk",nout);
		write_vtk_cells(filename, nX, nY, x, y, 
			zb, h, ux, uy);
		// Discharge data
		fprintf(qFile,"%.3lf\t %.6lf\t %.6lf\t %.6lf\n",
			t,
			Qin,
			Qout,
			Qbalance);
		
		nout += 1;	
	}	


}
//////////////////////////////////     END TIME LOOP     //////////////////////////////////////////////////////////
sprintf(filename, "outputFiles/celldata%d.vtk",nout);
write_vtk_cells(filename, nX, nY, x, y, 
	zb, h, ux, uy);



CPUtime = clock() - CPUtime ;
double time = ((float)CPUtime)/CLOCKS_PER_SEC;








// DISPLAY FINAL INFORMATION	
printf("\n\n>> Final Time: %.3lf seg",t);
printf("\n\n>> Computation time %.3lf seg",time);
printf("\n\n>> Simulation completed!");
printf("\n\n ");

fprintf(logFile,"\n>> Final Time: %.3lf seg",t);
fprintf(logFile,"\n\n>> Computation time %.3lf seg",time);
fprintf(logFile,"\n\n>> Simulation completed!");
fprintf(logFile,"\n\n ");






// CLOSING DATA STORAGE
fclose(logFile);
fclose(qFile);



} // End of main function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







