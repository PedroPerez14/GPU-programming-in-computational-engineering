//   Shallow Water Equations 2D Simulator (SWE2D)
//  ===========================================================================================
//   Square grid
//  ===========================================================================================
//   Version 1.0 - Mar 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza   
//  =========================================================================================== 


#include "shallow_water.h"

int read_raster_header(const char *filename, int *nX, int *nY, 
	double *Xmin, double *Ymin, double *dx){
	
	double nodata;

	FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return -1;
    }
	
	// Read header
    fscanf(file, "NCOLS %d\n", nX);
    fscanf(file, "NROWS %d\n", nY);
    fscanf(file, "XLLCORNER %lf\n", Xmin);
    fscanf(file, "YLLCORNER %lf\n", Ymin);
    fscanf(file, "CELLSIZE %lf\n", dx);
    fscanf(file, "NODATA_VALUE %lf\n", &nodata);

    fclose(file);
    return 0;

}

int read_raster(const char *filename, int nX, int nY, 
	double *data){
	
	char buffer[256];  

	FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return -1;
    }

	// Skip header lines  
	for (int i = 0; i < 6; i++) {  
        if (!fgets(buffer, sizeof(buffer), file)) {  
            perror("Error reading header");  
            fclose(file);  
            return -1;  
        }  
    }  
      
	// Read matrix  
	for (int j = 0; j < nY; j++) {     // Rows  
		for (int i = 0; i < nX; i++) { // Columns  
			if (fscanf(file, "%lf", &data[IDX(i, j, nX)]) != 1) {  
				perror("Error reading raster data");  
				fclose(file);  
				return -1;  
			}  
		}  
}  


    fclose(file);
    return 0;

}


int read_simulation_setup(char *tempfile, double *simTime, 
	double *Toutput, double *CFL, double *n1 ){		

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {

		item = "simTime:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*simTime) = (double)strtod(position, &endPtr); //double 
        }

		item = "Toutput:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Toutput) = (double)strtod(position, &endPtr); //double 
        }		

		item = "CFL:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*CFL) = (double)strtod(position, &endPtr); //double 
        }

		item = "n1:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*n1) = (double)strtod(position, &endPtr); //double 
        }									

    }
	fclose(fp);

	return 1;
}



int read_boundary_configuration(char *tempfile, 
	double *QIN, double *HIN, 
	double *HOUT, double *ZSOUT){		

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {

		item = "QIN:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*QIN) = (double)strtod(position, &endPtr); //double 
        }

		item = "HIN:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*HIN) = (double)strtod(position, &endPtr); //double 
        }		

		item = "HOUT:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*HOUT) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZSOUT:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZSOUT) = (double)strtod(position, &endPtr); //double 
        }									

    }
	fclose(fp);

	return 1;
}



void h_initialize_variables(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy, 
	double *n, 
	double *DU1, double *DU2, double *DU3){

    	for(int ic = 0; ic < nCells; ic++){
		h[ic] = 0.0;
		qx[ic] = 0.0;
		qy[ic] = 0.0;
		ux[ic] = 0.0;
		uy[ic] = 0.0;

		n[ic] = 0.0;

		DU1[ic] = 0.0;	
		DU2[ic] = 0.0;
		DU3[ic] = 0.0;

	}

}



int h_compute_initial_flow_variables(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy, 
	double *n, double n1){

    int ic;
	for(ic=0; ic<nCells; ic++){		
		// Minimum depth control
		if(h[ic] >= hmin){
			// Minimum depth control
			ux[ic]=qx[ic]/h[ic];
			uy[ic]=qy[ic]/h[ic];
		}else{
			qx[ic]=0.0;
			qy[ic]=0.0;
			ux[ic]=0.0;
			uy[ic]=0.0;
		}		
		// Cell Manning Coefficient
        n[ic] = n1;

	}			

    return 1;

}   



void h_compute_water_mass(int nCells,
	double *h, double area, double *massWt){

	(*massWt) = 0.0;
	for(int ic=0; ic<nCells; ic++){
		(*massWt) += h[ic]*area;		
	
	}

}

void h_compute_flow_time_step_2D(int nX, int nY,
    double *h, double *ux, double *uy, double dx,
    double *dtSW) {
    
    double uxROE, uyROE, cROE;
    double dt = 1e6;

    // Process horizontal interfaces
	for (int j = 0; j < nY; j++) {     // Rows  
		for (int i = 0; i < nX-1; i++) { // Columns  
            int left = IDX(i,j,nX);
            int right = IDX(i+1,j,nX);
			if (h[left] > tol9 || h[right] > tol9){
				cROE = sqrt(g*0.5*(h[left]+h[right]));
				uxROE = (sqrt(h[left])*ux[left]+sqrt(h[right])*ux[right]) / (sqrt(h[left])+sqrt(h[right]));
				if(fabs(uxROE)+cROE > tol9) {
					dt = MIN(dt, dx/(fabs(uxROE)+cROE));
				}
			}
        }
    }

    // Process vertical interfaces
	for (int j = 0; j < nY-1; j++) {     // Rows  
		for (int i = 0; i < nX; i++) { // Columns  
            int bottom =IDX(i,j,nX);
            int top = IDX(i,j+1,nX);
			if (h[bottom] > tol9 || h[top] > tol9){
				cROE = sqrt(g*0.5*(h[top]+h[bottom]));
				uyROE = (sqrt(h[top])*uy[top]+sqrt(h[bottom])*uy[bottom]) / (sqrt(h[top])+sqrt(h[bottom]));
            	if(fabs(uyROE)+cROE > tol9) {
             		dt = MIN(dt, dx/(fabs(uyROE)+cROE));
   		        }
			}
        }
    }

    (*dtSW) = dt;
}



void h_compute_x_fluxes(int nX, int nY,
	double *h, double *qx, double *qy, double *ux, double *uy, 
	double *zb, double *n,
  	double *DU1, double *DU2, double *DU3,
	double dx){

	double uxROE, uyROE, cROE;                		// ROE averaged variables
	double lambda1, lambda2, lambda3;				// Eigenvalues
	double e1[3], e2[3], e3[3];					// Eigenvectors

	double deltah0, deltaqx0, deltaqy0;				// Spatial increments of the conserved variables
	double alpha1, alpha2, alpha3;					// Wave-strenghts

	double So, Fp, b1P;			// Hydrostatic pressure force at bed
	double nav;								// 2D averaged Manning Coefficient
	double Sf, Ff, b1F;				// Fupwind force at bed
	double beta1, beta2, beta3;					// Source term coefficients

	double hstar,hx,hxx,Qx,Qstart;			// Non-negative area values control - Beta limitation

	double lb1l, lb1r, lb3l, lb3r, lambda1plus, lambda1minus, lambda3plus, lambda3minus;		// Entropy correction

	double hav;				

	double un, normalX, normalY;

	//normal direction in X
	normalX=1.0;
	normalY=0.0;

   	// Process horizontal interfaces
	for (int j = 0; j < nY; j++) {     // Rows  
		for (int i = 0; i < nX-1; i++) { // Columns  
            int left = IDX(i,j,nX);
            int right = IDX(i+1,j,nX);
			if(h[left] >= tol9 || h[right] >= tol9){			// Wet wall

				// Averaged quantities at the walls
				hav = (h[left]+h[right])/2.;

				// ROE averaged variables at the walls
				// Note about the notation:
				// ux=u (x-axis component of the velocity)
				// uy=v (y-axis component of the velocity)
				cROE = sqrt(g*hav);
				uxROE = (sqrt(h[left])*ux[left]+sqrt(h[right])*ux[right]) / (sqrt(h[left])+sqrt(h[right]));	
				uyROE = (sqrt(h[left])*uy[left]+sqrt(h[right])*uy[right]) / (sqrt(h[left])+sqrt(h[right]));	

				//normal projection of the velocity
				un = uxROE*normalX + uyROE*normalY;

				// Eigenvalues 
				lb1l = ux[left]*normalX + uy[left]*normalY-sqrt(g*h[left]);
				lb1r = ux[right]*normalX + uy[right]*normalY-sqrt(g*h[right]);
				lambda1 = un-cROE;
				lambda2 = un;
				lb3l = ux[left]*normalX + uy[left]*normalY+sqrt(g*h[left]);
				lb3r = ux[right]*normalX + uy[right]*normalY+sqrt(g*h[right]);
				lambda3 = un+cROE;

				// Eigenvectors
				e1[0] = 1.;      e1[1] = uxROE-cROE*normalX;      e1[2] = uyROE - cROE*normalY;
				e2[0] = 0.;      e2[1] = -cROE*normalY;           e2[2] = cROE*normalX;
				e3[0] = 1.;      e3[1] = uxROE+cROE*normalX;      e3[2] = uyROE + cROE*normalY;

				// Wave-strenght coefficients
				deltah0 = h[right]-h[left];
				deltaqx0 = qx[right]-qx[left];
				deltaqy0 = qy[right]-qy[left];
				alpha1 = deltah0/2.0 - (1./(2.*cROE))*(deltaqx0*normalX + deltaqy0*normalY - un*deltah0);
				alpha2 = (1./cROE)*((deltaqy0 - uyROE*deltah0)*normalX - (deltaqx0 - uxROE*deltah0)*normalY); 
				alpha3 = deltah0/2.0 + (1./(2.*cROE))*(deltaqx0*normalX + deltaqy0*normalY - un*deltah0);

				// Source term coefficients
				//Bed slope differential formulation

				#if bed_slope_integral
					double deltaz = zb[right] - zb[left];
					double l1 = zb[left] + h[left];
					double l2 = zb[right] + h[right];
					double dzp = deltaz;
					double hp;
					
					if (deltaz >= 0.0) {
						hp = h[left];
						if (l1 < zb[right]) {
							dzp = h[left];
						}
					} else {
						hp = h[right];
						if (l2 < zb[left]) {
							dzp = -h[right];
						}
					}
					
					// Calculate the bed slope source term
					b1P= (1./(2.*cROE)) * g * (hp - 0.5 * fabs(dzp)) * dzp;
				#else
					So = -(zb[right]-zb[left])/dx;
					Fp = g*hav*So*dx;
					b1P = (-1./(2.*cROE)) * Fp;
				#endif

				//friction
				nav = 0.5*(n[left]+n[right]);
				if(h[left] < tol9 || h[right] < tol9){
					nav=0.0;
				}
				if(friction == 1) {
					// Manning
					if(Fmodel == 1) { Sf = nav*nav*un*sqrt(uxROE*uxROE + uyROE*uyROE)/pow(hav,4./3.); } // Unit formulation
					//Viscous
					if(Fmodel == 2) { Sf = 3.*mu*un / (rhow*g*hav); } // Unit formulation
				} else {
					Sf = 0.0;
				}

				Ff = -g*hav*Sf*dx;

				b1F = (-1./(2.*cROE)) * Ff;
				
				// Friction fix
				Qstart = (qx[left]+alpha1*e1[1])*normalX + (qy[left]+alpha1*e1[2])*normalY - b1P;
				Qx = Qstart - b1F;
				if(fabs(Qstart)<tol9) Qstart=0.0;
				if(fabs(Qx)<tol9) Qx=0.0;
				
				if(Qx*Qstart < 0.0){
					b1F = Qstart;
				}

				// beta coefficient
				beta1 = b1P + b1F;
				beta2= 0.0;
				beta3 = -beta1;

				// Positivity fix
				hstar = h[left]+alpha1;		
				if(hstar<tol9) hstar = 0.0;			

				if( (h[left]>0.0 && h[right]>0.0) &&
					(hstar>0.0) &&
					(lambda1*lambda3 < 0.0) ){

					hx = h[left]+alpha1-beta1/lambda1;
					hxx = h[right]-alpha3+beta3/lambda3;
					if(fabs(hx)<tol9) hx = 0.0;
					if(fabs(hxx)<tol9) hxx = 0.0;

					if(hx < 0.0){
						beta1 = hstar*lambda1;
						beta3 = -beta1;
					}

					if(hxx < 0.0){
						beta1 = hstar*lambda3;
						beta3 = -beta1;
					}				
				}
				

				// Update contributions		
				hx = h[left]+alpha1-beta1/lambda1;
				hxx = h[right]-alpha3+beta3/lambda3;
				if(fabs(hx)<tol9) hx = 0.0;
				if(fabs(hxx)<tol9) hxx = 0.0;		

				// First wave
				if(h[left]<tol9 && hx < 0.0){	// dry-wet wall
					DU1[right] += (lambda1*alpha1-beta1)*e1[0];
					DU2[right] += 0.0;
					DU3[right] += 0.0;

				} else if(h[right]<tol9 && hxx < 0.0){ // wet-dry wall
					DU1[left] += (lambda1*alpha1-beta1)*e1[0];
					DU2[left] += 0.0;
					DU3[left] += 0.0;

				} else if(lb1l < 0.0 && lb1r > 0.0){
					lambda1minus = lb1l*(lb1r-lambda1)/(lb1r-lb1l);
					DU1[left] += (lambda1minus*alpha1-beta1)*e1[0];
					DU2[left] += (lambda1minus*alpha1-beta1)*e1[1];
					DU3[left] += (lambda1minus*alpha1-beta1)*e1[2];
					lambda1plus = lb1r*(lambda1-lb1l)/(lb1r-lb1l);
					DU1[right] += (lambda1plus*alpha1)*e1[0];
					DU2[right] += (lambda1plus*alpha1)*e1[1];
					DU3[right] += (lambda1plus*alpha1)*e1[2];

				} else if(lambda1 < 0.0){
					DU1[left] += (lambda1*alpha1-beta1)*e1[0];
					DU2[left] += (lambda1*alpha1-beta1)*e1[1];
					DU3[left] += (lambda1*alpha1-beta1)*e1[2];

				} else if(lambda1 >= 0.0){
					DU1[right] += (lambda1*alpha1-beta1)*e1[0];	
					DU2[right] += (lambda1*alpha1-beta1)*e1[1];
					DU3[right] += (lambda1*alpha1-beta1)*e1[2];
				}
				
				// Second wave
				if(h[left]<tol9 && hx < 0.0){ 		// dry-wet wall
					DU1[right] += (lambda2*alpha2-beta2)*e2[0];
					DU2[right] += 0.0;
					DU3[right] += 0.0;

				} else if(h[right]<tol9 && hxx < 0.0){	// wet-dry wall
					DU1[left] += (lambda2*alpha2-beta2)*e2[0];
					DU2[left] += 0.0;
					DU3[left] += 0.0;

				} else if(lambda2 < 0.0){
					DU1[left] += (lambda2*alpha2-beta2)*e2[0];
					DU2[left] += (lambda2*alpha2-beta2)*e2[1];
					DU3[left] += (lambda2*alpha2-beta2)*e2[2];

				} else if(lambda2 >= 0.0){
					DU1[right] += (lambda2*alpha2-beta2)*e2[0];
					DU2[right] += (lambda2*alpha2-beta2)*e2[1];
					DU3[right] += (lambda2*alpha2-beta2)*e2[2];
				}


				// Third wave
				if(h[left]<tol9 && hx < 0.0){ 		// dry-wet wall
					DU1[right] += (lambda3*alpha3-beta3)*e3[0];
					DU2[right] += 0.0;
					DU3[right] += 0.0;

				} else if(h[right]<tol9 && hxx < 0.0){	// wet-dry wall
					DU1[left] += (lambda3*alpha3-beta3)*e3[0];
					DU2[left] += 0.0;
					DU3[left] += 0.0;

				} else if(lb3l < 0.0 && lb3r > 0.0){
					lambda3minus=lb3l*(lb3r-lambda3)/(lb3r-lb3l);
					DU1[left] += (lambda3minus*alpha3)*e3[0];
					DU2[left] += (lambda3minus*alpha3)*e3[1];
					DU3[left] += (lambda3minus*alpha3)*e3[2];
					lambda3plus=lb3r*(lambda3-lb3l)/(lb3r-lb3l);
					DU1[right] += (lambda3plus*alpha3-beta3)*e3[0];
					DU2[right] += (lambda3plus*alpha3-beta3)*e3[1];
					DU3[right] += (lambda3plus*alpha3-beta3)*e3[2];

				} else if(lambda3 < 0.0){
					DU1[left] += (lambda3*alpha3-beta3)*e3[0];
					DU2[left] += (lambda3*alpha3-beta3)*e3[1];
					DU3[left] += (lambda3*alpha3-beta3)*e3[2];

				} else if(lambda3 >= 0.0){
					DU1[right] += (lambda3*alpha3-beta3)*e3[0];
					DU2[right] += (lambda3*alpha3-beta3)*e3[1];
					DU3[right] += (lambda3*alpha3-beta3)*e3[2];
				}

			} // End wet walls
		}
	} // End of fluxes calculation - Wall loop				

}


void h_compute_y_fluxes(int nX, int nY,
	double *h, double *qx, double *qy, double *ux, double *uy, 
	double *zb, double *n,
  	double *DU1, double *DU2, double *DU3,
	double dx){

	double uxROE, uyROE, cROE;                		// ROE averaged variables
	double lambda1, lambda2, lambda3;				// Eigenvalues
	double e1[3], e2[3], e3[3];					// Eigenvectors

	double deltah0, deltaqx0, deltaqy0;				// Spatial increments of the conserved variables
	double alpha1, alpha2, alpha3;					// Wave-strenghts

	double So, Fp, b1P;			// Hydrostatic pressure force at bed
	double nav;								// 2D averaged Manning Coefficient
	double Sf, Ff, b1F ;				// Fupwind force at bed
	double beta1, beta2, beta3;					// Source term coefficients

	double hstar,hx,hxx,Qx,Qstart;			// Non-negative area values control - Beta limitation

	double lb1l, lb1r, lb3l, lb3r, lambda1plus, lambda1minus, lambda3plus, lambda3minus;		// Entropy correction

	double hav;

	double un, normalX, normalY;

	//normal direction in Y
	normalX=0.0;
	normalY=1.0;
	

    // Process vertical interfaces
	for (int j = 0; j < nY-1; j++) {     // Rows  
		for (int i = 0; i < nX; i++) { // Columns  
			int bottom =IDX(i,j,nX);
			int top = IDX(i,j+1,nX);
			if(h[bottom] >= tol9 || h[top] >= tol9){			// Wet wall

				// Averaged quantities at the walls
				hav = (h[bottom]+h[top])/2.;

				// ROE averaged variables at the walls
				cROE = sqrt(g*hav);
				uxROE = (sqrt(h[bottom])*ux[bottom]+sqrt(h[top])*ux[top]) / (sqrt(h[bottom])+sqrt(h[top]));	
				uyROE = (sqrt(h[bottom])*uy[bottom]+sqrt(h[top])*uy[top]) / (sqrt(h[bottom])+sqrt(h[top]));	

				//normal projection of the velocity
				un = uxROE*normalX + uyROE*normalY;

				// Eigenvalues 
				lb1l = ux[bottom]*normalX + uy[bottom]*normalY-sqrt(g*h[bottom]);
				lb1r = ux[top]*normalX + uy[top]*normalY-sqrt(g*h[top]);
				lambda1 = un-cROE;
				lambda2 = un;
				lb3l = ux[bottom]*normalX + uy[bottom]*normalY+sqrt(g*h[bottom]);
				lb3r = ux[top]*normalX + uy[top]*normalY+sqrt(g*h[top]);
				lambda3 = un+cROE;

				// Eigenvectors
				e1[0] = 1.;      e1[1] = uxROE-cROE*normalX;      e1[2] = uyROE - cROE*normalY;
				e2[0] = 0.;      e2[1] = -cROE*normalY;           e2[2] = cROE*normalX;
				e3[0] = 1.;      e3[1] = uxROE+cROE*normalX;      e3[2] = uyROE + cROE*normalY;

				// Wave-strenght coefficients
				deltah0 = h[top]-h[bottom];
				deltaqx0 = qx[top]-qx[bottom];
				deltaqy0 = qy[top]-qy[bottom];

				alpha1 = deltah0/2.0 - (1./(2.*cROE))*(deltaqx0*normalX + deltaqy0*normalY - un*deltah0);
				alpha2 = (1./cROE)*((deltaqy0 - uyROE*deltah0)*normalX - (deltaqx0 - uxROE*deltah0)*normalY); 
				alpha3 = deltah0/2.0 + (1./(2.*cROE))*(deltaqx0*normalX + deltaqy0*normalY - un*deltah0);

				// Source term coefficients

				#if bed_slope_integral
					double deltaz = zb[top] - zb[bottom];
					double l1 = zb[bottom] + h[bottom];
					double l2 = zb[top] + h[top];
					double dzp = deltaz;
					double hp;
					
					if (deltaz >= 0.0) {
						hp = h[bottom];
						if (l1 < zb[top]) {
							dzp = h[bottom];
						}
					} else {
						hp = h[top];
						if (l2 < zb[bottom]) {
							dzp = -h[top];
						}
					}
					
					// Calculate the bed slope source term
					b1P= (1./(2.*cROE)) * g * (hp - 0.5 * fabs(dzp)) * dzp;
				#else
					So = -(zb[top]-zb[bottom])/dx;
					Fp = g*hav*So*dx;
					b1P = (-1./(2.*cROE)) * Fp;
				#endif


				nav = 0.5*(n[bottom]+n[top]);
				if(h[bottom] < tol9 || h[top] < tol9){
					nav=0.0;
				}
				if(friction == 1) {
					// Manning
					if(Fmodel == 1) { Sf = nav*nav*un*sqrt(uxROE*uxROE + uyROE*uyROE)/pow(hav,4./3.); } // Unit formulation
					//Viscous
					if(Fmodel == 2) { Sf = 3.*mu*un / (rhow*g*hav); } // Unit formulation
				} else {
					Sf = 0.0;
				}

				Ff = -g*hav*Sf*dx;

				b1F = (-1./(2.*cROE)) * Ff;

				// Friction fix
				Qstart = (qx[bottom]+alpha1*e1[1])*normalX + (qy[bottom]+alpha1*e1[2])*normalY - b1P;
				Qx = Qstart - b1F;
				if(fabs(Qstart)<tol9) Qstart=0.0;
				if(fabs(Qx)<tol9) Qx=0.0;

				if(Qx*Qstart < 0.0){
					b1F = Qstart;
				}

				// beta coefficients				
				beta1 = b1P + b1F;
				beta2= 0.0;
				beta3 = -beta1;

				// Positivity fix
				hstar = h[bottom]+alpha1;		
				if(hstar<tol9) hstar = 0.0;			

				if( (h[bottom]>0.0 && h[top]>0.0) &&
					(hstar>0.0) &&
					(lambda1*lambda3 < 0.0) ){

					hx = h[bottom]+alpha1-beta1/lambda1;
					hxx = h[top]-alpha3+beta3/lambda3;
					if(fabs(hx)<tol9) hx = 0.0;
					if(fabs(hxx)<tol9) hxx = 0.0;

					if(hx < 0.0){
						beta1 = hstar*lambda1;
						beta3 = -beta1;
					}

					if(hxx < 0.0){
						beta1 = hstar*lambda3;
						beta3 = -beta1;
					}				
				}

				// Update contributions		
				hx = h[bottom]+alpha1-beta1/lambda1;
				hxx = h[top]-alpha3+beta3/lambda3;
				if(fabs(hx)<tol9) hx = 0.0;
				if(fabs(hxx)<tol9) hxx = 0.0;		

				// First wave
				if(h[bottom]<tol9 && hx < 0.0){	// dry-wet wall
					DU1[top] += (lambda1*alpha1-beta1)*e1[0];
					DU2[top] += 0.0;
					DU3[top] += 0.0;

				} else if(h[top]<tol9 && hxx < 0.0){ // wet-dry wall
					DU1[bottom] += (lambda1*alpha1-beta1)*e1[0];
					DU2[bottom] += 0.0;
					DU3[bottom] += 0.0;

				} else if(lb1l < 0.0 && lb1r > 0.0){
					lambda1minus = lb1l*(lb1r-lambda1)/(lb1r-lb1l);
					DU1[bottom] += (lambda1minus*alpha1-beta1)*e1[0];
					DU2[bottom] += (lambda1minus*alpha1-beta1)*e1[1];
					DU3[bottom] += (lambda1minus*alpha1-beta1)*e1[2];
					lambda1plus = lb1r*(lambda1-lb1l)/(lb1r-lb1l);
					DU1[top] += (lambda1plus*alpha1)*e1[0];
					DU2[top] += (lambda1plus*alpha1)*e1[1];
					DU3[top] += (lambda1plus*alpha1)*e1[2];

				} else if(lambda1 < 0.0){
					DU1[bottom] += (lambda1*alpha1-beta1)*e1[0];
					DU2[bottom] += (lambda1*alpha1-beta1)*e1[1];
					DU3[bottom] += (lambda1*alpha1-beta1)*e1[2];

				} else if(lambda1 >= 0.0){
					DU1[top] += (lambda1*alpha1-beta1)*e1[0];	
					DU2[top] += (lambda1*alpha1-beta1)*e1[1];
					DU3[top] += (lambda1*alpha1-beta1)*e1[2];
				}
				
				// Second wave
				if(h[bottom]<tol9 && hx < 0.0){ 		// dry-wet wall
					DU1[top] += (lambda2*alpha2-beta2)*e2[0];
					DU2[top] += 0.0;
					DU3[top] += 0.0;

				} else if(h[top]<tol9 && hxx < 0.0){	// wet-dry wall
					DU1[bottom] += (lambda2*alpha2-beta2)*e2[0];
					DU2[bottom] += 0.0;
					DU3[bottom] += 0.0;

				} else if(lambda2 < 0.0){
					DU1[bottom] += (lambda2*alpha2-beta2)*e2[0];
					DU2[bottom] += (lambda2*alpha2-beta2)*e2[1];
					DU3[bottom] += (lambda2*alpha2-beta2)*e2[2];

				} else if(lambda2 >= 0.0){
					DU1[top] += (lambda2*alpha2-beta2)*e2[0];
					DU2[top] += (lambda2*alpha2-beta2)*e2[1];
					DU3[top] += (lambda2*alpha2-beta2)*e2[2];
				}


				// Third wave
				if(h[bottom]<tol9 && hx < 0.0){ 		// dry-wet wall
					DU1[top] += (lambda3*alpha3-beta3)*e3[0];
					DU2[top] += 0.0;
					DU3[top] += 0.0;

				} else if(h[top]<tol9 && hxx < 0.0){	// wet-dry wall
					DU1[bottom] += (lambda3*alpha3-beta3)*e3[0];
					DU2[bottom] += 0.0;
					DU3[bottom] += 0.0;

				} else if(lb3l < 0.0 && lb3r > 0.0){
					lambda3minus=lb3l*(lb3r-lambda3)/(lb3r-lb3l);
					DU1[bottom] += (lambda3minus*alpha3)*e3[0];
					DU2[bottom] += (lambda3minus*alpha3)*e3[1];
					DU3[bottom] += (lambda3minus*alpha3)*e3[2];
					lambda3plus=lb3r*(lambda3-lb3l)/(lb3r-lb3l);
					DU1[top] += (lambda3plus*alpha3-beta3)*e3[0];
					DU2[top] += (lambda3plus*alpha3-beta3)*e3[1];
					DU3[top] += (lambda3plus*alpha3-beta3)*e3[2];

				} else if(lambda3 < 0.0){
					DU1[bottom] += (lambda3*alpha3-beta3)*e3[0];
					DU2[bottom] += (lambda3*alpha3-beta3)*e3[1];
					DU3[bottom] += (lambda3*alpha3-beta3)*e3[2];

				} else if(lambda3 >= 0.0){
					DU1[top] += (lambda3*alpha3-beta3)*e3[0];
					DU2[top] += (lambda3*alpha3-beta3)*e3[1];
					DU3[top] += (lambda3*alpha3-beta3)*e3[2];
				}

			} // End wet walls
		}
	} // End of fluxes calculation - Wall loop				

}





void h_check_depth_positivity(int nCells, 
	double *h, double *DU1, double dx, double *dt){
   		for(int ic = 0; ic < nCells; ic++){
	double aux1;

		aux1 = h[ic] - (*dt)*DU1[ic]/dx;		
		if(fabs(aux1) < tol9) aux1 = 0.0;

		while(aux1 < 0.0) {
			(*dt) = 0.50 * (*dt);
			
			aux1 = h[ic] - (*dt)*DU1[ic]/dx;		
			if(fabs(aux1) < tol9) aux1 = 0.0;
		}
	}

} 



void h_update_cells_2D(int nCells, 
	double *h, double *qx, double *qy, double *ux, double *uy, 
	double *DU1, double *DU2, double *DU3,
	double dx, double dt){
	
        for(int ic = 0; ic < nCells; ic++){
		h[ic] = h[ic] - DU1[ic]*dt/dx;
		qx[ic] = qx[ic] - DU2[ic]*dt/dx;
		qy[ic] = qy[ic] - DU3[ic]*dt/dx;

		if(h[ic] >= tol9){
			// Minimum depth control
			if(h[ic] >= hmin){
				ux[ic] = qx[ic]/h[ic];
				uy[ic] = qy[ic]/h[ic];
			}else{
				qx[ic]=0.0;
				qy[ic]=0.0;
				ux[ic] = 0.0;
				uy[ic] = 0.0;
			}				
		} else {
			h[ic] = 0.0;
			qx[ic] = 0.0;
			qy[ic] = 0.0;
			ux[ic] = 0.0;
			uy[ic] = 0.0;
		}

		// Reset U fluxes differences
		DU1[ic]=0.0;
		DU2[ic]=0.0;
		DU3[ic]=0.0;
		
	}			


}  

void h_wet_dry_x(int nX, int nY,
	double *h, double *qx, double *ux,	double *zb){
	for (int j = 0; j < nY; j++) {     // Rows  
		for (int i = 0; i < nX-1; i++) { // Columns  
            int left = IDX(i,j,nX);
            int right = IDX(i+1,j,nX);
			if((h[right] < tol9) && (h[left] + zb[left] <zb[right])){
				qx[left]=0.0;
				ux[left]=0.0;
			}
			if((h[left] < tol9) && (h[right] + zb[right] <zb[left])){
				qx[right]=0.0;
				ux[right]=0.0;
			}
		}
	}
}

void h_wet_dry_y(int nX, int nY,
	double *h, double *qy, double *uy,	double *zb){

	// Process vertical interfaces
	for (int j = 0; j < nY-1; j++) {     // Rows  
		for (int i = 0; i < nX; i++) { // Columns  
			int bottom =IDX(i,j,nX);
			int top = IDX(i,j+1,nX);
			if((h[top] < tol9) && (h[bottom] + zb[bottom] <zb[top])){
				qy[bottom]=0.0;
				uy[bottom]=0.0;
			}
			if((h[bottom] < tol9) && (h[top] + zb[top] <zb[bottom])){
				qy[top]=0.0;
				uy[top]=0.0;
			}
		}
	}
}



void h_set_west_boundary(int nX, int nY, double *h, 
	double *qx, double *qy, double *ux, double *uy, 
	double QIN, double HIN){
	for(int j = 0; j < nY; j++) { //
		int idx = IDX(0,j,nX);
		
		if(QIN > 0.0){
			qx[idx] = QIN;
			qy[idx] = 0.0;
		}
		if(HIN > 0.0) {
			h[idx] = HIN;
			qy[idx] = 0.0;
		}

		if(QIN ==-1 && HIN ==-1){
			qx[idx]=0.0;
		}

		if(h[idx]>=hmin){
			ux[idx] = qx[idx]/h[idx];
			uy[idx] = qy[idx]/h[idx];
		}else{
			ux[idx] = 0.0;
			uy[idx] = 0.0;
		}
	}

}

void h_set_east_boundary(int nX, int nY, double *h, 
	double *qx, double *qy, double *ux, double *uy, double *zb,
	double HOUT, double ZSOUT){
		
	for(int j = 0; j < nY; j++) { //
		int idx = IDX(nX-1,j,nX);

		if(HOUT > 0.0) { 
			h[idx] = HOUT; 
			qy[idx] = 0.0;
		}
	
		if(ZSOUT > 0.0) { 
			h[idx] = ZSOUT-zb[idx]; 
			qy[idx] = 0.0;
		}

		if(HOUT ==-1 && ZSOUT ==-1){
			qx[idx]=0.0;
		}

		if(h[idx]<0.0) h[idx]=0.0;

		if(h[idx]>=hmin){
			ux[idx] = qx[idx]/h[idx];
			uy[idx] = qy[idx]/h[idx];
		}else{
			ux[idx] = 0.0;
			uy[idx] = 0.0;
		}

	}

}


void h_set_north_boundary(int nX, int nY, double *h, 
	double *qx, double *qy, double *ux, double *uy){
	for(int i = 0; i < nX; i++) {
		int idx = IDX(i,nY-1,nX);
		qy[idx] = 0.0;
		uy[idx]= 0.0;
	}

}

void h_set_south_boundary(int nX, int nY, double *h, 
	double *qx, double *qy, double *ux, double *uy){
		
	for(int i = 0; i < nX; i++) {
		int idx = IDX(i,0,nX);
		qy[idx] = 0.0;
		uy[idx]= 0.0;
	}

}



int write_vtk_cells(const char *filename, int nX, int nY, double *x, double *y, 
    double *zb, double *h, double *ux, double *uy) {

    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening output file");
        return -1;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "VTK file containing raster data\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(file, "POINTS %d float\n", (nX + 1) * (nY + 1));

    // Write coordinates
    for (int j = 0; j <= nY; j++) {
        for (int i = 0; i <= nX; i++) {
            fprintf(file, "%.6f %.6f 0.0\n", x[IDX(i,j,nX+1)], y[IDX(i,j,nX+1)]);
        }
    }

    // Define the cells (connectivity)
    fprintf(file, "CELLS %d %d\n", nX*nY, 5*nX*nY);
    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            // For structured grid, points are ordered by their indices
            int p0 = IDX(i, j, nX+1);
            int p1 = IDX(i+1, j, nX+1);
            int p2 = IDX(i+1, j+1, nX+1);
            int p3 = IDX(i, j+1, nX+1);
            fprintf(file, "4 %d %d %d %d\n", p0, p1, p2, p3);
        }
    }

    // Define the cell types
    fprintf(file, "CELL_TYPES %d\n", nX*nY);
    for (int i = 0; i < nX*nY; i++) {
        fprintf(file, "9\n");  // 9 is the VTK code for quad cells
    }

    // Write elevation values as CELL_DATA
    fprintf(file, "CELL_DATA %d\n", nX*nY);
    fprintf(file, "SCALARS z double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f\n", zb[IDX(i,j,nX)]);
        }
    }

    fprintf(file, "SCALARS h double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f\n", h[IDX(i,j,nX)]);
        }
    }
	fprintf(file, "SCALARS h+z double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f\n", h[IDX(i,j,nX)] + zb[IDX(i,j,nX)]);
        }
    }

    fprintf(file, "VECTORS v double\n");
    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f %.6f 0.0\n", ux[IDX(i,j,nX)], uy[IDX(i,j,nX)]);
        }
    }

    fclose(file);
    return 0;
}


