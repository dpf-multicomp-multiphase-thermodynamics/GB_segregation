#include "./General.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <unistd.h>
#include <map>
using namespace std;

std::vector<double> nearest_neighbor1(std::vector<double> coord) {
    std::vector<double> grid(2);
    
    for (int i = 0; i < 2; i++) {
            grid[i] = std::round(coord[i] * 1000) / 1000.0;
    }
    return grid;
}

std::vector<double> nearest_neighbor2(std::vector<double> coord) {
    std::vector<double> grid(2);
    for (int i = 0; i < 2; i++) {
        grid[i] = std::round(coord[i] * 1000) / 1000.0;
    }

    std::vector<double> dif(2);
    for (int i = 0; i < 2; i++) {
        dif[i] = coord[i] - grid[i];
    }

    std::vector<double> nearest(2);
    for (int i = 0; i < 2; i++) {
        if (dif[i] > 0) {
            nearest[i] = 0.001;
        }
        if (dif[i] < 0) {
            nearest[i] = -0.001;
        }
        if (dif[i] = 0) {
            nearest[i] = 0.0;
        }
    }

    std::vector<double> neighbor(2);
    for (int i = 0; i < 2; i++) {
            neighbor[i] = std::round((coord[i] + nearest[i]) * 1000) / 1000.0;
    }
    return neighbor;
}

double linear_interpolation(std::vector<double> coord1, std::vector<double> coord2, std::vector<double> coord21, std::vector<double> coord12, std::vector<double> coordx, std::vector<double> values) {
    double interpolated_value = 0;
    double threshold = 1e-5;

    double diff1020 = fabs(coord1[0]-coord2[0]);
    double diff1121 = fabs(coord1[1]-coord2[1]);
    double diff11x1 = fabs(coord1[1]-coordx[1]);
    double diff10x0 = fabs(coord1[0]-coordx[0]);
    double diff20x0 = (coord2[0]-coordx[0]);
    double diffx010 = (coordx[0]-coord1[0]);
    double diff2010 = (coord2[0]-coord1[0]);
    double diff21x1 = (coord2[1]-coordx[1]);
    double diffx111 = (coordx[1]-coord1[1]);
    double diff2111 = (coord2[1]-coord1[1]);



    if ((diff1020<threshold) && (diff1121<threshold)){
        interpolated_value = values[0];
    }

    if ((diff11x1<threshold) && (diff10x0>threshold)) {
        interpolated_value = (values[0]*diff20x0 + values[2]*diffx010)/(diff2010);
    }
 
    if ((diff10x0<threshold) && (diff11x1>threshold)) {
        interpolated_value = (values[0]*diff21x1 + values[2]*diffx111)/(diff2111);
    }

    if ((diff10x0>threshold) && (diff11x1>threshold)) {

    double denom = diff1020*diff1121;
    double sec1 = diff10x0*diff11x1*values[2]/denom;
    double sec2 = fabs(diff20x0)*fabs(diffx111)*values[3]/denom;
    double sec3 = fabs(diff20x0)*fabs(diff21x1)*values[0]/denom;
    double sec4 = (diff10x0)*fabs(diff21x1)*values[1]/denom;
    interpolated_value = sec1+sec2+sec3+sec4;
    }

    return interpolated_value;
}

// Function to find the maximum value in a double array
double findMax(const std::vector<double>& tempX, int size) {
    double max = tempX[0]; // Assume the first element is the maximum
    for (int i = 1; i < size; i++) {
        if (tempX[i] > max) {
            max = tempX[i]; // Update max if a larger value is found
        }
    }
    return max;
}

// Function to find the minimum value in a double array
double findMin(const std::vector<double>& tempX, int size) {
    double min = tempX[0]; // Assume the first element is the minimum
    for (int i = 1; i < size; i++) {
        if (tempX[i] < min) {
            min = tempX[i]; // Update min if a smaller value is found
        }
    }
    return min;
}


int main(int argc, char* argv[])
{
//******************************************** Constants
	double Pi = 3.1415;
	double R  = 8.315;

//******************************************** Folders
//******************************************** Parall:OpenMP
	//cout << endl << "#############################################" << endl << endl;
	//cout << "<<<<<<<<<<<<<<<<The number of processors is	" << omp_get_num_procs ( ) << endl;
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(10); // Use 4 threads for all consecutive parallel regions

	int id = 0;
	# pragma omp parallel
	{
	  id = omp_get_thread_num ( );
	  if ( id == 0 )
	  {
	    cout << "<<<<<<<<<<<<<<<<<<<<The maximum number of threads in use is " << omp_get_num_threads ( ) << endl << endl; 
	  }
	}

	//cout << "#############################################" << endl << endl;

// Input parameters
	string filename = "/home/jmurugan/NewPotentials_FCC723K/TCFE10_MOBFE5_FeNiCr_FCC_A1_723K_N1000.dat";
    ifstream file(filename);
    int count = 0;
    char c;
    while (file.get(c)) {
        if (c == '\n') {
            count++;
        }
    }
    file.close();
	int NUM_ROWS = count-1;

	Array3D tcdata(1, 21, NUM_ROWS); // TC DATA: (0)Free energy [J/mol]	(1)Chemical Potential of x1 [J/mol] (2)Chemical potential of x2 [J/mol] 
				     //          (3)Enthalpy of mixing BCC_A2 [J/mol]	(4)Derivative of mixing enthalpy of x1 [J/mol]	(5)Derivative of mixing enthalpy of x2 [J/mol]
				     //          (6)Mobility of base element	(7)Mobility of x1	(8)Mobility of x2	(9)Composition of x1
					 //			 (10)Composition of x2		(11)Stability Coefficient
	
	std::ifstream inFile("/home/jmurugan/NewPotentials_FCC723K/TCFE10_MOBFE5_FeNiCr_FCC_A1_723K_N1000.dat");
    std::string header;
	std::getline(inFile, header);
	std::string line;

	for(int j = 0; j < NUM_ROWS; j++){


        inFile >> tcdata(0, 0, j);
		inFile >> tcdata(0, 1, j);
		inFile >> tcdata(0, 2, j);
		inFile >> tcdata(0, 3, j);
		inFile >> tcdata(0, 4, j);
		inFile >> tcdata(0, 5, j);
		inFile >> tcdata(0, 6, j);
		inFile >> tcdata(0, 7, j);
		inFile >> tcdata(0, 8, j);
        inFile >> tcdata(0, 9, j);
		inFile >> tcdata(0, 10, j);
		inFile >> tcdata(0, 11, j);
	}
	  inFile.close();


  map<double, int> freq;
  for (int i = 0; i < NUM_ROWS; i++) {
    freq[tcdata(0, 9, i)]++;
  }


  vector<double> elements;
  for (auto it = freq.begin(); it != freq.end(); it++) {
    elements.push_back(it->first);
  }
  sort(elements.begin(), elements.end());

  int loops = 0;
  for (int i = 0; i < elements.size(); i++) {
		loops=loops+1;
  }

Array3D spindata(1, 13, NUM_ROWS);
Array3D spinoid(1, 5, NUM_ROWS);

//////////////////////////////// Input parameters /////////////////////////////////////// 

	//double rho1 = 0.9;
	//double rho2 = 0.8;
    	int resol_comp = 1/0.001;
	//double factor= 0.454;
	double E_PotFe = -39655.3;
	double G0_Fe = tcdata(0,0,0);
	double G0_X1 = tcdata(0,0,NUM_ROWS-1);
	double G0_X2 = tcdata(0,0,resol_comp);
	double E_PotX1 = -41199.1;
	double E_PotX2 = -38015.1;

////////////////////////////////// Calculation for effect of density ///////////////////////////////////
#pragma omp parallel
{
    #pragma omp for
	for(int j = 0; j < NUM_ROWS; j++){
		tcdata(0,12,j) = (E_PotX1 + (G0_X1-E_PotX1)) - (E_PotFe + (G0_Fe-E_PotFe)) + tcdata(0, 4, j) - tcdata(0, 1, j); // Entropy term Ni
		tcdata(0,13,j) = (E_PotX2 + (G0_X2-E_PotX2)) - (E_PotFe + (G0_Fe-E_PotFe)) + tcdata(0, 5, j) - tcdata(0, 2, j); // Entropy term Cr
	}
}

////////////////////////////////////// Segregation Isotherms ////////////////////////////////
    int num_comps = 3;
    double gbwidth = 1e-9;
    int sizex = 401;
    int sizex1 = sizex-1;
    int eta = 190; //(sizex-1)/2;
    double dt = 1.e-4;
    int totaltimestep = 15000000;
    int printstep = 5000000;
    double dx = 0.1*1.e-10;
    double ratecon = 1.e-4;  //rate of relaxation of concentration equation
    double raterho = 1.e-2;  //rate of relaxation of density equation
    double soft = 1.0/1.0; // to change the siffness between equations
    std::vector<double> rho(sizex);
    double rhocentre = 0.9;
    double tol = 1.e-12;
    double Xdot_max = 1e-7;
    double kappa_X1 = 1.e-16;
    double kappa_X2 = 1.6*1.e-17;
    std::vector<double> rhodot(sizex);
    std::vector<double> gradrho(sizex);
    std::vector<double> murho(sizex);
    double kappa_rho = 2*1.e-13;




Array3D muX1(1, sizex, NUM_ROWS);
Array3D muX2(1, sizex, NUM_ROWS);
Array3D X1dot(1, sizex, NUM_ROWS);
Array3D X2dot(1, sizex, NUM_ROWS);
Array3D X1(1, sizex, NUM_ROWS);
Array3D X2(1, sizex, NUM_ROWS);
Array3D XFe(1, sizex, NUM_ROWS);
Array3D gradX1(1, sizex, NUM_ROWS);
Array3D gradX2(1, sizex, NUM_ROWS);

//std::vector<double> tempX1(sizex);
//std::vector<double> tempX2(sizex);
//std::vector<double> temprho(sizex);

std::vector<double> gb_dist(sizex);
rho[0] = 1.;
rho[sizex-1] = 1.;
rho[(sizex)/2] = rhocentre;
gb_dist[0] = -2;
gb_dist[sizex-1] = -1*gb_dist[0];
int countii = 1;
for(int i = 1; i < sizex-1; i++){
    gb_dist[i] = gb_dist[0] + i * abs(2*gb_dist[0]) / (sizex - 1);
    if ((i > eta) && (i<(sizex-eta-1))){
        int ii = countii;
        rho[i] = ((1+rhocentre)/2)+((1-rhocentre)/2)*(cos((Pi*ii/((sizex-1-2*eta)/2.))));
        countii++;
    }
    else rho[i] = 1.;
}

double X_1 = 0.6; // initial X1 composition of the system  
double X_2 = 0.28; // initial X2 composition of the system  

int j = (X_1*1000/2)*(1001+(1001-(X_1*1000-1)))+(X_2*1000);

//double Hmix = 1251.65;

    for(int i = 0; i < sizex; i++){  
        X1(0,i,j) = tcdata(0,9,j); 
        X2(0,i,j) = tcdata(0,10,j); 
        XFe(0,i,j) = 1. - X1(0,i,j) - X2(0,i,j);
    }

//std::vector<double> X1dotdt(sizex);
//std::vector<double> X2dotdt(sizex);
//std::vector<double> rhodotdt(sizex);

    // Begin time loop
    for(int t = 0; t <= totaltimestep; t++){
 
    // Apply BC concecntration and density
    X1(0,0,j) = tcdata(0,9,j);
    X2(0,0,j) = tcdata(0,10,j);
    X1(0,sizex-1,j) = tcdata(0,9,j);
    X2(0,sizex-1,j) = tcdata(0,10,j);
    XFe(0,0,j) = 1. - tcdata(0,9,j) - tcdata(0,10,j);
    XFe(0,sizex-1,j) = 1. - tcdata(0,9,j) - tcdata(0,10,j);
    muX1(0,0,j) = tcdata(0,1,j);
    muX2(0,0,j) = tcdata(0,2,j);
    muX1(0,sizex-1,j) = tcdata(0,1,j);
    muX2(0,sizex-1,j) = tcdata(0,2,j);
    rho[0] = 1.;
    rho[sizex-1] = 1.;

    // calculate the gradient terms & Chemical potential
    #pragma omp parallel
    {    
        #pragma omp for
        for(int i = 1; i < sizex1; i++){  
            //int i = sizex1/2;
                //X1(0,i,j) = 0.30;
                //X2(0,i,j) = 0.213;

            double dcx1p = (X1(0,i+1,j) - X1(0,i,j));
			double dcx1m = (X1(0,i,j) - X1(0,i-1,j));
            gradX1(0,i,j) = ((dcx1p - dcx1m)) / dx / dx;

            double dcx2p = (X2(0,i+1,j) - X2(0,i,j));
			double dcx2m = (X2(0,i,j) - X2(0,i-1,j));
            gradX2(0,i,j) = ((dcx2p - dcx2m)) / dx / dx;

            double drhop = (rho[i+1] - rho[i]);
			double drhom = (rho[i] - rho[i-1]);
            gradrho[i] = ((drhop - drhom)) / dx / dx;

                double interpol_HX1 = 0;
                double interpol_HX2 = 0;
                double interpol_SX1 = 0;
                double interpol_SX2 = 0;
                double interpol_Hmix = 0;

                double epsilon = 1e-5;

                vector<double> targetPoint = {X1(0,i,j), X2(0,i,j)};
                vector<double> coord1 = nearest_neighbor1(targetPoint);
                vector<double> coord2 = nearest_neighbor2(targetPoint);
                vector<double> coord21 = {coord2[0], coord1[1]};
                vector<double> coord12 = {coord1[0], coord2[1]};


                if ((std::fmod(X1(0,i,j), 0.001) == 0) && (std::fmod(X2(0,i,j), 0.001) == 0)){
                
                    coord1 = {X1(0,i,j), X2(0,i,j)};
                    coord2 = {X1(0,i,j), X2(0,i,j)};
                    coord21 = {X1(0,i,j), X2(0,i,j)};
                    coord12 = {X1(0,i,j), X2(0,i,j)};
                }

                int line1 = (coord1[0]*1000/2)*(1001+(1001-(coord1[0]*1000-1)))+(coord1[1]*1000);
                int line2 = (coord2[0]*1000/2)*(1001+(1001-(coord2[0]*1000-1)))+(coord1[1]*1000);
                int line3 = (coord2[0]*1000/2)*(1001+(1001-(coord2[0]*1000-1)))+(coord2[1]*1000);
                int line4 = (coord1[0]*1000/2)*(1001+(1001-(coord1[0]*1000-1)))+(coord2[1]*1000);

                vector<double> HX1 = {tcdata(0, 4, line1), tcdata(0, 4, line2), tcdata(0, 4, line3), tcdata(0, 4, line4)};
                vector<double> HX2 = {tcdata(0, 5, line1), tcdata(0, 5, line2), tcdata(0, 5, line3), tcdata(0, 5, line4)};
                vector<double> SX1 = {tcdata(0, 12, line1), tcdata(0, 12, line2), tcdata(0, 12, line3), tcdata(0, 12, line4)};
                vector<double> SX2 = {tcdata(0, 13, line1), tcdata(0, 13, line2), tcdata(0, 13, line3), tcdata(0, 13, line4)};
                vector<double> Hmix = {tcdata(0, 3, line1), tcdata(0, 3, line2), tcdata(0, 3, line3), tcdata(0, 3, line4)};
                
                interpol_HX1 = linear_interpolation(coord1, coord2, coord21, coord12, targetPoint, HX1);
                interpol_HX2 = linear_interpolation(coord1, coord2, coord21, coord12, targetPoint, HX2);
                interpol_SX1 = linear_interpolation(coord1, coord2, coord21, coord12, targetPoint, SX1);
                interpol_SX2 = linear_interpolation(coord1, coord2, coord21, coord12, targetPoint, SX2);
                interpol_Hmix = linear_interpolation(coord1, coord2, coord21, coord12, targetPoint, Hmix);
                
                muX1(0,i,j) = (rho[i]*rho[i]*E_PotX1 + rho[i]*(G0_X1-E_PotX1)) - (rho[i]*rho[i]*E_PotFe + rho[i]*(G0_Fe-E_PotFe)) + rho[i]*rho[i]*interpol_HX1 - interpol_SX1 - kappa_X1*gradX1(0,i,j); // add gradient term // Chemical potential for Ni @ rho
		        muX2(0,i,j)  = (rho[i]*rho[i]*E_PotX2 + rho[i]*(G0_X2-E_PotX2)) - (rho[i]*rho[i]*E_PotFe + rho[i]*(G0_Fe-E_PotFe)) + rho[i]*rho[i]*interpol_HX2 - interpol_SX2 - kappa_X2*gradX2(0,i,j); // Chemical potential for Cr @ rho
        

                X1dot(0,i,j) = ratecon*((dx*dx)/(6.0*kappa_X1))*(tcdata(0, 1, j) - muX1(0,i,j)); // for Element X1 
                X2dot(0,i,j) = ratecon*((dx*dx)/(6.0*kappa_X2))*(tcdata(0, 2, j) - muX2(0,i,j)); // for Element X2 

                if (abs(X1dot(0,i,j)) < Xdot_max)
                    X1(0,i,j)  +=  (X1dot(0,i,j));
                else
                    X1(0,i,j) += Xdot_max*(X1(0,i,j)/abs(X1(0,i,j)));

                if (abs(X2dot(0,i,j)) < Xdot_max)
                    X2(0,i,j)  +=  (X2dot(0,i,j));
                else
                    X2(0,i,j) += Xdot_max*(X2(0,i,j)/abs(X2(0,i,j)));

            XFe(0,i,j) = 1. - X1(0,i,j) -X2(0,i,j);
            
            murho[i] = X1(0,i,j)*(2*rho[i]*E_PotX1 + G0_X1-E_PotX1) + X2(0,i,j)*(2*rho[i]*E_PotX2 + G0_X2-E_PotX2) - X1(0,i,j)*kappa_rho*gradrho[i] - X2(0,i,j)*kappa_rho*gradrho[i] + 2*rho[i]*interpol_Hmix; // add gradient term // Chemical potential for Ni @ rho
            rhodot[i]= raterho*((dx*dx)/(6.0*kappa_rho))*(0.0 - murho[i]);
            rho[i] +=  soft*rhodot[i];

            if (rho[i] > 1) {rho[i] = 1.;}
            if (rho[i] < 0) {rho[i] = 0;}
			rho[int(sizex/2.)] = rhocentre;
            rho[1] = 1.;
			rho[sizex-2] = 1.;
  }
}        

            if (t % printstep == 0) {
            double timename = static_cast<double>(t); 
            std::ostringstream filenameStream;
            filenameStream << "output_" << std::fixed << std::setprecision(0) << timename << ".dat";
            std::string filename = filenameStream.str();
            std::ofstream outputFile(filename);
            if (outputFile.is_open()) {
                for (int i = 0;  i < sizex; i++){
                    outputFile << gb_dist[i] <<"   "<< rho[i] <<"   "<< murho[i]<<"   "<< XFe(0,i,j) <<"   "<< X1(0,i,j) <<"   "<< X2(0,i,j) <<"   "<< muX1(0,i,j) <<"   "<< muX2(0,i,j) << endl;
                }   
                outputFile.close();
            }
        }
}     
return 0;
}



