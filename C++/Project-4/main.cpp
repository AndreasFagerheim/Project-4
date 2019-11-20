

/* Program for simulating a 2D lattice
 * with the Ising model by using the
 * Metropolis algorithm
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "lib.h"
#include "time.h"
#include <iomanip>

using namespace std;


ofstream ofile;
// ensuring boundary conditions by
// "coupling togheter" the elements
// on the edges
inline int periodic(int i, int dimension, int add) {
    return (i+dimension+add) % (dimension);
}

void Metropolis(int **spin_matrix, int n_spins, double &energy, double &magnetization, double *w, long &idum, int &accepted, int **neighbour);
void SaveResultsFromSeveralMCCycles(string FileName, double *HeatCapacity, double *Susceptibility, double *MeanMagnetization, double *MC_cycles, int N, double *, double *, double *, double *AverageEnergy);
void Initialize(int **, double, double & , double & , long &, bool, int **);
void ProbabilityForGivenEnergy(int **spin_matrix, int n_spins, double &energy, int **neighbour);
void ExpectationValues(double *values, double normalization, double &accepted_total, double &HeatCapacity, double &Susceptibility, double &MeanMagnetization, int L, double &AverageEnergy, double temperature);

int main(){
    //------------------------------------------------------
    // Main parameters
    // No. of MC-values we calculate for
    int MC_cycles = 100000;


    double start_temperature = 2.4;
    int L = 20; // Size of lattice
    double T = 2.4;

    //sets initial state: true = random state, false = uniform state.
    bool random = false;
    // oppgave b = "data_A.txt"
    string FileName = "configs2.txt";
    //------------------------------------------------------
    double values[5]; // Array for keeping our results
    double energy=0, magnetization=0;
    double w[17]; // Array containing all possible energies for the system kind of sth what
    long idum  = -1; // MC-seed
    cout<<"here";
    // Creating array which contain information about the neighbouring
    // spins when we use periodic boundary conditions.
    int **neighbour;
    neighbour = (int **) matrix(L,2,sizeof(L));
    for(int i=0;i<L;i++){
        neighbour[i][0] = periodic(i,L,-1);
        neighbour[i][1] = periodic(i,L,1);
    }


    // Initialising w
    for(int delta_E=-8;delta_E<=8;delta_E++) w[delta_E+8] = exp(-delta_E/start_temperature);
    //for(int i=0;i<MC_cycles;i++) accepted_total[i]=0, probability[i]=0;

    // Setting expectation values to zero
    for(int i=0;i<5;i++) values[i] = 0;
    clock_t start, finish;
    int**spin_matrix;
    // Allocating memory for state matrix and then filling it
    spin_matrix = (int **) matrix(L, L, sizeof(L));
    Initialize(spin_matrix, L, energy, magnetization, idum, random, neighbour);
    // loop over different number of MC-cycles
    // Using metropolis algorithm in order to reach an eiuilibrium

    ofile.open(FileName, ios::out);
    ofile << "# Initial state: " << "random" << endl;
    ofile << "# Temperature: " << T << endl;
    ofile << "# Lattice size L = " << L << endl;
    ofile << "# MCcycles  -  <E>/particle  -  <|M|>/particle   -  N_accepted" << endl;
    //start = clock();
    cout<<"here";
    double E_av = 0,M_av = 0, norm = MC_cycles;
    int confi_accept; // accepted energy states
    confi_accept = 0;
    for(int i=1;i<MC_cycles;i++){

        Metropolis(spin_matrix, L, energy, magnetization, w, idum, confi_accept, neighbour);
        // update expectation values
        values[0] += energy;    values[1] += energy*energy;
        values[2] += magnetization;    values[3] += magnetization*magnetization;  values[4] += fabs(magnetization);

        double norm = 1/((double) i);
        E_av = values[0]*norm;
        M_av = values[4]*norm;
        //cout<<confi_accept<<endl;
        //cout<<"-------------------"<<endl;
        ofile << setw(15) << setprecision(8) << i;
        ofile << setw(15) << setprecision(8) << E_av/L/L;     // <E>
        ofile << setw(15) << setprecision(8) << M_av/L/L;  // <|M|>
        ofile << setw(15) << setprecision(8) << confi_accept*norm<<endl;
        //outfile << setw(15) << setprecision(8) << calc_total_energy(S, L);
        //outfile << setw(15) << setprecision(8) << accepted/((double) cycles) << endl;
         //cout << MC_cycles_no[i] <<" Energy: "<<values[0]<<endl;


        // Calculate probabilities for having a certain state after reaching eiuilibrium
        //if( MC_cycles_no[i] > 100000 ) ProbabilityForGivenEnergy(spin_matrix, L, probability[i], neighbour);
        //finish = clock();
        //time_used[i] = ((double) (finish-start)/CLOCKS_PER_SEC);

    }
    cout<<"finished";
    free_matrix(((void **) spin_matrix)); // free memory

    ofile.close();
   // SaveResultsFromSeveralMCCycles(FileName, HeatCapacity,Susceptibility, MeanMagnetization, MC_cycles_no, N, time_used, accepted_total, probability, AverageEnergy);
    return 0;
}
void Metropolis(int **spin_matrix, int n_spins, double &energy, double &magnetization, double *w, long &idum, int &accepted, int **neighbour){
    for(int y=0;y<n_spins;y++){
        for(int x=0;x<n_spins;x++){
            // Picking random site in lattice
            int iy =    (int) (ran2(&idum)*(double)n_spins);
            int ix = (int) (ran2(&idum)*(double)n_spins);
            // calculating energy difference if one spin is flipped
            int delta_E = 2*spin_matrix[iy][ix]*
                    (spin_matrix[iy][neighbour[ix][0]] +
                    spin_matrix[neighbour[iy][0]][ix] +
                    spin_matrix[iy][neighbour[ix][1]] +
                    spin_matrix[neighbour[iy][1]][ix]);
            // comparing random number with energy difference
            if(  ran2(&idum) <= w[delta_E + 8]){//delta_E <= 0 ||
                accepted += 1;
                //cout<<accepted;
                spin_matrix[iy][ix] *= -1;
                energy += (double) delta_E;
                magnetization += (double) 2*spin_matrix[iy][ix];
            }
        }
    }
}

void Initialize(int **spin_matrix, double n_spins, double &energy, double &magnetization, long &idum, bool random, int **neighbour){
    // can have to different states: uniform or random

    if(random){
        for (int i = 0; i < n_spins; i++) {
            for (int j = 0; j < n_spins; j++) {
                spin_matrix[i][j] = round(ran2(&idum))*2 - 1;
                magnetization += (double) spin_matrix[i][j];
            }
        }
    }

    if(random != true){
        // Filling state matrix with spin one states
            for(int i=0;i<n_spins;i++){
                for(int j=0;j<n_spins;j++){
                    spin_matrix[i][j] = 1.;
                    magnetization += (double) spin_matrix[i][j];
                }
            }
    }


    // Calculating energy and magnetization of initial state
    for(int y=0;y<n_spins;y++){
        for(int x=0;x<n_spins;x++){
            energy -=  (double) spin_matrix[y][x]*(spin_matrix[neighbour[y][0]][x] +
                    spin_matrix[y][neighbour[x][0]]);
            magnetization += (double) spin_matrix[y][x];
        }
    }
}

void SaveResultsFromSeveralMCCycles(string FileName, double *HeatCapacity, double *Susceptibility, double *MeanMagnetization, double *MC_cycles,
                                    int N, double *time_used, double *accepted_total,
                                    double *probability, double *AverageEnergy){
    ofstream ofile;
    ofile.open(FileName);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) <<  "MC";
    ofile << setw(15) << setprecision(8) << "avg E";
    ofile << setw(15) << setprecision(8) << "C_v"; //heatcapacity
    ofile << setw(15) << setprecision(8) << "mean_m";//susceptibility;
    ofile << setw(15) << setprecision(8) << "X";
    ofile << setw(15) << setprecision(8) << "Time";
    ofile << setw(15) << setprecision(8) << "accepted states";
    ofile << setw(15) << setprecision(8) << "prob"<<endl;

    for (int i=0; i<N; i++){
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) <<  MC_cycles[i];
    ofile << setw(15) << setprecision(8) << AverageEnergy[i];
    ofile << setw(15) << setprecision(8) << HeatCapacity[i];
    ofile << setw(15) << setprecision(8) << MeanMagnetization[i];
    ofile << setw(15) << setprecision(8) << Susceptibility[i];
    ofile << setw(15) << setprecision(8) << time_used[i];
    ofile << setw(15) << setprecision(8) << accepted_total[i];
    ofile << setw(15) << setprecision(8) << probability[i]<< endl;
    }

    ofile.close();
}

void ProbabilityForGivenEnergy(int **spin_matrix, int n_spins, double &energy, int **neighbour){
    // Function for calculating the probability of a given energy
    energy = 0;
    for(int y=0;y<n_spins;y++){
        for(int x=0;x<n_spins;x++){
            energy -= (double) spin_matrix[y][x]*(spin_matrix[neighbour[y][0]][x] +
                                                    spin_matrix[y][neighbour[x][0]]);
        }
    }
}

void ExpectationValues(double *values, double normalization, double &accepted_total, double &HeatCapacity, double &Susceptibility, double &MeanMagnetization, int L, double &AverageEnergy, double temperature){
    // Every value is calculated per spin
    accepted_total = accepted_total/normalization;
    HeatCapacity = (values[1]/normalization - values[0]*values[0]/(normalization*normalization))/(L*L*temperature*temperature);
    Susceptibility = (values[3]/normalization - values[4]*values[4]/(normalization*normalization))/(L*L*temperature);
    MeanMagnetization = values[4]/(normalization*L*L);
    AverageEnergy = values[0]/(normalization*L*L);
}




