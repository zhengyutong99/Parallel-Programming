#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include "fdock2.h"
#include <cuda_runtime.h>

#define BLOCK_SIZE 16

__device__ double d_BOX[3][2]; //read_only
__device__ double d_FIT_HYDRO; // read_only
__device__ double d_FIT_CHARGE; // read_only
__device__ int d_MAXGEN; // read_only
__device__ int d_generation; // read_only
__device__ int d_HitKeyAtom;   
__device__ int d_HitPosiAtom; 
__device__ int d_ChargeContactCount;   

__device__ double outsideBox(Point p) {
	double   t[3] ,r=0,rate=0.95;
	int     i ;
	t[0]=p.x ; t[1]=p.y ; t[2]=p.z ;
	for(i=0 ; i <3 ; i++) {
					if((t[i]>d_BOX[i][1]) || (t[i]<d_BOX[i][0])) /* outside the box */
					{
									r=r+10000.0;
					}
					else if (t[i]>d_BOX[i][1]*rate)
								r= r+(t[i]-rate*d_BOX[i][1])*(t[i]-rate*d_BOX[i][1]);
					else if (t[i]<d_BOX[i][0]*rate) /* outside the box */
								r= r+(t[i]-rate*d_BOX[i][0])*(t[i]-rate*d_BOX[i][0]);
	}
	return r*r; /* inside the box */
}

__device__ double vectorLength(Point p1, Point p2) {
	double   length ;
	length=(p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)
					+(p1.z-p2.z)*(p1.z-p2.z) ;
	return(length) ;
}

__device__ double computeElecStatic(double q1, double q2, double r)
{
    double fit;
    if (r < 0.5) r=0.5;
    fit=332.0*((q1*q2)/(4.0*r*r));
    if (fit < -10) fit= -10;
    else if (fit > 10) fit=10;
    return(0.2*fit);
}

#define EngBR   0.2 //0.2
#define EngCR   0.9 // 0.9  // 1.0   // 0.5
#define EngDR   6.0 // 7.0  // 8.0  //0.6
__device__ double computeVDW965(Atoms *atom1, Atoms *atom2, double r, bool intra, bool *HbF) {
    int type1, type2;
    double temp;
    double Ener  = 0;
    bool   Hbond = false;
    double A,B,C,D,E, m1,m2, m3;

    type1 = atom1->HbType;
    type2 = atom2->HbType;
    *HbF=false;
    if (type1 > 0 && type2 > 0)
    {  // check H-bond (1: both, 2: acceptor, 3: donor)
        if ( type1 ==1 || type2 ==1 ||
             (type1 ==2 && type2 !=2) || (type1 ==3 && type2 !=3) )
        {
            if (atom2->WeiType != 'S') {//     
                A = 2.3;
                B = 2.6;
                C = 3.1; // 3.4 //3.1
                D = 3.6; // 4.5 //3.4
                E = 2.5; // 3.0 //2.0
            } else {    // strict hbond constraint
                A = 1.9; 
                B = 2.0; // 2.0 
                C = 2.2; // 2.3
                D = 2.8; // 2.8
                E = 3.0; // 3.0
            }    
            Hbond = true;
            *HbF=true;
        }
    }

    if(!Hbond || d_FIT_HYDRO <= -9.9)
    { // FIT_HYDRO <= -9.9 (consider only van der wrass Force)
        A= 3.3; // 3.4 3.6 4.5 5.5 -0.4    20
        B= 3.6; // 3.6
        C= 4.5; // 4.5
        D= 6.0; // 5.5
        E= 0.3; // e
        Hbond = false;
    }

    m1 = -20/A;
    m2 = -E/(B-A);
    m3 = E/(D-C);

    if (r <= A ) {
        Ener= 20 + r*m1;
    }
    else if ( r > A && r <= B ) {
        Ener= (r-A)*m2;
    }
    else if ( r > B && r <= C) {
        Ener=-E;
    }
    else if ( r > C && r <= D) {
        temp= r-C;
        Ener= -E+temp*m3;
    }
    else if ( r > D) {
        Ener=0;
    }

    if (Hbond) { // prefer N==O Hbond
        if ((atom1->name=='N' && atom2->name == 'O') ||
            (atom1->name=='O' && atom2->name == 'N') || 
            (atom1->name == 'M' || atom2->name == 'M'))
            Ener=1.4*Ener;
    }
    if (r <= 4.5) {
        if (atom2->WeiType == 'K') {
					// d_HitKeyAtom = d_HitKeyAtom+1;
					atomicAdd(&d_HitKeyAtom, 1);
				}

        if (atom2->WeiType == 'P') {
					// d_HitPosiAtom = d_HitPosiAtom+1;
					atomicAdd(&d_HitPosiAtom, 1);
					}
    }    

    if (atom2->WeiType == ' ' || atom2->WeiType == 'A' || (atom2->WeiType == 'H' && Hbond) || (atom2->WeiType == 'S' && Hbond) || (atom2->WeiType == 'V' && !Hbond)) {
        Ener=atom2->Weight*Ener;
   	}
   return Ener;

}

__global__ void calculateEnergyKernel(Atoms *m1_atoms, Atoms *m2_atoms, int m1_count, int m2_count, double *d_eVdrTotal, double *d_eHbondTotal, double *d_eElecTotal, double *d_Vtotal, double *d_TCharge, int *d_InRange) {
	
	// int row = blockIdy.y * blockDim.y + threadIdx.y;
	int i;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	bool HbF;
	double Vtotal=0.0, TCharge=0.0;
	double len=0.0, energy=0.0;
	double eVdr=0.0, eHbond=0.0, eElec=0.0;
	double eVdrTotal=0.0, eHbondTotal=0.0, eElecTotal=0.0;
	int InRange=0;

	eVdr=0; eHbond=0; eElec=0;
	if (j < m2_count) {
		for (i = 0; i < m1_count; i++) {
			Vtotal += outsideBox(m1_atoms[i].pos);
			len = sqrt(vectorLength(m1_atoms[i].pos, m2_atoms[j].pos));
			if (len <= 6.0) {
				if (d_FIT_CHARGE != 0) {
					if (m1_atoms[i].charge != 0 && m2_atoms[j].charge != 0) {
						energy = m2_atoms[j].Weight * computeElecStatic(m1_atoms[i].charge, m2_atoms[j].charge, len);
						eElec = eElec + energy;
						// d_ChargeContactCount++;
						atomicAdd(&d_ChargeContactCount, 1);
					}
				}
				InRange++;
				energy=computeVDW965(&m1_atoms[i], &m2_atoms[j], len, false, &HbF);
				if (HbF) {
					eHbond = eHbond + energy; 
				} else {
					eVdr = eVdr + energy;  
				}
			}
		}

			// printf("\n [serial] each protein atom on device eVdr: %f", eVdr); 

			eVdrTotal += eVdr;
			eHbondTotal += eHbond;
			eElecTotal += eElec;
			Vtotal = Vtotal + eVdr + eHbond;
			TCharge = TCharge + eElec;

			d_eVdrTotal[j] = eVdrTotal;
			d_eHbondTotal[j] = eHbondTotal;
			d_eElecTotal[j] = eElecTotal;
			d_Vtotal[j] = Vtotal;
			d_TCharge[j] = TCharge;
			d_InRange[j] = InRange;
	}
}

void hostFE(Molecule *m1, Molecule *m2, double *eVdrTotal, double *eHbondTotal, double *eElecTotal, double *Vtotal, double *TCharge, int *InRange, int *HitKeyAtom, int *HitPosiAtom, int *ChargeContactCount, double (*BOX)[3][2], double *FIT_HYDRO, double *FIT_CHARGE, int *generation, int *MAXGEN) {
	double *d_eVdrTotal, *d_eHbondTotal, *d_eElecTotal, *d_Vtotal, *d_TCharge;
	int *d_InRange;
	double h_eVdrTotal_sum = 0, h_eHbondTotal_sum = 0, h_eElecTotal_sum = 0, h_Vtotal_sum = 0, h_TCharge_sum = 0;
	int h_InRange_sum = 0;
	Atoms *d_m1_atoms, *d_m2_atoms;

	// int nDevices;
	// cudaGetDeviceCount(&nDevices);
	// for (int i = 0; i < nDevices; i++) {
	// 		cudaDeviceProp prop;
	// 		cudaGetDeviceProperties(&prop, i);
	// 		std::cout << "Device Number: " << i << std::endl;
	// 		std::cout << "  Device name: " << prop.name << std::endl;
	// 		std::cout << "  Compute capability: " << prop.major << "." << prop.minor << std::endl;
	// }

	cudaMemcpyToSymbol(d_BOX, &BOX, sizeof(BOX));
	cudaMemcpyToSymbol(d_FIT_HYDRO, &FIT_HYDRO, sizeof(FIT_HYDRO));
	cudaMemcpyToSymbol(d_FIT_CHARGE, &FIT_CHARGE, sizeof(FIT_CHARGE));
	cudaMemcpyToSymbol(d_MAXGEN, &MAXGEN, sizeof(MAXGEN));
	cudaMemcpyToSymbol(d_generation, &generation, sizeof(generation));
	cudaMemcpyToSymbol(d_HitKeyAtom, &HitKeyAtom, sizeof(HitKeyAtom));
	cudaMemcpyToSymbol(d_HitPosiAtom, &HitPosiAtom, sizeof(HitPosiAtom));
	cudaMemcpyToSymbol(d_ChargeContactCount, &ChargeContactCount, sizeof(ChargeContactCount));
	
	cudaMalloc(&d_m1_atoms, m1->count * sizeof(Atoms));
	cudaMalloc(&d_m2_atoms, m2->count * sizeof(Atoms));
	
	cudaMemcpy(d_m1_atoms, m1->atom, m1->count * sizeof(Atoms), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m2_atoms, m2->atom, m2->count * sizeof(Atoms), cudaMemcpyHostToDevice);

	// dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE);
	// dim3 gridSize((m1->count + (BLOCK_SIZE - 1)) / BLOCK_SIZE, (m2->count + (BLOCK_SIZE - 1)) / BLOCK_SIZE);
	int num_block = (m2->count + BLOCK_SIZE - 1) / BLOCK_SIZE;

	double *h_eVdrTotal = (double*)malloc(sizeof(double) * m2->count);
	double *h_eHbondTotal = (double*)malloc(sizeof(double) * m2->count);
	double *h_eElecTotal = (double*)malloc(sizeof(double) * m2->count);
	double *h_Vtotal = (double*)malloc(sizeof(double) * m2->count);
	double *h_TCharge = (double*)malloc(sizeof(double) * m2->count);
	int *h_InRange = (int*)malloc(sizeof(int) * m2->count);

	cudaHostRegister(h_eVdrTotal, sizeof(double) * m2->count, cudaHostRegisterDefault);
	cudaHostRegister(h_eHbondTotal, sizeof(double) * m2->count, cudaHostRegisterDefault);
	cudaHostRegister(h_eElecTotal, sizeof(double) * m2->count, cudaHostRegisterDefault);
	cudaHostRegister(h_Vtotal, sizeof(double) * m2->count, cudaHostRegisterDefault);
	cudaHostRegister(h_TCharge, sizeof(double) * m2->count, cudaHostRegisterDefault);
	cudaHostRegister(h_InRange, sizeof(int) * m2->count, cudaHostRegisterDefault);

	cudaHostGetDevicePointer(&d_eVdrTotal, h_eVdrTotal, 0);
	cudaHostGetDevicePointer(&d_eHbondTotal, h_eHbondTotal, 0);
	cudaHostGetDevicePointer(&d_eElecTotal, h_eElecTotal, 0);
	cudaHostGetDevicePointer(&d_Vtotal, h_Vtotal, 0);
	cudaHostGetDevicePointer(&d_TCharge, h_TCharge, 0);
	cudaHostGetDevicePointer(&d_InRange, h_InRange, 0);

	calculateEnergyKernel<<<num_block, BLOCK_SIZE>>>(d_m1_atoms, d_m2_atoms, m1->count, m2->count, d_eVdrTotal, d_eHbondTotal, d_eElecTotal, d_Vtotal, d_TCharge, d_InRange);
	cudaDeviceSynchronize();

  for (int i = 0; i < m2->count; i++){
		  // printf("\n each protein atom on device eVdr: %f", *(d_eVdrTotal + i)); 
			// printf("\n each protein atom on host eVdr: %f", *(h_eVdrTotal + i));
      h_eVdrTotal_sum = h_eVdrTotal_sum + *(h_eVdrTotal + i);
			h_eHbondTotal_sum = h_eHbondTotal_sum + *(h_eHbondTotal + i);
			h_eElecTotal_sum = h_eElecTotal_sum + *(h_eElecTotal + i);
			h_Vtotal_sum = h_Vtotal_sum + *(h_Vtotal + i);
			h_TCharge_sum = h_TCharge_sum + *(h_TCharge + i);
			h_InRange_sum = h_InRange_sum + *(h_InRange + i);
  }

	*eVdrTotal = h_eVdrTotal_sum;
	*eHbondTotal = h_eHbondTotal_sum;
	*eElecTotal = h_eElecTotal_sum;
	*Vtotal = h_Vtotal_sum;
	*TCharge = h_TCharge_sum;
	*InRange = h_InRange_sum;

	cudaMemcpyFromSymbol(&HitKeyAtom, d_HitKeyAtom, sizeof(HitKeyAtom));
	cudaMemcpyFromSymbol(&HitPosiAtom, d_HitPosiAtom, sizeof(HitPosiAtom));
	cudaMemcpyFromSymbol(&ChargeContactCount, d_ChargeContactCount, sizeof(ChargeContactCount));

	cudaHostUnregister(h_eVdrTotal);
	cudaHostUnregister(h_eHbondTotal);
	cudaHostUnregister(h_eElecTotal);
	cudaHostUnregister(h_Vtotal);
	cudaHostUnregister(h_TCharge);
	cudaHostUnregister(h_InRange);

	free(h_eVdrTotal);
	free(h_eHbondTotal);
	free(h_eElecTotal);
	free(h_Vtotal);
	free(h_TCharge);
	free(h_InRange);

}


