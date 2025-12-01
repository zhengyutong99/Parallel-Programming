void hostFE(Molecule *m1, Molecule *m2, double *eVdrTotal, double *eHbondTotal, double *eElecTotal, double *Vtotal, double *TCharge, int *InRange, int *HitKeyAtom, int *HitPosiAtom, int *ChargeContactCount, double (*BOX)[3][2], double *FIT_HYDRO, double *FIT_CHARGE, int *generation, int *MAXGEN);

double calculateEnergyCUDA (Molecule *m1, Molecule *m2) {
	
    double eVdrTotal = 0.0, eHbondTotal = 0.0, eElecTotal = 0.0; 
    double Vtotal = 0.0, TCharge = 0.0;
	double t, len;
	int InRange = 0;
    
    HitKeyAtom=0;
    HitPosiAtom=0;
    ChargeContactCount=0;

	ThresholdContactPairs = GetThresholdContactPairs(m1->count);
	hostFE(m1, m2, &eVdrTotal, &eHbondTotal, &eElecTotal, &Vtotal, &TCharge, &InRange, &HitKeyAtom, &HitPosiAtom, &ChargeContactCount, &BOX, &FIT_HYDRO, &FIT_CHARGE, &generation, &MAXGEN);

	t = ThresholdContactPairs - InRange/m1->count;
	if (t <= 0) {
		RangeRatio = 0;
	} else {
		RangeRatio = t;
	}
  
	if (TCharge !=0 && FIT_CHARGE > 0 ) { 
			t=RangeRatio/ThresholdContactPairs; 
			if (t <= 0) len=0.75;   
			else if (t < 0.05) len=0.5; 
			else if (t >= 0.05 && t < 0.1) len=0.25;
			else len =0.1; 
			TChargeEnergy= 1.5*FIT_CHARGE * len * TCharge; 
   }

	if (!DOCKING) { // constraint for drug screening
			if (FIT_CHARGE < 0 && (FIT_CHARGE+LigChargeAtomCount)>0)
					TChargeEnergy=10.0 * (LigChargeAtomCount); // penalty for charged Ligand

			len=(double)LigHbondCount/(double)(DrugAtomCount);
			if (FIT_HYDRO < 0 && (LigHbondCount >= 6 && (len+FIT_HYDRO)>0))
					Vtotal=Vtotal+ 5.0*LigHbondCount;// penalty for Hydrophlic ligand
			if (KeyAtomCount >= 1 && HitKeyAtom < KeyAtomCount && generation > 40)
					Vtotal=Vtotal+10*(KeyAtomCount-HitKeyAtom);//penalty for violate key-atom constraint
			if (PosiAtomCount >= 1 && generation > 40) {
					len=(double)HitPosiAtom/(double)PosiAtomCount;
					if (len < 0.5) //penalty for violate position-atom constraint
							Vtotal=Vtotal+(PosiAtomCount-HitPosiAtom)*5.0*(0.5-len);
			}
	}

	Vtotal = Vtotal + TChargeEnergy;

	return Vtotal;

}

/*---------------------------------------------------------------
     para[0~2]: axis coordinator
     para[3~5]: axis rotation angles
     para[6~ ]: rotation angles of single bonds
     0,...,0: original X-ray ligand and protein structure
 *---------------------------------------------------------------*/
void ChangeConformation(Chromosome *ind, int idx)
{
    int     i;
       if (generation <= 1) {
           dummy[0].count=mole[0].count ;
           for(i=0 ; i<mole[0].count ;i++)
               dummy[0].atom[i]=mole[0].atom[i];
       }
       else
           for(i=0 ; i<mole[0].count ;i++)
               dummy[0].atom[i].pos=mole[0].atom[i].pos;
       // === transform the rotation angle
      Transformation(&dummy[0],MatrixMultMatrix(
                      UpdateRotate_M(Xaxis,ind->para[3]),
                      UpdateRotate_M(Yaxis,ind->para[4]),
                      UpdateRotate_M(Zaxis,ind->para[5])
                      )) ;
       //=== transform the single bond
       for(i=0 ; i<singleB[0].numSB ;i++) {
           if (idx < UnRotSBNum)
               ind->para[6+i]=0; // single bond keep constand
           else
               PartialTransform(&dummy[0],0,i,ind->para[6+i]) ;
       }
       // === transform to new axis
       Transformation(&dummy[0],UpdateTran_M(ind->para[0],
                               ind->para[1],ind->para[2])) ;
}
/*---------------------------------------------------------------
 * F u n c t i o n E v a l u a t i o n
     Evaluate function values. mole[0]: drug mole[1]:protein
 *---------------------------------------------------------------*/
double FunctionEvaluation(Chromosome *ind, int idx)
{
    int     i;
    double   s,fit ;

    if (idx!=-1) /* -1 is native state */
        ChangeConformation(ind,idx);

    if  (FIT_TYPE=='R') {
        fit=RMSD(&dummy[0]);
        //printf("\n %f %f",fit, BEST_IND.value);getch();
        ind->RangeRatio=0;
        ind->CE=0;
        ind->value=fit;
        ind->singleE[0]=0;
        ind->singleE[1]=0;
    }
    else {
        if (idx ==-1) /* native state */
            // fit=CalculateEnergy(&mole[0],&mole[1]);
            fit=calculateEnergyCUDA(&mole[0],&mole[1]);

        else
            // fit=CalculateEnergy(&dummy[0],&mole[1]);   // inter energy
            fit=calculateEnergyCUDA(&mole[0],&mole[1]);

        ind->RangeRatio=RangeRatio;
        ind->CE=TChargeEnergy;

        // intra energy
        SBE[0]=0;
        if (idx!=-1)  /* native state */
            SBE[0]=0.01*CalculateIntraEnergy(&dummy[0], ind); // 0.01

        fit=fit+SBE[0];
#if 0 // important for correct docking pose
        if (FIT_BOX=='1') { /* important */
         //  fit=fit-RangeRatio*1000;
           if (DrugAtomCount >= 40 || Protein.acidNum >= 100)
            s=0.03;
           else if (DrugAtomCount >= 30)
                s=0.075;
           else if ( DrugAtomCount >= 12)
                s=0.12;
           else
                s=0.15;
           if (RangeRatio < s) // 0.05 penality for ligand out receptor
                fit=fit+(s-RangeRatio)*1000.0; // 1000
        }
#endif
#if 1 // important for correct docking pose
        fit=fit+RangeRatio*20;
#endif
        ind->value=fit;
        ind->singleE[0]=SBE[0] ;
        ind->singleE[1]=SBE[1] ;
    }
    if (idx == -1) { /* native state */
        return(fit);
    }
    if (fit < BEST_IND.value)
    {
        s=0;
        if (DOCKING)
            s = RMSD(&dummy[0]);
        LastGen       = generation;
        LastTime      = time(NULL)-InitTime;
        BEST_IND      = *ind;
        BEST_IND.rmsd = s;
        BEST_IND.idx = idx;
if (FORCERECORD && MAXGEN-generation <=10) {
        for (i=0;i<CurtHBondCount;i++)
           DrugHBond[BEST][i] = DrugHBond[CURT][i];

        for (i=0;i<CurtElectCount;i++)
            DrugElect[BEST][i] = DrugElect[CURT][i];
        BestHBondCount = CurtHBondCount;
        BestElectCount = CurtElectCount;
}

    if (DOCKING ) { // single docking
#ifdef WEB
        if (generation >= 1)
#else
        if (generation > 30 && generation%10 == 0 )
#endif
        {
            for (i=0;i<DrugAtomCount;i++)
                PredDrug3D[i]=dummy[0].atom[i].pos;
        //-------- Modification --------
        //    WritePDBFile();
        }
     } //  if (! DOCKING)
    } // if (fit

if (FORCERECORD && MAXGEN-generation <=10) {
    CurtHBondCount = 0;
    CurtElectCount = 0;
}

    CountE++;
    return(fit);
}

double GetLigBaseEnergy()
{
    int i;
    Chromosome ind;

    for(i=0 ; i<singleB[0].numSB ; i++)
        ind.para[6+i]=0;
    return(CalculateIntraEnergy(&mole[0],&ind));
}

/*--------------------------------------------------------------------
 * O u t s i d e B o x
 *      Check whether an atom of the drug is outside the bounding box of
 * protein. If so, return 1 o/w return 0.
 *--------------------------------------------------------------------*/
double OutsideBox(Point p)
//Point p ;
{
        double   t[3] ,r=0,rate=0.95;
        int     i ;

        t[0]=p.x ; t[1]=p.y ; t[2]=p.z ;
        for(i=0 ; i <3 ; i++) {
                if((t[i]>BOX[i][1]) || (t[i]<BOX[i][0]))/* outside the box */
                {
                        r=r+10000.0;
                       // printf("\ni=%d L=%f U=%f p=%f",i,BOX[i][1],BOX[i][0], t[i]);
                       // getch();
                }
                else if (t[i]>BOX[i][1]*rate)
                     r= r+(t[i]-rate*BOX[i][1])*(t[i]-rate*BOX[i][1]);
                else if (t[i]<BOX[i][0]*rate) /* outside the box */
                     r= r+(t[i]-rate*BOX[i][0])*(t[i]-rate*BOX[i][0]);
        }
        return(r*r) ;             /* inside the box */
}

/*---------------------------------------------------------------------
 * C a l c u l a t e Intra E n e r g y
      Calculate original Lennard-Jones 6-12 potential function.
 *--------------------------------------------------------------------*/
double CalculateIntraEnergy(Molecule *m, Chromosome *ind)
{
   int     i,j,k, AdjF;
   double   len;
   double   Vtotal=0.0;
   bool   HbF;

   for(i=0 ; i<m->count ;i++) {             /* each atom of drug */
        for(j=i+1 ; j<m->count; j++) {      /* each atom of drug */
            AdjF=0;
            for (k=0;k<DrugAtomConnNum[i];k++) // check adjust links
                if (DrugAtomAdj[i][k]==j) {
                    AdjF=1;
                    break;
                }
           if (AdjF) continue;
           len=sqrt(VectorLength(m->atom[i].pos,m->atom[j].pos));
           if(len<6.0) {
                if ( len < 1.45) //2.35 //1.8 // (1.45 is important)
                    Vtotal=Vtotal+10000.0;
                else {
                    if (FIT_TYPE == '5') // 95
                        Vtotal+=ComputeVDW95(&m->atom[i],&m->atom[j],len);
                    else if (FIT_TYPE == '6') // 965
                        Vtotal+=ComputeVDW965(&m->atom[i],&m->atom[j],len,true, &HbF);
                    else if (FIT_TYPE == 'A') // AMBER function
                        Vtotal+=ComputeVDW96(&m->atom[i],&m->atom[j],len);
                }
           } /* end if : len */
       } /* for : j */
   } /* for : i */
    /* calculate torsonal energy */
   if (FIT_SB=='1') {
       for(i=0 ; i<singleB[0].numSB ; i++)
          Vtotal+=ComputeTorEnergy(ind->para[6+i], singleB[0].type[i]);
   }
   return(Vtotal) ;
}


/*===========================================
*  Compute the torsion energy
================================================*/
double ComputeTorEnergy(double dang, char type)
{
    double Ener=0.0, A=0, n=0, ang=0;

    if (type == 1) { // sp3-sp3
        A=3.0; n=3; ang=M_PI;
    }
    else if (type == 2) { // sp3-sp2
        A=1.5; n=6; ang=M_PI;// ang=0.0;
    }
    if (A > 1)
        Ener=A*(1.0+cos(n*dang-ang));

    return(Ener);
}

/*---------------------------------------------------------------------
 * C a l c u l a t e E n e r g y: m1: drug  m2:protein
      Calculate original Lennard-Jones 6-12 potential function.
 *--------------------------------------------------------------------*/
double CalculateEnergy(Molecule *m1, Molecule *m2)
{
   int      i,j ;
   bool     HbF;
   double   len,eVdr,eHbond,eElec,energy,t;
   double   eVdrTotal=0.0,eHbondTotal=0.0,eElecTotal=0.0;
   double   Vtotal=0.0,TCharge=0.0;
   double   InRange=0,RangeTot=0;       /* count the distance between ligand
                                 and recptor within 8 and out 8*/
   HitKeyAtom=0;
   HitPosiAtom=0;
   ChargeContactCount=0;
   ThresholdContactPairs=GetThresholdContactPairs(m1->count);

   for(j=0 ; j<m2->count ;j++) {        /* protein */  
      eVdr=0; eHbond=0; eElec=0;
      for(i=0 ; i<m1->count; i++) {       /* each atom of drug */
         Vtotal+=OutsideBox(m1->atom[i].pos);
         RangeTot++;
         len=sqrt(VectorLength(m1->atom[i].pos,m2->atom[j].pos)) ;
         if(len<=6.0) {          /* original: 8 can't cutoff */
            if (FIT_CHARGE != 0)
                if (m1->atom[i].charge != 0 && m2->atom[j].charge != 0)
                {
                    energy = m2->atom[j].Weight*ComputeElecStatic(m1->atom[i].charge,
                                                      m2->atom[j].charge,len);
                    eElec=eElec+energy;                                  

if (FORCERECORD && MAXGEN-generation <=10) {
                    if (len <= 5.0) {
                    ForceRecord *elect = &DrugElect[CURT][CurtElectCount];

                    elect->ResAcidSeqId = m2->atom[j].AcidSeqId;
                    elect->ResAcidName  = (char*)
                                    ( (m2->atom[j].AcidType == -1) ?
                                      //-------- Modification (Optional) : C++ type conversion for array index --------
                                      m2->atom[j].OtherName : acidName[(int)m2->atom[j].AcidType] );
                    elect->ResAtomSeqId = m2->atom[j].AtomSeqId;
                    elect->ResAtomName  = m2->atom[j].name;
                    elect->DrgAtomSeqId = m1->atom[i].AtomSeqId;
                    elect->DrgAtomName  = m1->atom[i].name;
                    elect->Energy       = energy;

                    CurtElectCount++;
                    }
                   // printf("\n AA-ID=%d AA-Na=%s", elect->ResAcidSeqId,elect->ResAcidName);
}
                    ChargeContactCount++;
                }
             InRange++;
             energy=ComputeVDW965(&m1->atom[i],&m2->atom[j],len,false, &HbF);
             if (HbF) {
                eHbond=eHbond+energy; 
                //printf("\n after %d",HbF);getchar();
             }
             else
                eVdr=eVdr+energy;   
         } /* end if : len */
      } /* for : j: drug*/
      eVdrTotal += eVdr;
      eHbondTotal += eHbond;
      eElecTotal += eElec;
        // printf("\n [serial] each protein atom eVdr: %f", eVdr); 
      Vtotal=Vtotal+eVdr+eHbond;
      TCharge = TCharge+eElec;
      if (GenInterProfile){ //do not run docking: generate protein-ligand profile for statistic only
            //-------- Modification --------
            //old:fprintf(dlogfp,"%11.4f %11.4f %11.4f",eVdr,eHbond,eElec);
            fprintf(dlogfp,"\t%.4f\t%.4f\t%.4f",eVdr,eHbond,eElec);
            oss_profile << tab << eVdr << tab << eHbond << tab << eElec;
            //if (eVdr || eHbond ) printf("a");
      }      
    } /* for : i: Protein */

 //  RangeRatio=InRange/RangeTot;
 /*
   if (TCharge !=0 && FIT_CHARGE > 0 ) { // calculate the electroatcic energy
        len=RangeRatio-0.05; // important for correct docking pose
        if (len < 0) len=0.1;   // 0.1
        else if (len < 0.04) len=0.25; // 0.25
        else if (len >= 0.04 && len < 0.01) len=0.5;
        else len =0.75; // 0.4 (0.75 is important for 3dfr, 1hvr ...)
        TChargeEnergy=2.0 * FIT_CHARGE * len * TCharge; // electrostatic energy
   }
*/
    t=ThresholdContactPairs-InRange/m1->count;
    if (t <= 0) RangeRatio=0;
    else RangeRatio=t;
    if (TCharge !=0 && FIT_CHARGE > 0 ) { // calculate the electroatcic energy
        t=RangeRatio/ThresholdContactPairs; // important for correct docking pose
        if (t <= 0) len=0.75;   // 0.1
        else if (t < 0.05) len=0.5; // 0.25
        else if (t >= 0.05 && t < 0.1) len=0.25;
        else len =0.1; // 0.4 (0.75 is important for 3dfr, 1hvr ...)
        TChargeEnergy= 1.5*FIT_CHARGE * len * TCharge; // electrostatic energy
   }

   if (!DOCKING) { // constraint for drug screening
        if (FIT_CHARGE < 0 && (FIT_CHARGE+LigChargeAtomCount)>0)
            TChargeEnergy=10.0 * (LigChargeAtomCount); // penalty for charged Ligand

        len=(double)LigHbondCount/(double)(DrugAtomCount);
        if (FIT_HYDRO < 0 && (LigHbondCount >= 6 && (len+FIT_HYDRO)>0))
            Vtotal=Vtotal+ 5.0*LigHbondCount;// penalty for Hydrophlic ligand

        if (KeyAtomCount >= 1 && HitKeyAtom < KeyAtomCount && generation > 40)
            Vtotal=Vtotal+10*(KeyAtomCount-HitKeyAtom);//penalty for violate key-atom constraint
        if (PosiAtomCount >= 1 && generation > 40) {
            len=(double)HitPosiAtom/(double)PosiAtomCount;
            if (len < 0.5) //penalty for violate position-atom constraint
                Vtotal=Vtotal+(PosiAtomCount-HitPosiAtom)*5.0*(0.5-len);
        }
   }
   Vtotal=Vtotal+TChargeEnergy;
   if (GenInterProfile) {//do not run docking: generate protein-ligand profile for statistic only
        fprintf(dlogfp,"\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                eVdrTotal,eHbondTotal,eElecTotal,InRange/m1->count,Vtotal);
        oss_profile << endl;
        oss_fitness << tab << eVdrTotal+eHbondTotal+eElecTotal<< tab << eVdrTotal << tab << eHbondTotal << tab << eElecTotal << tab << InRange/m1->count << endl;;
   }

   return(Vtotal) ;
}

/* ==============================
AUTODOCK: Journal of Computer-Aided Molecular Design 10, 1996,
            pp. 293-204
atom1: atom of drug in mole[0]
atom2: atom of receptor in mole[1]
r: distance between atom1 atom2
=========================================*/
double ComputeElecStatic(double q1, double q2, double r)
{
    double fit;
    if (r < 0.5) r=0.5;
    fit=332.0*((q1*q2)/(4.0*r*r));
    //fit=332.0*((q1*q2)/(4.0*r));
    if (fit < -10) fit= -10;
    else if (fit > 10) fit=10;
    return(0.2*fit);
}

/* ==============================
AUTODOCK: Journal of Computer-Aided Molecular Design 10, 1996,
            pp. 293-204
atom1: atom of drug in mole[0]
atom2: atom of receptor in mole[1]
r: distance between atom1 atom2
AMBER-base function
=========================================*/
double ComputeVDW96(Atoms *atom1, Atoms *atom2, double r)
{
    int type1, type2, E1, E2,HF=0;
    double temp, temp1,c1=0,c2=0;
    double Ener=0;

    type1=atom1->HbType;
    type2=atom2->HbType;
    E1=atom1->EngType;
    E2=atom2->EngType;

    if (type1 > 0 && type2 > 0) // check H-bond (1: both, 2: acceptor, 3: donor)
       if (type1 ==1 || type2 ==1 || (type1 ==2 && type2 !=2) ||
          (type1 ==3 && type2 !=3)) {
           HF=1;
           if (E1==3 || E2==3) { // Atom S constructing H-bond
              c1=2.5; c2=1.0; //c1=2.5
           }
           else {
              c1=1.9; c2=5.0; //c1=1.9 , c2=5.0
           }
       }

    if  (HF) { // h-bond
        if (r>0.5) {
            temp=(c1/r);
            temp1=5.0*Power(temp,12); //^r3
            Ener=6.0*Power(temp, 10); //2.0
            Ener=c2*(temp1-Ener);     // increasing factor 10.0
        }
        else
            Ener=100000.0;
    }
    else {
        c1= 0.5*(EngR[E1]+EngR[E2]); // r0
        c2= sqrt(EngE[E1]*EngE[E2]); // e
        if (r>0.5) {
            temp=(c1/r);
            temp1=temp*temp*temp; //^r3
            temp=temp1*temp1; //^6
            Ener=2.0*temp; //2.0
            temp=temp*temp;// ^12
            Ener=c2*(temp-Ener);
        }
        else
            Ener=100000.0;
    }
    if (Ener > 30) {
        Ener=30;
    }
    return(Ener);
}

/* ===============================
Combine 96 and 95
atom1: Ligand atom2:receptor
======================================== */
#define EngBR   0.2 //0.2
#define EngCR   0.9 // 0.9  // 1.0   // 0.5
#define EngDR   6.0 // 7.0  // 8.0  //0.6
double ComputeVDW965(Atoms *atom1, Atoms *atom2, double r, bool intra, bool *HbF)
{
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

    if(!Hbond || FIT_HYDRO <= -9.9)
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
        if (atom2->WeiType == 'K') HitKeyAtom=HitKeyAtom+1;
        if (atom2->WeiType == 'P') HitPosiAtom=HitPosiAtom+1;
    }    

    if (atom2->WeiType == ' ' || atom2->WeiType == 'A' || (atom2->WeiType == 'H' && Hbond)
        || (atom2->WeiType == 'S' && Hbond) || (atom2->WeiType == 'V' && !Hbond))
    {
        Ener=atom2->Weight*Ener;
        if( Hbond && Ener != 0 && r>= 1.95 && r<=3.4 )
        {
           if (FORCERECORD && MAXGEN-generation <=10) {
           ForceRecord *hb = &DrugHBond[CURT][CurtHBondCount];

           hb->ResAcidSeqId = (intra || atom2->AcidSeqId == -1 ) ?
                              99999 : atom2->AcidSeqId;
           hb->ResAcidName  = (char*)
                              (intra ? "PRE" : (atom2->AcidType == -1) ?
                               //Modification
                               //atom2->OtherName : acidName[atom2->AcidType]);
                               atom2->OtherName : acidName[(int)atom2->AcidType]);
           hb->ResAtomSeqId = atom2->AtomSeqId;
           hb->ResAtomName  = atom2->name;
           hb->ResHbondType = (atom2->HbType == 1) ? 'B' :
                              (atom2->HbType == 2) ? 'A' : 'D' ;
           hb->DrgAtomSeqId = atom1->AtomSeqId;
           hb->DrgAtomName  = atom1->name;
           hb->DrgHbondType = (atom1->HbType == 1) ? 'B' :
                              (atom1->HbType == 2) ? 'A' : 'D' ;
           hb->Energy       = Ener;

           HBindingCount++;
           CurtHBondCount++;
         }
        }
        else
             NonHBindingCount++;
   }
   return Ener;
}

/* ==============================
atom1: atom number of drug in mole[0]
atom2: atom number of receptor in mole[1]
r: distance between atom1 atom2
=========================================*/
double ComputeVDW(Atoms *atom1, Atoms *atom2, double r)
{
    int type1, type2, E1, E2;
    double temp, temp1,c1=0,c2=0;
    double Ener=0;

    type1=atom1->HbType;
    type2=atom2->HbType;
    E1=atom1->EngType;
    E2=atom2->EngType;

    if (type1 > 0 && type2 > 0) { // check H-bond (1: both, 2: acceptor, 3: donor)
       if (type1 ==1 || type2 ==1 || (type1 ==2 && type2 !=2) ||
          (type1 ==3 && type2 !=3))  {
          c1=2.9;
          c2=3.0;
          //printf("1");
        }
    }

    if ( c1==0) {
        c1= 0.5*(EngR[E1]+EngR[E2]); // r0
        c2= sqrt(EngE[E1]*EngE[E2]); // e
        //printf("2");
    }
    if (r>1) {
         temp=(c1/r);
         temp1=temp*temp*temp; //^r3
         temp=temp1*temp1; //^6
         Ener=2.0*temp; //2.0
         temp=temp*temp;// ^12
         Ener=c2*(temp-Ener);
    }
    else
         Ener=100000000.0;
    if (Ener > 20) Ener=20;
    return(Ener);
}


/* ===============================
Chemistry and Biology, 1995, vol2. no.5
===================================== */
double a[2][6]={{3.4, 3.6, 4.5, 5.5, -0.4, 20},
                {2.3, 2.6, 3.1, 3.4, -2.0, 20}};
double m[2][3]={{-5.882352941,    -2,  0.4},
                {-8.695652174,    -6.666666667, 6.666666667}
                        };      /* value of slope */
double ComputeVDW95(Atoms *atom1, Atoms *atom2, double r)
{
    int type1, type2;
    double temp;
    int    Hbond=0;
    double Ener=0;

    type1=atom1->HbType;
    type2=atom2->HbType;
    if (type1 > 0 && type2 > 0)  // check H-bond (1: both, 2: acceptor, 3: donor)
       if (type1 ==1 || type2 ==1 || (type1 ==2 && type2 !=2) ||
          (type1 ==3 && type2 !=3))
              Hbond=1;

     if (r <= a[Hbond][0]) {
          Ener= a[Hbond][5]+r*(m[Hbond][0]);
     }
     else if ( r > a[Hbond][0] && r <= a[Hbond][1]) {
          temp= r-a[Hbond][0];
          Ener= temp*(m[Hbond][1]);
     }
      else if ( r > a[Hbond][1] && r <= a[Hbond][2]) {
          Ener=a[Hbond][4];
     }
     else if ( r > a[Hbond][2] && r <= a[Hbond][3]) {
          temp= r-a[Hbond][2];
          Ener= a[Hbond][4]+temp*(m[Hbond][2]);
     }
     else if ( r > a[Hbond][3]) {
          Ener=0;
     }
     if (Hbond ==1) HBindingCount++;
     else NonHBindingCount++;
    return(Ener);
}

double RMSD(Molecule *m)
{
    int i;
    double rmsd=0;
    Point p;
    for (i=0;i<DrugAtomCount;i++) {
        p.x=m->atom[i].pos.x+ProteinCenter[0];
        p.y=m->atom[i].pos.y+ProteinCenter[1];
        p.z=m->atom[i].pos.z+ProteinCenter[2];
        if (RefCrossDrug) /* reference drug is the docking drug */
            rmsd=rmsd+VectorLength(p,DrugOrgPosi[i]);
        else /* reference drug is the ligand in the docking protein */
            rmsd=rmsd+VectorLength(p,RefDrugOrgPosi[i]);
    }
    rmsd=sqrt(rmsd/DrugAtomCount);
    return(rmsd);
}

double GetThresholdContactPairs(int LNum)
{
       double rv=0;                             
       if (LNum <=12) rv=18.97;
       else if (LNum > 12 && LNum <=20) rv=18.46;
       else if (LNum > 20 && LNum <=30) rv=14.15;
       else if (LNum > 30 && LNum <=40) rv=13.26;
       else if (LNum > 40 && LNum <=55) rv=12.9;
       else rv=8.36;
       return(rv);
}       
