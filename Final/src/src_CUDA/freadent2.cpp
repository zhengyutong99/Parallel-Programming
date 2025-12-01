int GetDrugIdx(int serId)
{
    int i,Tid;
    Tid = serId- DrugATOMS[0].AtomSeqId;
    if (Tid < 0) {
        return (-1);
    }
    //-------- Modification (Critical) -------- 
    if (DrugAtomCount>Tid)
    if (DrugATOMS[Tid].AtomSeqId==serId) return(Tid);

    for (i=0;i<DrugAtomCount;i++)
        if (DrugATOMS[i].AtomSeqId==serId) return(i);
    return(-1);
}

/* ====
Generate the connection of peptide residues
======*/
void GenPeptideConn(int SeqId, int ConnId, char SBType)
{
    int  Id1, Id2,S,C, count, i, Rep1=0, Rep2=0;

    Id1   = GetDrugIdx(SeqId);
    Id2   = GetDrugIdx(ConnId);
    if (Id1 < 0 || Id2 < 0)
        return;     // connections not in the Drug Molecular

    S=DrugATOMS[Id1].numB;
    C=DrugATOMS[Id2].numB;
    for (i=0;i<S;i++)
        if (DrugATOMS[Id1].adj[i] == Id2) Rep1=1;
    if (Rep1 == 0) {
        DrugATOMS[Id1].adj[S]=Id2; DrugATOMS[Id1].numB++;
    }
    for (i=0;i<C;i++)
        if (DrugATOMS[Id2].adj[i] == Id1) Rep2=1;
    if (Rep2 == 0) {
        DrugATOMS[Id2].adj[C]=Id1; DrugATOMS[Id2].numB++;
    }
/*
    printf("\n%d Id=%d count=%d, c1=%d c2=%d c3=%d c4=%d",SeqId, Id1, DrugATOMS[Id1].numB,
       DrugATOMS[Id1].adj[0],DrugATOMS[Id1].adj[1],
       DrugATOMS[Id1].adj[2],DrugATOMS[Id1].adj[3]);
    printf("\n%d Id=%d count=%d, c1=%d c2=%d c3=%d c4=%d",ConnId, Id2, DrugATOMS[Id2].numB,
       DrugATOMS[Id2].adj[0],DrugATOMS[Id2].adj[1],
       DrugATOMS[Id2].adj[2],DrugATOMS[Id2].adj[3]);
    getch();
*/

    if (SBType != 'N' && Rep1 == 0 && Rep2 == 0) {
        count=singleB[0].numSB;
        singleB[0].atom[count][0]=Id1;
        singleB[0].atom[count][1]=Id2;
        singleB[0].type[count]=SBType;
        singleB[0].numSB++;
        /*
        printf("\n C=%3d A1=%3d A2=%3d Ty=%3d",count,singleB[0].atom[count][0],
            singleB[0].atom[count][1], singleB[0].type[count]);
        getch();
        */
    }
}
/*==================================================
    Generate the connection and SB of peptide drug
=========================================================*/
void GenPeptideResidueConn(int idx, int curAcidId)
{
    char NSB='N';
    int S=0,Id1;

    // (N=>CA) (CA=>C) (C=>O)
    S=PeptideDrug[idx].Start;
    if (S==0) { // complete (N=>CA, CA=>C)
        if (idx==0) { // first Residue && no connection
            Id1   = GetDrugIdx(PeptideDrug[idx].SeqId[0]);
            if (DrugATOMS[Id1].numB==0) { // no connection
                DrugATOMS[Id1].charge=0.5;
                LigChargeAtomCount++;
            }
        }
        GenPeptideConn(PeptideDrug[idx].SeqId[0],PeptideDrug[idx].SeqId[1],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[1-S],PeptideDrug[idx].SeqId[2-S],0);
    }
    else if (S==1) { // lack of N (CA=>C))
         GenPeptideConn(PeptideDrug[idx].SeqId[1-S],PeptideDrug[idx].SeqId[2-S],0);
    }
    if (idx==PepResCount-1) { // last residue && no connection
        Id1   = GetDrugIdx(PeptideDrug[idx].SeqId[3-S]);
        printf("\n 11111 %d %c",DrugATOMS[Id1].numB,DrugATOMS[Id1].name);
        if (DrugATOMS[Id1].numB==0) {// no connection
            printf("\n 11111 %d %d %c",Id1,DrugATOMS[Id1].numB,DrugATOMS[Id1].name);
            DrugATOMS[Id1].charge=-0.5;
            LigChargeAtomCount++;
        }
    }
        // C=>O
    GenPeptideConn(PeptideDrug[idx].SeqId[2-S],PeptideDrug[idx].SeqId[3-S],NSB);
    // CA ==> CB and not GLY
    if (curAcidId >= 1 && S <= 1)
        GenPeptideConn(PeptideDrug[idx].SeqId[1-S],PeptideDrug[idx].SeqId[4-S],0);
    switch (curAcidId)
    {
    case 1: // ALA
        break;
    case 2: // VAL (CB==> CG1 and CG2)
    case 7: // THR
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[6-S],NSB);
        break;
    case 3: // LEU ( CB ==> CG) (CG ==> CD1 and CD2)
    case 16: // ASP
    case 18: // ASN
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[6-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[7-S],NSB);
        break;
    case 4: // ILE (CB==> CG1 and CG2) (CG1 ==> CD1)
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[6-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[7-S],NSB);
        break;
    case 5: // SER
    case 6: // CYS (CB ==> CG)
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],0);
        break;
    case 8: // MET (CB=>CG) (CG==>SD) (SD==>CE)
    case 14: // LYS
    case 15: // ARG
    case 17: // GLU
    case 19: // GLN
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[6-S],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[6-S],PeptideDrug[idx].SeqId[7-S],0);
        if (curAcidId==14 || curAcidId==15) // ( CE => NZ)
            GenPeptideConn(PeptideDrug[idx].SeqId[7-S],PeptideDrug[idx].SeqId[8-S],NSB);
        if (curAcidId==15) { // ( CZ => NH1, CZ=>NH2)
            GenPeptideConn(PeptideDrug[idx].SeqId[8-S],PeptideDrug[idx].SeqId[9-S],NSB);
            GenPeptideConn(PeptideDrug[idx].SeqId[8-S],PeptideDrug[idx].SeqId[10-S],NSB);
        }
        if (curAcidId==17 || curAcidId==19) // (CD=> OE2)
            GenPeptideConn(PeptideDrug[idx].SeqId[6-S],PeptideDrug[idx].SeqId[8-S],NSB);
        break;
    case 9: // PRO (CB=>CG) (CG==>CD) (CD==>N)
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[6-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[6-S],PeptideDrug[idx].SeqId[0],NSB);
        break;
    case 10: // PHE (CB==>CG) (CG=>CD1 CG=>CD2) (CD1=>CE1) (CD2=>CE2)
    case 11: // TYR
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[6-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[7-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[6-S],PeptideDrug[idx].SeqId[8-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[7-S],PeptideDrug[idx].SeqId[9-S],NSB);
        // (CE1=>CZ) (CE2=>CZ)
        GenPeptideConn(PeptideDrug[idx].SeqId[8-S],PeptideDrug[idx].SeqId[10-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[9-S],PeptideDrug[idx].SeqId[10-S],NSB);
        if (curAcidId==11)
            GenPeptideConn(PeptideDrug[idx].SeqId[10-S],PeptideDrug[idx].SeqId[11-S],NSB);
        break;
    case 12: // TRP (CB==>CG) (CG=>CD1 CG=>CD2) (CD1=>NE1) (CD2=>CE2) (CD2=>CE3)
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[6-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[7-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[6-S],PeptideDrug[idx].SeqId[8-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[7-S],PeptideDrug[idx].SeqId[9-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[7-S],PeptideDrug[idx].SeqId[10-S],NSB);
        //  (NE1=>CE2) (CE2=>CZ2) (CE3=>CZ3)
        GenPeptideConn(PeptideDrug[idx].SeqId[8-S],PeptideDrug[idx].SeqId[9-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[9-S],PeptideDrug[idx].SeqId[11-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[10-S],PeptideDrug[idx].SeqId[12-S],NSB);
        // CZ2=>CH2 CZ3=>CH2
        GenPeptideConn(PeptideDrug[idx].SeqId[11-S],PeptideDrug[idx].SeqId[13-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[12-S],PeptideDrug[idx].SeqId[13-S],NSB);
        break;
    case 13: // HIS (CB==>CG) (CG=>ND1 CG=>CD2) (ND1=>CE1) (CD2=>NE2) (CE1=>NE2)
        GenPeptideConn(PeptideDrug[idx].SeqId[4-S],PeptideDrug[idx].SeqId[5-S],0);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[6-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[5-S],PeptideDrug[idx].SeqId[7-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[6-S],PeptideDrug[idx].SeqId[8-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[7-S],PeptideDrug[idx].SeqId[9-S],NSB);
        GenPeptideConn(PeptideDrug[idx].SeqId[8-S],PeptideDrug[idx].SeqId[9-S],NSB);
        break;
    }
}


/* ==================================
   delete right spaces in s
 =====================================*/
char* trimSpace(char *s)
{
    int i,L;
    L=strlen(s)-1;
    for (i=L;i>=0;i--)
        if (*(s+i)!=' ')
            break;
    *(s+i+1)='\0';
    return(s);
}

/* ==================================
   delete all spaces in s
 =====================================*/
char* trimAllSpace(char *s)
{
    int i,j=0;
    for (i=0;s[i];i++)
        if (s[i] != ' ') {
            s[j]=s[i];
            j++;
        }
    s[j]='\0';
    return(s);
}

/********************************************************/
/*              Read Amino-acid Model of AMBER                   */
/********************************************************/
void ReadAMBERModel()
{
    int i, j;
    //char buf[90],buf1[30];
    char buf[90];
    FILE * inputFile;

    if((inputFile = fopen("all_atom.txt", "r" ))==NULL) {     /* open config file */
                printf("Amber model file open error !\n") ;
                return;
                //exit(0) ;
    }

    j = 0;

    while (fscanf(inputFile, "%s ", buf) != EOF)
    {
        fscanf(inputFile, "%s \n", buf);
        modelAtomNum[j] = atoi(buf);
        for (i = 1; i <= modelAtomNum[j]; i ++)
        {
        fscanf(inputFile, "%s \n", buf);
        model[j][i].num = atoi(buf);

        fscanf(inputFile, "%s \n", buf);
        sprintf(model[j][i].atomName," %s",buf);

        fscanf(inputFile, "%s \n", buf);
        strcpy(model[j][i].atomType, buf);

        fscanf(inputFile, "%s \n", buf);
        model[j][i].treeType = buf[0];

        fscanf(inputFile, "%s \n", buf);        /* nA,nB,nC: the atoms forming the diAngle */
        model[j][i].nA = atoi(buf);

        fscanf(inputFile, "%s \n", buf);
        model[j][i].nB = atoi(buf);

        fscanf(inputFile, "%s \n", buf);
        model[j][i].nC = atoi(buf);

        fscanf(inputFile, "%s \n", buf);
        model[j][i].length = atof(buf);

        fscanf(inputFile, "%s \n", buf);
        model[j][i].angle = RADIAN(atof(buf));

        fscanf(inputFile, "%s \n", buf);
        model[j][i].diAngle = atof(buf);

        fscanf(inputFile, "%s \n", buf);
        model[j][i].charge = atof(buf);
/*
        printf("\n %d %d %d", model[j][i].nA, model[j][i].nB,
          model[j][i].nC );
        printf(" l=%f a=%f c=%f ang=%f",model[j][i].length,
          model[j][i].angle, model[j][i].diAngle, model[j][i].nC);

        printf(\n %s %f",,model[j][i].,model[j][i].charge);
        printf("  i=%d j=%d",i,j);
        getch();
*/
     }
     j++;
    }
    fclose(inputFile);
}

/*==================================================
    Get atom charge of Receptor Protein from AMBER model
=========================================================*/
double AtomCharge(char * uniqName, int curAcidId)
{
    int i;
    for (i=1; i <= modelAtomNum[curAcidId]; i++)
        if (strcmp(model[curAcidId][i].atomName,uniqName)==0)
            return(model[curAcidId][i].charge);
    return(0);
}


/*==================================================
    Get atom charge of Receptor Protein according to formal charge
=========================================================*/
double AtomFormCharge(char * uniqName, int curAcidId)
{
    char buf[4];
    strncpy( buf, &uniqName[1] , 3);
    buf[3] = '\0';
    if (strcmp(buf,"OXT")==0) return -1;
    //printf("|%s|%s|",uniqName,buf);
    switch (curAcidId)
    {
    case 13: // hisidine
        if (strcmp(buf,"ND1")==0 || strcmp(buf,"NE2")==0 )
            return 0.5;
        else
            return 0.0;
    case 14: // Lysine
        if (strcmp(buf,"NZ")==0)
            return 0.5;
        else
            return 0.0;
    case 15: // Arg
        if (strcmp(buf,"NH1")==0 || strcmp(buf,"NH2")==0)
            return 0.5;
        else
            return 0.0;
    case 16: // Asp
        if (strcmp(buf,"OD1")==0 || strcmp(buf,"OD2")==0)
            return -0.5;
        else
            return 0.0;
    case 17: // Glu
        if (strcmp(buf,"OE1")==0 || strcmp(buf,"OE2")==0)
            return -0.5;
        else
            return 0.0;
     default:
        return 0;
    }
}

/********************************************************/
/*Hydrogen  bond candidate for energy function
1: both
2: acceptor
3: donor
*********************************************************/
int HBond(char *uniqName, int curAcidId)
{
    // if not N, S, O return
    if ( uniqName[1] != 'N' && uniqName[1] != 'O'
      && uniqName[1] != 'S' )
        return 0;

    if (strcmp(uniqName," N")== 0)
            return 3; // donor   (give a Hedrogen atom
    if (strcmp(uniqName," O")== 0)
            return 2; // acceptor (get a Hedrogen atom)
    if (strcmp(uniqName," OXT")== 0)
            return 2; // acceptor (get a Hedrogen atom)

    switch (curAcidId)
    {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 9:
    case 8:
    case 10:
        break;                  // no candidate
    case 5:  // SER
    case 7:  // THR
    case 11: // TYR
    case 16: // ASP
        return 1;             // 1: return -OH  both
    case 6: // Cysteine
        return 5;             // return Cys -S-  disulfide
    case 12:
        return 3;               //   =N- donor
    case 13:  // HIS
        return 1;               // both
     case 14:
        return 3;             // return -N- donor
    case 15:
        if (uniqName[2] == 'E')
            return 3;
        else
            return 3;
//    case 16: // ASP
    case 17: // GLU
        return 2;             // return O=c< acceptor
    case 18:
    case 19:
        if (uniqName[1] == 'O')
            return 2;
        else
            return 3;
    default:
        return 0;
    }

    return 0;
}

/* get amino-acid id from file */
int CatchAcidID(char *string ,int startPos)
{
    char *buf;
    int i;

    buf = substr(string,startPos,3);

    for (i=0;i<26;i++)
        if(strncmp(acidName[i],buf,3)== 0)
            return i;
    if (strcmp("NLE",buf)== 0)
        return 8;  // modified MET
    //if (strcmp("TYR",buf)== 0)
    //    return 
    return -1;
}

#if 1
/********************************************************/
/* ATOM
COLUMNS   DATA TYPE     FIELD       DEFINITION
---------------------------------------------------------------------------------
 1 -  6   Record name   "ATOM  "
 7 - 11   Integer       serial      Atom serial number.
13 - 16   Atom          name        Atom name.
17        Character     altLoc      Alternate location indicator.
18 - 20   Residue name  resName     Residue name.
22        Character     chainID     Chain identifier.
23 - 26   Integer       resSeq      Residue sequence number.
27        AChar         iCode       Code for insertion of residues.
31 - 38   Real(8.3)     x           Orthogonal coordinates for X in Angstroms.
39 - 46   Real(8.3)     y           Orthogonal coordinates for Y in Angstroms.
47 - 54   Real(8.3)     z           Orthogonal coordinates for Z in Angstroms.
55 - 60   Real(6.2)     occupancy   Occupancy.
61 - 66   Real(6.2)     tempFactor  Temperature factor.
73 - 76   LString(4)    segID       Segment identifier, left-justified.
77 - 78   LString(2)    element     Element symbol, right-justified.
79 - 80   LString(2)    charge      Charge on the atom.

Read data into receptor protein (temp).
These data are temprate ane will be stored into
    ReceptorATOMS[MAXATOM] for binding procedure.
*/
/********************************************************/
void RecordAtom(char * line)
{
    char *buf, uniqName[5];
    int curAcidSer;             // (23:26 resSeq)sequential residue number in a protein
    int curAcidId;              // (18:20) current residue ID (0:19, AcidName)
    ProteinAtom  *curAtom;

    buf = substr(line,0,4);
    if (strcmp(buf, "ATOM") !=0) // non "ATOM" label
        return;
    if (line[17-1] == 'B' || line[17 -1] == '2' || line[17-1] == 'U')
        return;

    curAcidId = CatchAcidID(line, 17);
    if (curAcidId == -1) return; // non-standard residues

    // atom Name
    strncpy(uniqName,line+13-1,4);  uniqName[4] = '\0';

    // eliminate H, D atom for now
    if (uniqName[1] == 'H' || uniqName[1] == 'D')
        return;

    // ignore A, T, G, C (i.e. non-standard 20 amino acids
    if (line[18-1] == ' ')
        return;

    trimSpace(uniqName); // delete spaces in uniqname

    // get # of atom and # of residue in protein
    curAcidSer = atoi(substr(line,23-1,4));

    if ( curAcidSer != lastAcidSer) // # residue is changed ?
    { // yes: store data into Acid array (AcidList)

        if (lastAcidSer != -1)      // initial value = -1
            Protein.AcidList[acidCount-1].AtomNum = atomCount;
        acidCount ++;               // total residue count increase 1
        atomCount = 0;              // reset total residue atom to 0
        lastAcidSer = curAcidSer;   // set the current # residue
        Protein.AcidList[acidCount-1].acidID = curAcidId;
        Protein.AcidList[acidCount-1].flag = 0;
    }

    // begin: store the pdb data into Atom array (AtomList)
    atomIndex = (acidCount-1) * MAX_ACID_ATOM + atomCount;
    if(atomIndex >= 900000){ /*printf("\natomIndex:%d\n",atomIndex);*/
        cout<<"Excess max atom";exit(0);
    }
    curAtom = &Protein.AtomList[atomIndex];    // assign address to curAtom
    curAtom->atomName = line[13];
    curAtom->acidID=curAcidId;
    switch(curAtom->atomName)
    {
        case 'O': curAtom->EngType = 0; break;
        case 'N': curAtom->EngType = 1; break;
        case 'C': curAtom->EngType = 2; break;
        case 'S': curAtom->EngType = 3; break;
        case 'P': curAtom->EngType = 4; break;
        default:  curAtom->EngType = 5; // for hydorgen
    }

    curAtom->acidIdx   = acidCount;
    curAtom->AcidSeqId = atoi(substr(line,23-1,4));
    curAtom->AtomSeqId = atoi(substr(line, 7-1,5));

    //strncpy(curAtom->UniN, line+13-1, 4);    curAtom->UniN[4]='\0';
    //trimSpace(curAtom->UniN);


    // store atom position
    curAtom->position.x = atof(substr(line,31-1,8));
    curAtom->position.y = atof(substr(line,39-1,8));
    curAtom->position.z = atof(substr(line,47-1,8));

    curAtom->HbType = HBond(uniqName, curAcidId);

    curAtom->flag=0;
//    strcpy(curAtom->orgStr,line);

    if (curAcidId >= 13 && curAcidId <= 17)  // only for charge residues
        curAtom->charge = AtomFormCharge(uniqName, curAcidId);
    else
        curAtom->charge = 0;

    atomCount++;

#if 0
    printf("\n %4d %s %c %8.3f %3d %8.3f %4s %d",
            acidCount, acidName[curAcidId], curAtom->atomName,
            curAtom->position.x,curAtom->HbType,curAtom->charge,uniqName,curAcidId);
    getchar();
#endif

}


/*******************************************************
Read data into ProteinHET (temp).
These data are temprate ane will be stored into
    ReceptorATOMS[MAXATOM] for binding procedure.
     1: both 2: acceptor 3: donor
*******************************************************/
void RecordProHet(char *line)
{
    Atoms *curAtom;
    char  *buf, *uniqName;

    buf = substr(line, 0, 6);

    if (strcmp(buf,"HETATM") !=0 )
        return;
    if (line[17-1] == 'B' || line[17 -1] == '2' || line[17-1] == 'U')
        return;
    if (line[14-1] == 'H' || line[14-1] == 'D')
        return;

    uniqName = substr(line,17, 3);
    if (!ConsiderWater && strcmp(uniqName,"HOH") == 0)
        return;
        
    uniqName = substr(line,13-1, 4);
    trimAllSpace(uniqName); // delete spaces in uniqname
    if (!ConsiderHetal &&
        (strcmp(uniqName,"MG") == 0 ||
         strcmp(uniqName,"ZN") == 0 ||
         strcmp(uniqName,"MN") == 0 ||
         strcmp(uniqName,"FE") == 0 ||
         strcmp(uniqName,"CA") == 0 ||
         strcmp(uniqName,"CU") == 0)  )
        return;

    curAtom = &ProteinHET[ProHETAtomCount];    // assign address to curAtom
    curAtom->name = line[13];

    switch(curAtom->name)
    {
        case 'O': curAtom->EngType = 0; curAtom->HbType = 2; break;
        case 'N': 
        {        
            if(line[12]=='Z' || line[12]=='M')
            {
                  curAtom->name    = 'M';
                  curAtom->EngType = 6; // for metal atom
            }
            else
            {
                  curAtom->EngType = 1; 
                  curAtom->HbType  = 3;        
            }
            break; 
        }
        case 'C': curAtom->EngType = 2; break;
        case 'S': curAtom->EngType = 3; break;
        default:
                  curAtom->name    = 'A';
                  curAtom->EngType = 6; // for metal atom
    }

    curAtom->AtomSeqId  = atoi(substr(line, 7-1,5));
    curAtom->AcidSeqId  = atoi(substr(line,23-1,4));
    curAtom->AcidType   = -1;
    curAtom->OtherName  = (char *)malloc(5*sizeof(char));
    strcpy(curAtom->OtherName,substr(line,17,3));
    curAtom->pos.x = atof(substr(line,31-1, 8));
    curAtom->pos.y = atof(substr(line,39-1, 8));
    curAtom->pos.z = atof(substr(line,47-1, 8));
    curAtom->charge     = 0;

    uniqName = substr(line,17, 3);
    if (strcmp(uniqName,"HOH") == 0)  // decide the hbond type
    {
        curAtom->HbType = 1; // 1: both 2: acceptor 3: donor
        curAtom->charge = 0;
    }

    uniqName = substr(line,13-1, 4);
    trimAllSpace(uniqName); // delete spaces in uniqname
    if (line[13-1]!=' ' && 
        (strcmp(uniqName,"MG") == 0 ||
         strcmp(uniqName,"ZN") == 0 ||
         strcmp(uniqName,"MN") == 0 ||
         strcmp(uniqName,"FE") == 0 ||
         strcmp(uniqName,"CA") == 0 ||
         strcmp(uniqName,"CU") == 0    ))
    {
        curAtom->HbType     = 3; // 1: both 2: acceptor 3: donor
        curAtom->charge     = 2;
        //printf("\n have metal atom in Protein"); getchar();
        ProMetalAtomCount++;
    }    
#if 0
        printf("\n%c %5s %d %d %d %8.3f %3d %8.3f\n",
               curAtom->name,uniqName,curAtom->AcidType,curAtom->AcidSeqId,curAtom->AtomSeqId,
               curAtom->pos.x,curAtom->HbType,curAtom->charge); getchar();
#endif
 
#ifdef HAVE_PROSTR
    strcpy(OrgProHetSTR[ProHETAtomCount],line);
#endif
    ProHETAtomCount++;
#if 0
    printf(" %s %c %5s %d %8.3f %3d %8.3f\n", line,
        curAtom->name,uniqName,curAtom->AcidType,curAtom->pos.x,curAtom->HbType,curAtom->charge);
#endif
 }
#endif


/********************************************************/
/* 1. read pdb information into  Protein and DRUGATOMS   */
/*    . position, atom, paramerers of energy function   */
/* 2. set the h-bond and disulfide bond                */
/********************************************************/
void ClearPDB();
void ReadPDB(char *fileName)
{
    char line[256],fpath[256], *buf;
    FILE *inFile;
    int len, RRigFlag=0;

    sprintf(fpath,"%s%s",PDBPath,fileName);

    if((inFile=fopen(fpath,"r"))==NULL)    /* open config file */
    {
        printf("Pdb file %s : Open Error !\n", fpath) ;
        exit(1) ;
    }

    //-------- Modification (Optional)--------
    ClearPDB();

    atomIndex = 0;          /* atom index in whole protein */
    acidCount = 0;          /* total residues of receptor protein  */
    atomCount = 0;          /* total atoms of a residue in receptor protein*/
    lastAcidSer       = -1;
    RefLigAtomCount   =  0;
    ProHETAtomCount   =  0; /* total atoms of hetatm (water,etc.) in receptor*/
    ProMetalAtomCount =  0;

    while(fgets(line, 90, inFile)!=NULL)
    {
        buf = substr(line,0,3);

        if (strcmp(buf, "TER")==0)  RRigFlag = 1;

        RecordAtom(line);

        if (RRigFlag == 1)
        {
            len= strlen(line);
            if (ConsiderHetATOM && len >= 82 && line[81] == 'P')
            {
                // read special atom (metal, wator...)
               // printf("\n%s",line); getchar();
                RecordProHet(line);
            }
            if ( len >= 82 && line[81] != 'P' )
            {
                // read reference ligand hetatm for computing RMSD
                RecordRefLig(line);
            }
        }
        line[81]=' ';
    }

    Protein.AcidList[acidCount-1].AtomNum = atomCount;
    Protein.acidNum = acidCount;
    fclose(inFile);
}

/*=============================================
1. read "HETATM" with original ligand  in original
pdb as reference (original X-ray)
=====================================================*/
void RecordRefLig(char *line)
{
    char buf[30];

      strncpy( buf, &line[0] , 6);
      buf[6] = '\0';
      if ( (strcmp(buf, "HETATM") !=0 && strcmp(buf, "ATOM  ") !=0)
         || RefLigAtomCount >= MAXDRUGATOM)
             return; // read "HETATM" in original pdb as reference

    if (line[17-1] == 'B' || line[17 -1] == '2' || line[17-1] == 'U')
        return;

    if (line[14-1] == 'H' || line[14-1] == 'D')
        return;

    strncpy(buf ,&line[31 - 1] , 8);
    buf[8] = '\0';
    RefDrugOrgPosi[RefLigAtomCount].x=atof(buf);

    strncpy(buf ,&line[39 - 1] , 8);
    buf[8] = '\0';
    RefDrugOrgPosi[RefLigAtomCount].y=atof(buf);

    strncpy(buf ,&line[47 - 1] , 8);
    buf[8] = '\0';
    RefDrugOrgPosi[RefLigAtomCount].z=atof(buf);

    strcpy(RefLigAtom[RefLigAtomCount],line);

/*
    printf("\n %d x=%8.3f %8.3f %8.3f ", RefLigAtomCount,
       RefDrugOrgPosi[RefLigAtomCount].x,
       RefDrugOrgPosi[RefLigAtomCount].y,RefDrugOrgPosi[RefLigAtomCount].z);
    getch();
*/
     RefLigAtomCount++;
}

/* =======================================================
 BEGIN: denote the circle atoms on array DrugCircle[]
 type=0: atom is not a circle atom
      1: atom is a a circle atom
=======================================================*/
void DecideDrugCircle()
{
    int i, sbT=0;
    for (i=0;i<DrugAtomCount;i++) {
        DrugCircle[i].type=0; 
        DrugCircle[i].cid=-1;
    }
    for (i=0;i<DrugAtomCount;i++) {
        if (DrugCircle[i].type >= 1 || DrugATOMS[i].numB == 1) continue;
        TestDrugCircle(i);
    }

    if (DrugFormat >= 2) { // update the single bond of mol-format
         for(i=0 ; i<singleB[0].numSB ;i++) {
            if (DrugCircle[singleB[0].atom[i][0]].type==0 ||
                DrugCircle[singleB[0].atom[i][1]].type==0 ||
                DrugCircle[singleB[0].atom[i][0]].cid!=DrugCircle[singleB[0].atom[i][1]].cid ) { // Single bond not on circle
                singleB[0].atom[sbT][0] = singleB[0].atom[i][0];
                singleB[0].atom[sbT][1] = singleB[0].atom[i][1];
                singleB[0].type[sbT]    = singleB[0].type[i];
               // printf("\n %d (%d ==>%d)",sbT,singleB[0].atom[sbT][0],
               //     singleB[0].atom[sbT][1]); getchar();
                sbT++; 
            }
        } // for i
        singleB[0].numSB=sbT;
    } // if
    if (DelTerminalSB) { // delete terminal SB
        sbT=0;
        for(i=0 ; i<singleB[0].numSB ;i++) {
            if (DrugATOMS[singleB[0].atom[i][0]].numB>=2 && 
                        DrugATOMS[singleB[0].atom[i][1]].numB>=2 ) { // delete terminal SB
                singleB[0].atom[sbT][0] = singleB[0].atom[i][0];
                singleB[0].atom[sbT][1] = singleB[0].atom[i][1];
                singleB[0].type[sbT]    = singleB[0].type[i];
               // printf("\n %d (%d ==>%d)",sbT,singleB[0].atom[sbT][0],
               //     singleB[0].atom[sbT][1]); getchar();
                sbT++;
            }
        } // for i
        singleB[0].numSB=sbT;
    }    
}

int TestDrugCircle(int idx)
{
    int c[9],p[9];
    for (c[0]=0;c[0]<DrugATOMS[idx].numB;c[0]++) { // first node
        p[0]=DrugATOMS[idx].adj[c[0]]; // first node id
        if (DrugATOMS[p[0]].numB ==1) continue;
        for (c[1]=0;c[1]<DrugATOMS[p[0]].numB;c[1]++) { // second node
            p[1]=DrugATOMS[p[0]].adj[c[1]]; // second node id
            if (p[1]== idx || DrugATOMS[p[1]].numB ==1) continue;
            for (c[2]=0;c[2]<DrugATOMS[p[1]].numB;c[2]++) { // third node
                p[2]=DrugATOMS[p[1]].adj[c[2]]; // third node id
                if (p[2]== p[0] || DrugATOMS[p[2]].numB ==1) continue;
                for (c[3]=0;c[3]<DrugATOMS[p[2]].numB;c[3]++) { // fourth node
                    p[3]=DrugATOMS[p[2]].adj[c[3]]; // fouth node id
                    if (p[3]== p[1] || DrugATOMS[p[3]].numB ==1) continue;
                    if (TestCircle(idx,4,p)) return(1);
                    for (c[4]=0;c[4]<DrugATOMS[p[3]].numB;c[4]++) { // fifth node
                        p[4]=DrugATOMS[p[3]].adj[c[4]]; // fifth node id
                        if (p[4]== p[2] || DrugATOMS[p[4]].numB ==1) continue;
                        if (TestCircle(idx,5,p)) return(1);
                        for (c[5]=0;c[5]<DrugATOMS[p[4]].numB;c[5]++) { // sixth node
                            p[5]=DrugATOMS[p[4]].adj[c[5]]; // sixth node id
                            if (p[5]== p[3] || DrugATOMS[p[5]].numB ==1) continue;
                            if (TestCircle(idx,6,p)) return(1);
                            for (c[6]=0;c[6]<DrugATOMS[p[5]].numB;c[6]++) { // seventh node
                                p[6]=DrugATOMS[p[5]].adj[c[6]]; // seventh node id
                                if (p[6]== p[4] || DrugATOMS[p[6]].numB ==1) continue;
                                if (TestCircle(idx,7,p)) return(1);
                                for (c[7]=0;c[7]<DrugATOMS[p[6]].numB;c[7]++) { // eightth node
                                    p[7]=DrugATOMS[p[6]].adj[c[7]]; // eight node id
                                    if (p[7]== p[5] || DrugATOMS[p[7]].numB ==1) continue;
                                    if (TestCircle(idx,8,p)) return(1);
                                } // for eightth
                            } // for seventh
                        } // for sixth
                    } // for fifth
                } // for forth
            } // for third
        } // for second
    } // for first
    return(0);
}

int TestCircle(int idx, int count, int aNum[])
{
    int i;
    if (aNum[0] == aNum[count-2]) return(0); // connect to self
    if (idx==aNum[count-1]) {
        //printf("\n BEGIN (%d %d) (%d %d)",DrugATOMS[idx].SeqId,idx,DrugATOMS[aNum[count-1]].SeqId,aNum[count-1]); getch();
        for (i=0;i<count;i++) {
            DrugCircle[aNum[i]].type=count;
            DrugCircle[aNum[i]].cid=idx;
            //printf("\n (%d %d) (%d %d)",DrugATOMS[idx].SeqId,idx,DrugATOMS[aNum[i]].SeqId,aNum[i]); getch();
        }
        //return(1);
    }
    return(0);
}

/* =======================================================
 END: denote the circle atoms on array DrugCircle[]
=======================================================*/


/*===============================================
BEGIN: decide atom charge and H-BOND of Ligand
================================================*/
int GetEndAtomNum(int idx,char AType)
{
    int j,k, Count=0;
    if (DrugATOMS[idx].numB >=2)
        for (j=0;j<DrugATOMS[idx].numB;j++) {
            k=DrugATOMS[idx].adj[j];
            if (DrugATOMS[k].name==AType && DrugATOMS[k].numB==1)
                Count++;
        }
   return(Count);
}


int GetCircleAtomSer(int idx)
{
    char buf[10];
    int j,k=0;
    strncpy( buf , &OrgDrugSTR[idx][13 - 1] , 4);
    buf[4] = '\0';
    for (j=0;j<4;j++)
        if ( buf[j] >= '0' && buf[j] <= '9') {
            k= buf[j]-'0';
            break;
        }
    return(k);
}

int GetCircleAtomType(int idx, char AType)
{
    int j,k, p=0;
    if (DrugATOMS[idx].numB >=2)
        for (j=0;j<DrugATOMS[idx].numB;j++) {
            k=DrugATOMS[idx].adj[j];
            if (DrugATOMS[k].name==AType) {
                p=GetCircleAtomSer(idx);
                break;
            }
        }
   return(p);
}
/*==================================================
   decide the h-bond type and formal charge of an atom of Ligand
1: both 2: acceptor 3: donor
======================================================*/
void DecideDrugAtomHBondType()
{
    int i,idx, p3,PO1,PO2;
    //double len,ang;
    double len;
    char buf[10];

    LigHbondCount=0;
    LigChargeAtomCount=0;
    for (i=0;i<DrugAtomCount;i++) {
        DrugATOMS[i].charge=0;
        DrugATOMS[i].HbType=0;
    }
    for (i=0;i<DrugAtomCount;i++) {
        strncpy(buf, &OrgDrugSTR[i][0] , 4);
        buf[4] = '\0';
        if (strcmp(buf,"ATOM")==0) {
            idx = CatchAcidID( OrgDrugSTR[i], 17);
            strncpy( buf , &OrgDrugSTR[i][13 - 1] , 4);
            buf[4] = '\0';
            trimSpace(buf);
            DrugATOMS[i].HbType= HBond(buf, idx);
            DrugATOMS[i].charge=AtomFormCharge(buf, idx);
        } else {
        /*
        if (DrugATOMS[i].name == 'I') {
            DrugATOMS[i].charge=0.2;
            DrugATOMS[i].HbType=1;
        } */
        if (DrugATOMS[i].name == 'N') {
            switch(DrugATOMS[i].numB) {
                case 4:
                    DrugATOMS[i].charge=1;
                    break;
                case 3:
                    DrugATOMS[i].HbType=2; // acceptor
                  //  printf("\n hit"); getchar();
                    break;
                case 1:
                    DrugATOMS[i].HbType=3; // donor
                    idx=DrugATOMS[i].adj[0];
                    p3=GetEndAtomNum(idx,'N');
                    //printf("|%d|%d|%d|",i,idx,p3);
                    if (p3 >= 2)
                        DrugATOMS[i].charge=1;
                    break;
                case 2: // on TextBook(Mathews): BioChemistry p. 87
                    //printf("|%d|",DrugCircle[i].type);
                    DrugATOMS[i].HbType=3;
                    if (DrugCircle[i].type ==0) { // not on circle
                        PO1=DrugATOMS[i].adj[0];
                        PO2=DrugATOMS[i].adj[1];
                        if (DrugATOMS[PO1].name=='N' && DrugATOMS[PO2].name == 'N'
                            && (DrugATOMS[PO1].numB==1 || DrugATOMS[PO2].numB==1) ) { // N=N=N
                            if (DrugATOMS[PO1].numB==1) DrugATOMS[PO1].charge=-1;
                            else  DrugATOMS[PO2].charge=-1;
                            DrugATOMS[i].charge=1.0;
                        }
                        else
                            DrugATOMS[i].HbType=3; // Donor
                    }
                    else if (DrugCircle[i].type ==6) {// 6-atom circle
                        p3=GetCircleAtomSer(i);
                        PO1=GetCircleAtomType(DrugATOMS[i].adj[0], 'O');
                        PO2=GetCircleAtomType(DrugATOMS[i].adj[1], 'O');
                       // printf("|%d|%d|%d|",p3,PO1,PO2);
                        if (p3==1) {
                            if (PO1==0 && PO2==0) //neighbor C not connecting to "O"
                                DrugATOMS[i].HbType=2; // acceptor
                            else
                                DrugATOMS[i].HbType=3; // Donor
                        }
                        if (p3==3) {
                            if (PO1==0 || PO2==0) //neighbor C connecting to "O"
                                DrugATOMS[i].HbType=2; // acceptor
                            else
                                DrugATOMS[i].HbType=3; // Donor
                        }
                        if (p3==5 || p3 == 8) DrugATOMS[i].HbType=2; // Acceptor
                    }
                    else if (DrugCircle[i].type ==5) { // 5-atom circle
                        p3=GetCircleAtomSer(i);
                        if (p3==7 || p3 == 3) DrugATOMS[i].HbType=2; // Acceptor
                        else DrugATOMS[i].HbType=3; // Donor
                        //if (p3==9 || p3 == 1) DrugATOMS[i].HbType=3; // Donor
                    }
                    break;
            }
        } // if "N"
        if (DrugATOMS[i].name == 'O') {
            if (DrugATOMS[i].numB == 2) DrugATOMS[i].HbType=2; // acceptor
            if (DrugATOMS[i].numB == 1) {
                DrugATOMS[i].HbType=2; // acceptor
                idx=DrugATOMS[i].adj[0];
                len= sqrt(VectorLength(DrugATOMS[i].pos,DrugATOMS[idx].pos));
                p3=GetEndAtomNum(idx,'O');
                switch(DrugATOMS[idx].name) {
                    case 'C':
                        if (p3==2) {
                            DrugATOMS[i].HbType=2;
                            DrugATOMS[i].charge=-1;
                        }
                        else if (len > 1.35) // SB
                            DrugATOMS[i].HbType=1; // both
                        else DrugATOMS[i].HbType=2; // acceptor
                        break;
                    case 'P':
                        if (p3>=2) { // PO2, PO3, PO4
                            DrugATOMS[i].HbType=2;
                            DrugATOMS[i].charge=-1;
                        }
                        else
                            DrugATOMS[i].HbType=2; // acceptor
                        break;
                    case 'N':
                        if (p3>=2) {
                            DrugATOMS[i].HbType=2;
                            DrugATOMS[idx].charge=1.0; // N02
                            //printf("\n \n UUUU"); getch();
                        }
                        else
                            DrugATOMS[i].HbType=2; // acceptor
                        break;
                    case 'S':
                        if (p3>=3) { // SO3, SO4
                            DrugATOMS[i].HbType=2;
                            DrugATOMS[i].charge=-1;
                        }
                        else
                            DrugATOMS[i].HbType=2; // acceptor
                        break;
                } // switch
            } // if Num==1
         } // if "O"
       } // else
        if (DrugATOMS[i].HbType>0) LigHbondCount++;
        if (DrugATOMS[i].charge !=0) LigChargeAtomCount++;
    } // for
/*
    for (i=0;i<DrugAtomCount;i++) {
        printf("\n %d %c %d Adj=%d %d",i, DrugATOMS[i].name,
            DrugATOMS[i].numB, DrugATOMS[i].adj[0], DrugATOMS[i].HbType);
        getch();
    }
*/
}

/*
======================================================================
store data into Drug structure (DrugATOMS[]) directly
the format is the same as the "ATOM" Label
COLUMNS     DATA TYPE       FIELD          DEFINITION
----------------------------------------------------------------------
 1 -  6     Record name     "HETATM"
 7 - 11     Integer         serial         Atom serial number.
13 - 16     Atom            name           Atom name.
17          Character       altLoc         Alternate location indicator.
18 - 20     Residue name    resName        Residue name.
22          Character       chainID        Chain identifier.
23 - 26     Integer         resSeq         Residue sequence number.
27          AChar           iCode          Code for insertion of residues.
31 - 38     Real(8.3)       x              Orthogonal coordinates for X.
39 - 46     Real(8.3)       y              Orthogonal coordinates for Y.
47 - 54     Real(8.3)       z              Orthogonal coordinates for Z.
55 - 60     Real(6.2)       occupancy      Occupancy.
61 - 66     Real(6.2)       tempFactor     Temperature factor.
73 - 76     LString(4)      segID          Segment identifier;
77 - 78     LString(2)      element        Element symbol; right-justified.
79 - 80     LString(2)      charge         Charge on the atom.
======================================================================
*/
void RecordHetAtm(char *line)
{
    char  *buf;
    Atoms *curAtom;
    int CAcidId;

    buf = substr(line, 0, 6);

    if ( strcmp(buf, "HETATM") !=0 &&
         strcmp(buf, "ATOM  ") !=0)
        return; // non "HETATM" and "ATOM  " label

    if (line[17-1] == 'B' || line[17-1] == '2' || line[17-1] == 'U')
        return;
    if (line[14-1] == 'H' || line[14-1] == 'D')
        return;

    curAtom = &DrugATOMS[DrugAtomCount];  // assign address to curAtom

    curAtom->name = line[13];

    switch(curAtom->name)
    {
        case 'O':
        case 'I': curAtom->EngType = 0; break;

        case 'N': curAtom->EngType = 1; break;
        case 'C': curAtom->EngType = 2; break;
        case 'S': curAtom->EngType = 3; break;
        case 'P': curAtom->EngType = 4; break;
        default:  curAtom->EngType = 5;
    }

    //////////////////////////////
    //
    //  SeqId;          /* Sequence atom number in whole molecule */
    //  AtomSeqId;      /* Sequence atom number */
    //  AcidSeqId;      /* Sequence acid number in a protein */
    //

    //curAtom->SeqId     = DrugAtomCount+1;
    //buf = substr(line+ 6, 5);   curAtom->AtomSeqId = atoi(buf);
    //buf = substr(line+22, 4);   curAtom->AcidSeqId = atoi(buf);

    curAtom->SeqId     = DrugAtomCount+1;
    curAtom->AtomSeqId = atoi(substr(line, 7-1, 5));
    curAtom->AcidSeqId = atoi(substr(line,23-1, 4));
    curAtom->AcidType  = (CAcidId = CatchAcidID( line, 17));

    // store atom position
    curAtom->pos.x  = atof(substr(line,31-1, 8));
    curAtom->pos.y  = atof(substr(line,39-1, 8));
    curAtom->pos.z  = atof(substr(line,47-1, 8));

    curAtom->charge = 0.0;
    curAtom->numB   = 0;

    DrugOrgPosi[DrugAtomCount].x=curAtom->pos.x;
    DrugOrgPosi[DrugAtomCount].y=curAtom->pos.y;
    DrugOrgPosi[DrugAtomCount].z=curAtom->pos.z;

    strncpy(OrgDrugSTR[DrugAtomCount]  ,line,80);
    strcpy(OrgDrugSTR[DrugAtomCount]+80,"\n");

    buf = substr(line, 7-1, 5);
    if( strcmp(buf, "ATOM  ") ==0 && CAcidId !=-1 )
    {
        char *uniqName;
        int  curAcidSer,k;

        //////////////////////////////
        //
        // peptide drug
        //
        uniqName = substr(line,13-1, 5);
        trimSpace(uniqName);    // delete spaces in uniqname

        curAtom->charge = AtomFormCharge(uniqName, CAcidId);
        curAtom->HbType = HBond(uniqName, CAcidId);

        curAcidSer = atoi(substr(line,23-1, 4));

        if( curAcidSer != lastAcidSer)  // # residue is changed ?
        {
            lastAcidSer = curAcidSer;   // set the current # residue

            PeptideDrug[PepResCount].AcidType = CAcidId;
            PeptideDrug[PepResCount].AcidSer  = curAcidSer;
            PeptideDrug[PepResCount].AtomNum  = 0;

            trimAllSpace(uniqName);

            switch(uniqName[0])
            {
                case 'N':
                    PeptideDrug[PepResCount].Start = 0;
                    break;
                case 'C':
                    PeptideDrug[PepResCount].Start = (uniqName[1]) ? 1 : 2;
                    // 1 : (such 1hef) lacking N
                    // 2 : lacking N,CA
                    break;
                default:
                    //-------- Modification --------
                    printf("\n\n *** peptide drug error ***");
                    exit(1);
            }

            // total peptide residue count increase 1
            PepResCount++;
        }

        k = PeptideDrug[PepResCount-1].AtomNum;
        PeptideDrug[PepResCount-1].Atomname[k] = curAtom->name;
        PeptideDrug[PepResCount-1].SeqId[k]    = curAtom->AtomSeqId;
        PeptideDrug[PepResCount-1].AtomNum++;
    }

    DrugAtomCount++;
    OrgDrugSTRIndex = DrugAtomCount;

}

void RecordHetAtm2(char *line)
{
    char  *buf;
    Atoms *curAtom;
    int CAcidId;

    buf = substr(line, 0, 6);

    if ( strcmp(buf, "HETATM") !=0 &&
         strcmp(buf, "ATOM  ") !=0)
        return; // non "HETATM" and "ATOM  " label

    if (line[17-1] == 'B' || line[17-1] == '2' || line[17-1] == 'U')
        return;
    if (line[14-1] == 'H' || line[14-1] == 'D')
        return;

//-------- Modification --------
    if( DrugAtomCount > 100) 
    {
        isExceedMAXDRUGATOM = true;
        ostringstream oss;
        oss << "Error: Atom Number of Drug Exceeding Limit (" << DrugAtomCount << ")" << endl
            << " Skip docking drug - " << currentDrug << endl;
        UniversalErrorMessageLogger(oss.str());
        return;
    }

    curAtom = &DrugATOMS[DrugAtomCount];  // assign address to curAtom

    curAtom->name = line[13];

    switch(curAtom->name)
    {
        case 'O':
        case 'I': curAtom->EngType = 0; break;

        case 'N': curAtom->EngType = 1; break;
        case 'C': curAtom->EngType = 2; break;
        case 'S': curAtom->EngType = 3; break;
        case 'P': curAtom->EngType = 4; break;
        default:  curAtom->EngType = 5;
    }

    //////////////////////////////
    //
    //  SeqId;          /* Sequence atom number in whole molecule */
    //  AtomSeqId;      /* Sequence atom number */
    //  AcidSeqId;      /* Sequence acid number in a protein */
    //

    //curAtom->SeqId     = DrugAtomCount+1;
    //buf = substr(line+ 6, 5);   curAtom->AtomSeqId = atoi(buf);
    //buf = substr(line+22, 4);   curAtom->AcidSeqId = atoi(buf);

    curAtom->SeqId     = DrugAtomCount+1;
    curAtom->AtomSeqId = atoi(substr(line, 7-1, 5));
    curAtom->AcidSeqId = atoi(substr(line,23-1, 4));
    curAtom->AcidType  = (CAcidId = CatchAcidID( line, 17));

    // store atom position
    curAtom->pos.x  = atof(substr(line,31-1, 8));
    curAtom->pos.y  = atof(substr(line,39-1, 8));
    curAtom->pos.z  = atof(substr(line,47-1, 8));

    curAtom->charge = 0.0;
    curAtom->numB   = 0;

    DrugOrgPosi[DrugAtomCount].x=curAtom->pos.x;
    DrugOrgPosi[DrugAtomCount].y=curAtom->pos.y;
    DrugOrgPosi[DrugAtomCount].z=curAtom->pos.z;

    strncpy(OrgDrugSTR[DrugAtomCount]  ,line,80);
    strcpy(OrgDrugSTR[DrugAtomCount]+80,"\n");

    buf = substr(line, 7-1, 5);
    if( strcmp(buf, "ATOM  ") ==0 && CAcidId !=-1 )
    {
        char *uniqName;
        int  curAcidSer,k;

        //////////////////////////////
        //
        // peptide drug
        //
        uniqName = substr(line,13-1, 5);
        trimSpace(uniqName);    // delete spaces in uniqname

        curAtom->charge = AtomFormCharge(uniqName, CAcidId);
        curAtom->HbType = HBond(uniqName, CAcidId);

        curAcidSer = atoi(substr(line,23-1, 4));

        if( curAcidSer != lastAcidSer)  // # residue is changed ?
        {
            lastAcidSer = curAcidSer;   // set the current # residue

            PeptideDrug[PepResCount].AcidType = CAcidId;
            PeptideDrug[PepResCount].AcidSer  = curAcidSer;
            PeptideDrug[PepResCount].AtomNum  = 0;

            trimAllSpace(uniqName);

            switch(uniqName[0])
            {
                case 'N':
                    PeptideDrug[PepResCount].Start = 0;
                    break;
                case 'C':
                    PeptideDrug[PepResCount].Start = (uniqName[1]) ? 1 : 2;
                    // 1 : (such 1hef) lacking N
                    // 2 : lacking N,CA
                    break;
                default:
                    //-------- Modification --------
                    ostringstream oss;
                    oss << "Error: Peptide drug " << currentDrug << " Reading Error" << endl;
                    UniversalErrorMessageLogger(oss.str());
                    //printf("\n\n *** peptide drug error ***");
                    exit(1);
            }

            // total peptide residue count increase 1
            PepResCount++;
        }

        k = PeptideDrug[PepResCount-1].AtomNum;
        PeptideDrug[PepResCount-1].Atomname[k] = curAtom->name;
        PeptideDrug[PepResCount-1].SeqId[k]    = curAtom->AtomSeqId;
        PeptideDrug[PepResCount-1].AtomNum++;
    }

    DrugAtomCount++;
    OrgDrugSTRIndex = DrugAtomCount;

}

/*
===================== CONN ==========================================
COLUMNS     DATA TYPE   FIELD       DEFINITION
---------------------------------------------------------------------
 1 -  6     Record name "CONECT"
 7 - 11     Integer     serial      Atom serial number
12 - 16     Integer     serial      Serial number of bonded atom
17 - 21     Integer     serial      Serial number of bonded atom
22 - 26     Integer     serial      Serial number of bonded atom
27 - 31     Integer     serial      Serial number of bonded atom
32 - 36     Integer     serial      Serial number of hydrogen bonded
37 - 41     Integer     serial      Serial number of hydrogen bonded
42 - 46     Integer     serial      Serial number of salt bridged
47 - 51     Integer     serial      Serial number of hydrogen bonded
52 - 56     Integer     serial      Serial number of hydrogen bonded
57 - 61     Integer     serial      Serial number of salt bridged
* construct the connections of drug molecular
======================================================================
*/
void RecordConn(char * line)
{
    char *buf;
    //int  Tid, baseId, Id, count, s;
    int  Tid, Id, count, s;

    buf = substr(line, 0, 6);
    if (strcmp(buf, "CONECT") !=0) return; // non "CONECT" label

    Tid = atoi(substr(line, 7-1, 5));
    Id  = GetDrugIdx(Tid);
    if (Id < 0)
        return;     // connections not in the Drug Molecular

    s     = 11;         // the start position of a connection
    count =  0;
    while (1)
    {
        Tid = atoi(substr(line, s, 5));
        if (Tid <=0 || s > 80) break;

        Tid = GetDrugIdx(Tid);
        if (Tid < 0)
        {
            s=s+5;
            continue;
        }

        if (count > 5) printf("connection error");

        DrugATOMS[Id].adj[count++] = Tid;
        s=s+5;
    }

    DrugATOMS[Id].numB = count;

    strncpy(OrgDrugSTR[OrgDrugSTRIndex],&line[0],80);
    strcpy(OrgDrugSTR[OrgDrugSTRIndex]+80,"\n");
    OrgDrugSTRIndex++;
}


bool ReadSingleBond(char *line)
{
    char *buf;
    //int baseId,Tid,Tid2, count;
    int Tid,Tid2, count;

    buf = substr(line, 0, 6);
    if (strcmp(buf, "SINGLE") !=0) return false; // non "SINGLE" label

    Tid  = atoi(substr(line, 7-1, 6));
    Tid2 = atoi(substr(line,12-1, 6));

    if (Tid  < DrugATOMS[0].AtomSeqId ||
        Tid  > DrugATOMS[DrugAtomCount-1].AtomSeqId ||
        Tid2 < DrugATOMS[0].AtomSeqId ||
        Tid2 > DrugATOMS[DrugAtomCount-1].AtomSeqId)
    {
        printf( "\n\n\n**** Connction assigned :"
                "Ligand single bond error %d ==> %d****", Tid, Tid2);
        //exit(1);     // connections not in the Drug Molecular
        return false;     // connections not in the Drug Molecular
    }

    count = singleB[0].numSB;

    Tid  = GetDrugIdx(Tid);
    Tid2 = GetDrugIdx(Tid2);
    singleB[0].atom[count][0] = Tid;
    singleB[0].atom[count][1] = Tid2;

    buf = substr(line,18-1, 3);

    if      (strcmp(buf,"3 3") == 0)    singleB[0].type[count]=1;
    else if (strcmp(buf,"2 3") == 0 ||
             strcmp(buf,"3 2") == 0)    singleB[0].type[count]=2;
    else if (strcmp(buf,"2 2") == 0)    singleB[0].type[count]=3;
    else                                singleB[0].type[count]=0;
    /*
       printf("\n (SB) C=%3d A1=%3d A2=%3d Ty=%3d",
            count,singleB[0].atom[count][0],
       singleB[0].atom[count][1], singleB[0].type[count]);
       getch();
    */

    singleB[0].numSB++;

    strncpy(OrgDrugSTR[OrgDrugSTRIndex],&line[0],80);
    strcpy(OrgDrugSTR[OrgDrugSTRIndex]+80,"\n");
    OrgDrugSTRIndex++;
    return true;
}


/*
======================================================================
read ligand informtion: Atom, CONNCT, Single bond
Atom:        ATOM or HETATOM
Connect:     CONECT 3040 3039
single Bond: SINGLE 3033 3031 3 2
======================================================================
*/

bool ReadLigand(FILE *fp)
{
    char line[100];
    int i;

    DrugAtomCount    = 0;   /* total atoms of ligand molecular */
    singleB[0].numSB = 0;   /* total number of single bond */

    lastAcidSer = -1;
    PepResCount =  0;
    while(fgets(line, 95, fp)!=NULL)
    {
        RecordHetAtm(line);
        RecordConn(line);
        ReadSingleBond(line);

    }
    mole[0].count = DrugAtomCount;      // total atoms of drug
    mole[0].atom  = DrugATOMS;
    dummy[0].atom = DummyDrugATOMS;


    if (PepResCount > 0)
    { /* add the connections and SB of peptide residue */
      //  printf("kakak"); getchar();
        for (i=0;i<PepResCount;i++)
        {
            GenPeptideResidueConn(i,PeptideDrug[i].AcidType);

            if( i>0 &&
                PeptideDrug[i].AcidSer-PeptideDrug[i-1].AcidSer ==1 &&
                PeptideDrug[i].SeqId[0]-
                PeptideDrug[i-1].SeqId[PeptideDrug[i-1].AtomNum-1]==1)
            {
                // (N=>C)
                GenPeptideConn(
                        PeptideDrug[i].SeqId[0],
                        PeptideDrug[i-1].SeqId[2-PeptideDrug[i-1].Start],0
                        );
            }
        }
    }

#if 0
    for(int l=0;i<OrgDrugSTRIndex;l++)
        printf("%s",OrgDrugSTR[l]);
#endif
    return true;
}

bool ReadLigand2(FILE *fp)
{
    char line[100];
    int i;

    DrugAtomCount    = 0;   /* total atoms of ligand molecular */
    singleB[0].numSB = 0;   /* total number of single bond */

    lastAcidSer = -1;
    PepResCount =  0;
    isExceedMAXDRUGATOM = false;
    while(fgets(line, 95, fp)!=NULL)
    {
        RecordHetAtm2(line);
        RecordConn(line);
        ReadSingleBond(line);

        //-------- Modification --------
        if( isExceedMAXDRUGATOM ) return false;
    }
    mole[0].count = DrugAtomCount;      // total atoms of drug
    mole[0].atom  = DrugATOMS;
    dummy[0].atom = DummyDrugATOMS;


    if (PepResCount > 0)
    { /* add the connections and SB of peptide residue */
      //  printf("kakak"); getchar();
        for (i=0;i<PepResCount;i++)
        {
            GenPeptideResidueConn(i,PeptideDrug[i].AcidType);

            if( i>0 &&
                PeptideDrug[i].AcidSer-PeptideDrug[i-1].AcidSer ==1 &&
                PeptideDrug[i].SeqId[0]-
                PeptideDrug[i-1].SeqId[PeptideDrug[i-1].AtomNum-1]==1)
            {
                // (N=>C)
                GenPeptideConn(
                        PeptideDrug[i].SeqId[0],
                        PeptideDrug[i-1].SeqId[2-PeptideDrug[i-1].Start],0
                        );
            }
        }
    }

#if 0
    for(int l=0;i<OrgDrugSTRIndex;l++)
        printf("%s",OrgDrugSTR[l]);
#endif
    return true;
}




/*
======================================================================
    read Protein[] into ReceptorATOMS[MAXATOM] based on BIND_DIST
======================================================================
*/
void DecideBindingAtoms()
{
    int j,k;
    //double d,minD;
    Atoms *recAtom;

    ReceptorAtomCount = 0;
    for (j=0;j<Protein.acidNum;j++)
    {
        for (k=0;k<Protein.AcidList[j].AtomNum;k++)
        {
            int t = j*MAX_ACID_ATOM+k;
            recAtom = &ReceptorATOMS[ReceptorAtomCount];

            recAtom->name       = Protein.AtomList[t].atomName;
            recAtom->pos        = Protein.AtomList[t].position;
            recAtom->AtomSeqId  = Protein.AtomList[t].AtomSeqId;
            recAtom->AcidSeqId  = Protein.AtomList[t].AcidSeqId;
//          recAtom->SeqId      = Protein.AtomList[t].acidIdx;
            recAtom->AcidType   = Protein.AtomList[t].acidID;
            recAtom->HbType     = Protein.AtomList[t].HbType;
            recAtom->EngType    = Protein.AtomList[t].EngType;
            recAtom->charge     = Protein.AtomList[t].charge;
            recAtom->Weight     = 1;
            recAtom->WeiType    = ' ';

            ReceptorAtomCount++;
            Protein.AtomList[t].flag = 1;
            Protein.AcidList[j].flag = 1;
        }
    }

    // BEGIN: add hetatm atoms (metal, wator,...) of protein

    for (j=0;ConsiderHetATOM && j<ProHETAtomCount;j++)
    {
        recAtom = &ReceptorATOMS[ReceptorAtomCount];

        recAtom->name       = ProteinHET[j].name;
        recAtom->pos        = ProteinHET[j].pos;
        recAtom->AtomSeqId  = ProteinHET[j].AtomSeqId;
        recAtom->AcidSeqId  = ProteinHET[j].AcidSeqId;
        recAtom->AcidType   = ProteinHET[j].AcidType;
        recAtom->OtherName  = ProteinHET[j].OtherName;
        recAtom->HbType     = ProteinHET[j].HbType;
        recAtom->EngType    = ProteinHET[j].EngType;
        recAtom->charge     = ProteinHET[j].charge;
        recAtom->Weight     = 1;
        recAtom->WeiType    = ' ';

        ReceptorAtomCount++;
    }


    // END: add hetatm atoms (metal, wator,...) of protein

#if 0
    for (i=0;i<ReceptorAtomCount;i++) {
            printf("\n%4d %c %s %8.3f %8.3f %8.3f %8.3f H=%d V=%d",
            ReceptorATOMS[i].AcidSeqId, ReceptorATOMS[i].name,
            acidName[ReceptorATOMS[i].AcidType],ReceptorATOMS[i].pos.x,
            ReceptorATOMS[i].pos.y,ReceptorATOMS[i].pos.z,
            ReceptorATOMS[i].charge,ReceptorATOMS[i].HbType,
            ReceptorATOMS[i].AcidType);
    }
    getchar();
    printf("P%d PH%d D%d RA%d\n",
            Protein.acidNum,ProHETAtomCount,DrugAtomCount,ReceptorAtomCount);
#endif

    mole[1].count=ReceptorAtomCount;
    mole[1].atom=ReceptorATOMS;
    dummy[1].atom=DummyRecATOMS;

    //printf("\n total atoms of recptor=%d Drug=%d",ReceptorAtomCount,DrugAtomCount);
}

/*==========================================================
    assigning search bounded area by set range[][3] and BOX[][3]
    range[0~2][]: the x,y,z of search box
    [][0]: lower bounded  [][1]: upper bouned
    range[3~max][]: other ([][0]:0 and [][1]: 2PI)

=================================================== */
void DecideBoundArea()
{
    int i;
    double t;

    // =====BEGIN=====Set the search area ==============
    BOX[0][0]=ReceptorATOMS[0].pos.x;     // lower bound
    BOX[0][1]=ReceptorATOMS[0].pos.x;     // upper bound
    BOX[1][0]=ReceptorATOMS[0].pos.y;
    BOX[1][1]=ReceptorATOMS[0].pos.y;
    BOX[2][0]=ReceptorATOMS[0].pos.z;
    BOX[2][1]=ReceptorATOMS[0].pos.z;

  //  NO_VARS=6+singleB[0].numSB;  /* total NO_VARS */
    for (i=0;i<ReceptorAtomCount;i++) {
        if (ReceptorATOMS[i].pos.x < BOX[0][0])
            BOX[0][0]=ReceptorATOMS[i].pos.x;
        else if (ReceptorATOMS[i].pos.x > BOX[0][1])
            BOX[0][1]=ReceptorATOMS[i].pos.x;
        if (ReceptorATOMS[i].pos.y < BOX[1][0])
            BOX[1][0]=ReceptorATOMS[i].pos.y;
        else if (ReceptorATOMS[i].pos.y > BOX[1][1])
            BOX[1][1]=ReceptorATOMS[i].pos.y;
        if (ReceptorATOMS[i].pos.z < BOX[2][0])
            BOX[2][0]=ReceptorATOMS[i].pos.z;
        else if (ReceptorATOMS[i].pos.z > BOX[2][1])
            BOX[2][1]=ReceptorATOMS[i].pos.z;
    }
    //=======END===== Set the search area ==========


    for(i=0 ; i<3 ; i++) {  //
        BOX[i][1]=BOX[i][1]+BoundedLen;
        BOX[i][0]=BOX[i][0]-BoundedLen;
        ProteinCenter[i] =(BOX[i][1]+BOX[i][0])/2.0;
        t= fabs(BOX[i][1]-BOX[i][0])/2.0;
       // printf("\n%2d %8.3f %8.3f (%8.3f)",i, BOX[i][0],BOX[i][1],C[i]);
        /*== adjust range[] to boundary of original point */
        range[i][0]=-t;
        range[i][1]=t;
       // printf("(%8.3f %8.3f)",range[i][0],range[i][1]); getch();
    }

    for(i=3 ; i<6 ; i++) {
        range[i][0]=-M_PI ;
        range[i][1]=M_PI;
    }

}

void AdjustReceptorPosi() {
    int i;
    // adjust the center of Protein to origin point
    for (i=0;i<ReceptorAtomCount;i++) {
        ReceptorATOMS[i].pos.x=ReceptorATOMS[i].pos.x-ProteinCenter[0];
        ReceptorATOMS[i].pos.y=ReceptorATOMS[i].pos.y-ProteinCenter[1];
        ReceptorATOMS[i].pos.z=ReceptorATOMS[i].pos.z-ProteinCenter[2];
       // printf("\n After (Portein) ==> (%8.3f %8.3f %8.3f)",ProteinCenter[0],ProteinCenter[1], ProteinCenter[2]);
       // printf(" %8.3f %8.3f %8.3f",ReceptorATOMS[i].pos.x,ReceptorATOMS[i].pos.y,ReceptorATOMS[i].pos.z); getch();
    }
    /*== adjust the BOX[] to original point */
    for (i=0;i<3;i++) {
        BOX[i][0]=range[i][0];
        BOX[i][1]=range[i][1];
        range[i][0]=range[i][0]+3;
        range[i][1]=range[i][1]-3;
    }
}


void SetSBRange(){
    int i;
    NO_VARS=6+singleB[0].numSB;  /* total NO_VARS */
    for(i=6 ; i<NO_VARS ; i++) {
        range[i][0]=-M_PI/2.0;
        range[i][1]=M_PI/2.0;
    }
}


/*==================================================
    Denote the atoms connecting to a atom of durg molecular within 3 links
    DrugAtomAdj[][]: denote the atom list
    DrugAtomConnNum[]: denote the number of atoms
    for computing ligand energy
======================================================*/
void DecideDrugConn()
{
    int i,j,k,m,tot,id,id1;
    for (i=0;i<DrugAtomCount;i++) {
        tot=0;
        for (j=0;j<DrugATOMS[i].numB;j++) {
            id=DrugATOMS[i].adj[j];
            if (! CheckConnRep(&DrugAtomAdj[i][0],tot,i,id)) {
                DrugAtomAdj[i][tot]=id;
                tot++;
            }
            for (k=0;k<DrugATOMS[id].numB;k++) {
                id1=DrugATOMS[id].adj[k];
                if (! CheckConnRep(&DrugAtomAdj[i][0],tot,i,id1)) {
                    DrugAtomAdj[i][tot]=id1;
                    tot++;
                }

                for (m=0;m<DrugATOMS[id1].numB;m++) {
                    if (! CheckConnRep(&DrugAtomAdj[i][0],tot,i,DrugATOMS[id1].adj[m])) {
                        DrugAtomAdj[i][tot]=DrugATOMS[id1].adj[m];
                        tot++;
                    }
                }

            }
        }
        DrugAtomConnNum[i]=tot;
       // printf("\n tot=%d",tot ); getch();
    }
/*
    for (i=0;i<DrugAtomCount;i++) {
        printf("\n %d %3d",i, DrugAtomConnNum[i]);
        for (j=0;j<DrugAtomConnNum[i];j++)
            printf("%3d",DrugAtomAdj[i][j]);
        getch();
    }
*/
}

/*==================================================
Change centers of drug and protein to the original point
======================================================*/
void AdjustDrugPosi()
{
    int i;
    double rect[3][2];

    rect[0][0]=DrugATOMS[0].pos.x;     // lower bound
    rect[0][1]=DrugATOMS[0].pos.x;     // upper bound
    rect[1][0]=DrugATOMS[0].pos.y;
    rect[1][1]=DrugATOMS[0].pos.y;
    rect[2][0]=DrugATOMS[0].pos.z;
    rect[2][1]=DrugATOMS[0].pos.z;

    for (i=0;i<DrugAtomCount;i++) {
        if (DrugATOMS[i].pos.x < rect[0][0])
            rect[0][0]=DrugATOMS[i].pos.x;
        else if (DrugATOMS[i].pos.x > rect[0][1])
            rect[0][1]=DrugATOMS[i].pos.x;
        if (DrugATOMS[i].pos.y < rect[1][0])
            rect[1][0]=DrugATOMS[i].pos.y;
        else if (DrugATOMS[i].pos.y > rect[1][1])
            rect[1][1]=DrugATOMS[i].pos.y;
        if (DrugATOMS[i].pos.z < rect[2][0])
            rect[2][0]=DrugATOMS[i].pos.z;
        else if (DrugATOMS[i].pos.z > rect[2][1])
            rect[2][1]=DrugATOMS[i].pos.z;
    }

    for(i=0 ; i<3 ; i++) {  //
        //c[i] =(BOX[i][1]+BOX[i][0])/2.0; // 0: x, 1:y 2:z
        LigandCenter[i]=(rect[i][1]+rect[i][0])/2.0;
    }
    /* Adjust the center of Liagnd to Origin Point*/
    for (i=0;i<DrugAtomCount;i++) {
        DrugATOMS[i].pos.x=DrugATOMS[i].pos.x-LigandCenter[0];
        DrugATOMS[i].pos.y=DrugATOMS[i].pos.y-LigandCenter[1];
        DrugATOMS[i].pos.z=DrugATOMS[i].pos.z-LigandCenter[2];
        //printf("\n After(Ligand) (%8.3f %8.3f %8.3f)",LigandCenter[0],LigandCenter[1], LigandCenter[2]);
        //printf("\%8.3f %8.3f %8.3f",DrugATOMS[i].pos.x,DrugATOMS[i].pos.y, DrugATOMS[i].pos.z); getch();
    }

}


int CheckConnRep(int *ptr, int id, int own, int v)
{
    int i;
    if (own == v) return 1;
    for (i=0;i<id;i++)
        if (*(ptr+i)==v) return 1;
    return 0;
}
/*==========================================================
    Write Results as PDB Format
=================================================== */
void WritePDBFile()
{
    int i;
    char buf[100];
    FILE *OutF;

   sprintf(buf,"%s%s-%s-%s.pdb",PreStrPath,PDBName,LigandName,RunId);

    if( (OutF=fopen(buf,"w"))==NULL )
    {   /* open config file */
        printf("\n Tempfile open error ! ==> %s\n", buf) ;
        exit(1) ;
    }

    fprintf(OutF,
        "REMARK\n"
        "REMARK Protein : %-20.20s Res. Num : %4d  Receptor Atom Num : %4d \n"
        "REMARK Drug    : %-20.20s Atom Num : %4d  Single Bond Num   : %4d \n"
        //-------- Modification (Critical) : recorrect wrong word
        //"REMARK FitnessValue = HB&van deel Waal + Elect + IntraEnergy\n"
        "REMARK FitnessValue = HB&van der Waals + Elect + IntraEnergy\n"
        "REMARK %12.6f = %12.6f + %12.6f + %12.6f\n"
        "REMARK\n",
        PDBName, Protein.acidNum, ReceptorAtomCount,
        LigandName, DrugAtomCount, singleB[0].numSB,
        BEST_IND.value, BEST_IND.value-BEST_IND.singleE[0]-BEST_IND.CE,
        BEST_IND.CE, BEST_IND.singleE[0]
        );

    for (i=0;i<DrugAtomCount;i++)
    {   // write prediction drug structure
        char crood[100];

        strcpy(buf,OrgDrugSTR[i]);
        // Modification: add code
        // set residuce number of poses on HETATM to 999
        strncpy(buf+22," 999",4);        

        if(DrugATOMS[i].AcidType > 20 || DrugATOMS[i].AcidType < 0 )
            strncpy(buf+17,"PRE", 3);

        buf[21]= 'Z';  // assign Z chain for predicted chain

        sprintf(crood,"%8.3f%8.3f%8.3f",
                PredDrug3D[i].x+ProteinCenter[0],
                PredDrug3D[i].y+ProteinCenter[1],
                PredDrug3D[i].z+ProteinCenter[2]);
        strncpy(buf+30,crood, 24);
        fprintf(OutF,"%s",buf);
    }

    for(i=DrugAtomCount;i<OrgDrugSTRIndex;i++)
        fprintf(OutF,"%s",OrgDrugSTR[i]);

#ifdef FORCERECORD
    for(i=0;i<BestHBondCount;i++)
    {
        sprintf(buf,"DHBOND P%4d %3s %4d %c %c <-> D%4d %c %c %12.6f\n",
                DrugHBond[BEST][i].ResAcidSeqId,
                DrugHBond[BEST][i].ResAcidName,
                DrugHBond[BEST][i].ResAtomSeqId,
                DrugHBond[BEST][i].ResAtomName,
                DrugHBond[BEST][i].ResHbondType,
                DrugHBond[BEST][i].DrgAtomSeqId,
                DrugHBond[BEST][i].DrgAtomName,
                DrugHBond[BEST][i].DrgHbondType,
                DrugHBond[BEST][i].Energy
               );
        fprintf(OutF,"%s",buf);
    }

    for(i=0;i<BestElectCount;i++)
    {
        sprintf(buf,"DELECT P%4d %3s %4d %c <-> D%4d %c %12.6f\n",
                DrugElect[BEST][i].ResAcidSeqId,
                DrugElect[BEST][i].ResAcidName,
                DrugElect[BEST][i].ResAtomSeqId,
                DrugElect[BEST][i].ResAtomName,
                DrugElect[BEST][i].DrgAtomSeqId,
                DrugElect[BEST][i].DrgAtomName,
                DrugElect[BEST][i].Energy
               );
        fprintf(OutF,"%s",buf);
    }
#endif

    fclose(OutF);
}
void WritePDBFile2()
{
    int i;
    char buf[100];
    FILE *OutF;

    //sprintf(buf,"%s%s-%s-%s.pdb",PreStrPath,PDBName,LigandName,RunId);
    sprintf(buf,"%s%s-%s-%s.pdb",PreStrPath2.c_str(),PDBName,LigandName,RunId);

    //-------- Display Process --------
    //printf("Write docking ligand conformation to %s\n",buf);
    cout << endl; cout.flush(); 
    cout << "Write docking ligand conformation to " << buf << endl;
    cout << endl; cout.flush(); 
    
    ostringstream oss;
    if( (OutF=fopen(buf,"w"))==NULL )
    {   /* open config file */
        //printf("\n Tempfile open error ! ==> %s\n", buf) ;
        oss.str("");
        oss << "Predicted Pose" << buf << " open error" << endl;
        UniversalErrorMessageLogger(oss.str());
        exit(1) ;
    }

    fprintf(OutF,
        "REMARK\n"
        "REMARK Protein : %-20.20s Res. Num : %4d  Receptor Atom Num : %4d \n"
        "REMARK Drug    : %-20.20s Atom Num : %4d  Single Bond Num   : %4d \n"
        //-------- Modification (Critical) : recorrect wrong word
        //"REMARK FitnessValue = HB&van deel Waal + Elect + IntraEnergy\n"
        "REMARK FitnessValue = HB&van der Waals + Elect + IntraEnergy\n"
        "REMARK %12.6f = %12.6f + %12.6f + %12.6f\n"
        "REMARK\n",
        PDBName, Protein.acidNum, ReceptorAtomCount,
        LigandName, DrugAtomCount, singleB[0].numSB,
        BEST_IND.value, BEST_IND.value-BEST_IND.singleE[0]-BEST_IND.CE,
        BEST_IND.CE, BEST_IND.singleE[0]
        );

    for (i=0;i<DrugAtomCount;i++)
    {   // write prediction drug structure
        char crood[100];

        strcpy(buf,OrgDrugSTR[i]);
        // Modification: add code
        // set residuce number of poses on HETATM to 999
        strncpy(buf+22," 999",4);        

        if(DrugATOMS[i].AcidType > 20 || DrugATOMS[i].AcidType < 0 )
            strncpy(buf+17,"PRE", 3);

        buf[21]= 'Z';  // assign Z chain for predicted chain

        sprintf(crood,"%8.3f%8.3f%8.3f",
                PredDrug3D[i].x+ProteinCenter[0],
                PredDrug3D[i].y+ProteinCenter[1],
                PredDrug3D[i].z+ProteinCenter[2]);
        strncpy(buf+30,crood, 24);
        fprintf(OutF,"%s",buf);
    }

    for(i=DrugAtomCount;i<OrgDrugSTRIndex;i++)
        fprintf(OutF,"%s",OrgDrugSTR[i]);

#ifdef FORCERECORD
    for(i=0;i<BestHBondCount;i++)
    {
        sprintf(buf,"DHBOND P%4d %3s %4d %c %c <-> D%4d %c %c %12.6f\n",
                DrugHBond[BEST][i].ResAcidSeqId,
                DrugHBond[BEST][i].ResAcidName,
                DrugHBond[BEST][i].ResAtomSeqId,
                DrugHBond[BEST][i].ResAtomName,
                DrugHBond[BEST][i].ResHbondType,
                DrugHBond[BEST][i].DrgAtomSeqId,
                DrugHBond[BEST][i].DrgAtomName,
                DrugHBond[BEST][i].DrgHbondType,
                DrugHBond[BEST][i].Energy
               );
        fprintf(OutF,"%s",buf);
    }

    for(i=0;i<BestElectCount;i++)
    {
        sprintf(buf,"DELECT P%4d %3s %4d %c <-> D%4d %c %12.6f\n",
                DrugElect[BEST][i].ResAcidSeqId,
                DrugElect[BEST][i].ResAcidName,
                DrugElect[BEST][i].ResAtomSeqId,
                DrugElect[BEST][i].ResAtomName,
                DrugElect[BEST][i].DrgAtomSeqId,
                DrugElect[BEST][i].DrgAtomName,
                DrugElect[BEST][i].Energy
               );
        fprintf(OutF,"%s",buf);
    }
#endif

    fclose(OutF);
}



// -------- Modification (Optional) clear all data reads by ReadPDB --------
void ClearPoint(Point* point, int SIZE)
{
    for(int i=0;i<SIZE;i++) {
        point[i].x = 0.0; point[i].y = 0.0; point[i].z = 0.0;
    }
}
 
void ClearAtoms(Atoms* ATOMS, int MAX)
{
    for(int i=0;i<MAX;i++) {
        ATOMS[i].name = ' ';
        ATOMS[i].pos.x = 0.0;
        ATOMS[i].pos.y = 0.0;
        ATOMS[i].pos.z = 0.0;
        ATOMS[i].charge=0.0; 
        ATOMS[i].numB = 0;
        for( int j=0;j<6;j++) ATOMS[i].adj[j] = 0;
        ATOMS[i].SeqId     = 0;
        ATOMS[i].AtomSeqId = 0;
        ATOMS[i].AcidSeqId = 0;
        ATOMS[i].AcidType  = ' ';
        ATOMS[i].OtherName = NULL;
        ATOMS[i].HbType    = ' ';
        ATOMS[i]. EngType  = 0;
        ATOMS[i]. Weight   = 0.0;
        ATOMS[i]. WeiType  = ' ';
    }
}

void ClearProteinStr(ProteinStr& protein)
{
    protein.acidNum = 0;
    for(int i=0;i<MAX_ACID*MAX_ACID_ATOM;++i) {
        protein.AtomList[i].atomName = ' ';
        protein.AtomList[i].charge = 0.0;
        protein.AtomList[i].position.x = 0.0;
        protein.AtomList[i].position.y = 0.0;
        protein.AtomList[i].position.z = 0.0;
        protein.AtomList[i].EngType = 0;
        protein.AtomList[i].HbType = 0;
        protein.AtomList[i].acidIdx = 0;
        protein.AtomList[i].AcidSeqId = 0;
        protein.AtomList[i].AtomSeqId = 0;
        protein.AtomList[i].acidID = 0;
        protein.AtomList[i].flag = ' ';
    }
    for(int i=0;i<MAX_ACID;i++) {
        protein.AcidList[i].acidID = 0;
        protein.AcidList[i].flag = 0;
        protein.AcidList[i].AtomNum = 0;
    }
}

void ClearOrgDrugSTR()
{
    for(int i=0;i<MAXDRUGATOM*2;i++)
        for(int j=0;j<83;j++)
            OrgDrugSTR[i][j] = ' ';
}

void ClearPDB()
{
    ClearProteinStr(Protein);   // major 
    ClearOrgDrugSTR();
    ClearPoint(RefDrugOrgPosi,100);
    ClearPoint(PredDrug3D,MAXDRUGATOM);
    ClearAtoms(DrugATOMS,MAXDRUGATOM);
    ClearAtoms(DummyDrugATOMS,MAXDRUGATOM);
    ClearPoint(DrugOrgPosi,100);
    ClearAtoms(ReceptorATOMS,MAXATOM);
    ClearAtoms(DummyRecATOMS,MAXATOM);
    ClearAtoms(ProteinHET,MAXPROTEINHET);
}


//-------- ReadCavity --------

void ReadPDB2(const string& fileName)
{
    char line[256],fpath[256], *buf;
    FILE *inFile;
    int len, RRigFlag=0;

    sprintf(fpath,"%s%s",PDBPath,fileName.c_str());

    if((inFile=fopen(fpath,"r"))==NULL)    /* open config file */
    {
        //printf("Pdb file %s : Open Error !\n", fpath) ;
        sprintf(fpath,"%s",fileName.c_str());
        if((inFile=fopen(fpath,"r"))==NULL) 
        {
            fprintf(stderr,"Error: Pdb file %s : Open Error !\n", fpath);
            exit(0);
        }
    }

    //-------- Modification (Optional)--------
    ClearPDB();

    atomIndex = 0;          /* atom index in whole protein */
    acidCount = 0;          /* total residues of receptor protein  */
    atomCount = 0;          /* total atoms of a residue in receptor protein*/
    lastAcidSer       = -1;
    RefLigAtomCount   =  0;
    ProHETAtomCount   =  0; /* total atoms of hetatm (water,etc.) in receptor*/
    ProMetalAtomCount =  0;

    while(fgets(line, 90, inFile)!=NULL)
    {
        buf = substr(line,0,3);

        if (strcmp(buf, "TER")==0)  RRigFlag = 1;

        RecordAtom(line);

        if (RRigFlag == 1)
        {
            len= strlen(line);
            if (ConsiderHetATOM && len >= 82 && line[81] == 'P')
            {
                // read special atom (metal, wator...)
               // printf("\n%s",line); getchar();
                RecordProHet(line);
            }
            if (len >= 82 && line[81] != 'P' )
            {
                // read reference ligand hetatm for computing RMSD
                RecordRefLig(line);
            }
        }
        line[81]=' ';
    }

    Protein.AcidList[acidCount-1].AtomNum = atomCount;
    Protein.acidNum = acidCount;
    fclose(inFile);
}
