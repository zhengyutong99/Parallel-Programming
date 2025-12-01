/*
 *  function name ended with 2 means new function used by new_main function
 *  this is for compatiable between new and old version
 *  for example: 
        WritePDB() --> WritePDB2()
        ���� 
        bool ReadLigand2(char *fileName);
        void ReadPDB2(const string& fileName);
        void WritePDBFile2();
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cctype>
#include <ctime>
#include <vector>
#include <string>
using namespace std;

/* BEGIN: SPECIAL PARAMETERS  */
#define  WST        /* if define WST run in LINUX environment else Windows*/
#undef   WEB         /* run in web environment */ 
#define  RefineDockedConf  0 /* 1: assign initial docked ligand conformation for refinement*/ 
#define  DelTerminalSB   1   /* 1: not consider the terminal SB (freadent.cpp)*/
bool  GenInterProfile=false;  /* true: not run docking: give docked ligand conf. to generate protein-ligand interaction profile for post-analysis*/
#undef  HAVE_PROSTR  //  output Receptor with ligand
#define SHOWFORCE    //  display h-bond informations (GA_ES)
#define FORCERECORD  1 // 1: record (output) h-bond and electrostatic interaction to ligand conformation file
/* END: SPECIAL PARAMETERS*/

int RUNTIMES = 3;                   /* number of independent runs */

#ifdef WST
    #define     PDBPath     "../WCavPDB/"
    #define     LigPath     "../Ligand/"
    #define     CavPath     "../WCavPDB/"
    #define     DrgPath     "../Drug/"
    #define     PreStrPath  "PrePDB/"
    #include <unistd.h>
#else
    #define     PDBPath     "..\\WCavPDB\\"
    #define     LigPath     "..\\Ligand\\"
    #define     CavPath     "..\\WCavPDB\\"
    #define     DrgPath     "..\\Drug\\"
    #define     PreStrPath  "PrePDB\\"
#endif

bool DOCKING    = true;  // false==> run screening (no .rlt but generate run.log and dock.log, dock.log records best run for each ligand)
bool DrugList   = false; // false: single drug True: Multiple drugs and ID lised in the file
int  DrugFormat = 1;     // 1:PDB 2. MOL format 3. Mol2 format

char    FIT_TYPE     = '6';     /* 5: 95 6: 965 A: Amber R:rmsd */
double  FIT_CHARGE   = 1;     /* 0: no 1: including charge, < 0: penalty for charged ligand */
double  FIT_HYDRO    = 1;       /* 1: neutral binding site, < 0: penalty for polar ligand */
char    FIT_BOX      = '1';     /* 0: no 1: including BOX penalty (important) */
char    FIT_SB       = '1';     /* 0: no 1: include single bond */
double  InitSizeRate = 0.0125;  /* initial size of step size */
int     ADASearchLen = 2;
int     DECSearchLen = 2;
char    OPER_DEC     = '1';
char    OPER_CAU     = '1';
char    OPER_DE      = '0';
char    OPER_GAU     = '0';
char    SerialNo[20] = "0001";

#define NUMPOP          4001    /* maximum number of population */
int  MAXGEN=70;      /* maximun generation */

#define ConsiderHetATOM 1   /* 1: consider the structure water or metal as part of protein */
#define ConsiderWater   0   /* 1: consider the structure water, metal as part of protein */
#define ConsiderHetal   1   /* 1: consider metal atoms as part of protein */
int RefCrossDrug=0;         /* 1: reference drug is the docking drug */

#define BoundedLen      2.5 /* the boundary of binding area (1 or 3 )*/
char LigandName[256];        /* (3 char) which ligand binding to protein */
char PDBName[256];           /* which protein (4 char) */
char ListName[256];
char RunId[10];             /* which run for generate predicated structure only */
#define BIND_DIST     1200  /* distance of binding atoms between protein and ligand */

#define MAXSB           150              /* maximal number of single bond */
#define MAXGRID         1               /* table lookup */
#define MAXATOM         90000            /* maximum receptor atom number */
#define MAXDRUGATOM     200             /* maximum drug atom number */
#define NUMP            70              /* maximum number of parameter */
#define MINIMUM         -1.0E32
#define MAXIMUM         1.0E32
double  LigBaseEnergy;          /* base energy of x-ray ligand */
int UnRotSBNum=10;      /* number of ind. don't change single bond */
double  RangeRatio;     /* ratio=(number of distance < 8)/all number of distance between ligand and protein */
double  TChargeEnergy;  /* electron energy */
double  LigandCenter[3]; /* Original center of Ligand*/
double  ProteinCenter[3]; /* original center of the Protein */

// O, N , C, S, P, H, metal
// energy parameters of AutoDock
//double EngR[6]={ 1.60, 1.75, 2.00, 2.00, 2.10, 1.00};
double EngE[7]={ 0.20, 0.16, 0.15, 0.20, 0.20, 0.02, 0.20};
double EngR[7]={ 3.20, 3.50, 4.00, 4.00, 4.20, 2.00, 4.20};
double EngA[7]={ 2.85, 3.15, 3.55, 3.55, 3.75, 1.75, 3.75}; // good
//double EngA[7]={ 3.20, 3.40, 3.8, 3.8, 4.0, 2.00, 4.0}; // bad
//double EngAE[7]={ 0.40, 0.32, 0.3, 0.40, 0.40, 0.04, 0.40}; bad
double EngAE[7]={ 0.20, 0.16, 0.15, 0.20, 0.20, 0.02, 0.20}; // good

#define SP3SP3Angle     2.0*(M_PI/3.0)
#define SP2SP3Angle     2.0*(M_PI/6.0)
#define SP3SP3Start     -2.0*M_PI/3.0
#define SP2SP3Start     -M_PI
typedef struct{
        double   x,y,z ;         /* 3D coordinate */
} Point ;

/*
    switch(curAtom->name) {
        case 'O':
        case 'I':
                  curAtom->EngType = 0; break;
        case 'N': curAtom->EngType = 1; break;
        case 'C': curAtom->EngType = 2; break;
        case 'S': curAtom->EngType = 3; break;
        case 'P': curAtom->EngType = 4; break;
        default:  curAtom->EngType = 5;
    }
*/
/*==================================================================
    BEGIN: ATOM structure for binding drug and receptor
=============================================================== */
typedef struct {
    char    name ;          /* atom name, N, O, C, ... */
    Point   pos ;           /* 3D coordinate */
    double  charge ;        /* eletronic charge */
    int     numB ;          /* number of bond belonging to me */
    int     adj[6] ;        /* adjacent index of atom */

    int     SeqId;          /* Sequence atom number in whole molecule */
    int     AtomSeqId;      /* Sequence atom number */
    int     AcidSeqId;      /* Sequence acid number in a protein */
    char    AcidType;       /* Which Acid type: 20 types */
    char   *OtherName;
    
    char    HbType;         /* Type of Hydrogen bond:energy function*/
    int     EngType;        /* vander War Parameter: energy function*/

    float    Weight;         /* Assign imporant (or weak) Atom Weight, default=1.0*/
    char     WeiType;        /* 'K': key atom, 'P': Position constraint, default=' ' */
} Atoms,*AtomPtr;

int   KeyAtomCount=0;        /* number of key atoms defined in  WeiType */
int   PosiAtomCount=0;        /* number of position atoms defined in  WeiType */
int   HitKeyAtom;          /* number of ligand conformation hit Key atoms */
int   HitPosiAtom;         /* number of ligand conformation hit Position atoms */
int   RefLigAtomCount;
char  RefLigAtom[MAXDRUGATOM][83];  /* store original(HETATM) in cavity PDB as reference*/
Point RefDrugOrgPosi[100];          /* save posi of referecne Drug in Cavity PDB */

typedef struct{
    int type;       /* circle type (4, 5, 6,..) */
        int cid;        /* belong to which circle */
} CircleType;

CircleType  DrugCircle[MAXDRUGATOM];      /* denote circle 1: circle atom */


char  OrgDrugSTR[MAXDRUGATOM*2][83];  /* store original pdb string (HETATM) */
int   OrgDrugSTRIndex = 0;
Point PredDrug3D[MAXDRUGATOM];      /* store drug original 3D coordinate */
Atoms DrugATOMS[MAXDRUGATOM];       /* store ligand atom list */
Atoms DummyDrugATOMS[MAXDRUGATOM];  /* store temp Drug */
Point DrugOrgPosi[100];             /* save posi of the Drug */



Atoms ReceptorATOMS[MAXATOM];       /* store receptor atom list */
Atoms DummyRecATOMS[MAXATOM];       /* store temp receptor atom list */
int   ReceptorAtomCount;            /* total binding atoms of receptor */
int   DrugAtomCount;                /* total atoms of ligand molecular */
int   LigHbondCount      = 0;       /* total atoms of forming h-bond */
int   LigChargeAtomCount = 0;       /* total atoms of formal charge >0 */
int   HBindingCount      = 0;       /* total h-bond binding*/
int   NonHBindingCount   = 0;       /* total nonh-bond binding*/
//float  HBindingCount=0;           /* total h-bond binding*/
//float  NonHBindingCount=0;        /* total nonh-bond binding*/

int   DrugAtomConnNum[MAXDRUGATOM]; /* work with DrugAtomAdj[]: total atoms connecting to a atom within 3 links */
int   DrugAtomAdj[MAXDRUGATOM][200]; /* link id within 3 links connecting to a atom for computing intraEnergy */


//== record the information of hetatm of receptor (such as metal, wator)
int ProHETAtomCount;            /* total atoms(for metal and water) of receptor */
int ProMetalAtomCount;
int ChargeContactCount;

#define MAXPROTEINHET     1000  /* maximum atoms(for metal and water) number */
Atoms ProteinHET[MAXPROTEINHET];    /* store ligand atom list */
#ifdef HAVE_PROSTR
char  OrgProHetSTR[MAXPROTEINHET][200];    /* store original pdb string (HETATM) */
#endif
void  RecordProHet(char * line);
//== END of record the information of hetatm of receptor (such as metal, wator)



typedef struct{
    AtomPtr   atom ;      /* atom's data list */
    int       count;      /* total number of atoms */
    double    box ;       /* maximal bounding box */
} Molecule;

Molecule mole[2];       /* molecular data mole[0]:drug mole[1]:protein */
                        /* link to arrays DrugATOMS and ReceptorATOMS */
Molecule dummy[2];      /* temp store: DummyDrug and TempATOMS */

/*==================================================================
    BEGIN: Single Bond structure for binding drug and receptor
=============================================================== */
typedef struct{
        int     numSB ;         /* total number of single bond */
        int     atom[MAXSB][2] ;/* index of two atoms */
        char    type[MAXSB];    /* single bond sp2-sp2, sp3-sp2, sp3-sp3 */
        double  baseE[MAXSB] ;  /* base energy of the two group */
} SingleBond ;

SingleBond      singleB[2] ;    /* single bond 0: ligand 1: receptor*/

/* ==== Peptide drug =====*/
typedef struct {
        int     AcidSer;        /* serial residue number of peptide drug */
        int     AtomNum;        /* atom count of a residue */
        char    AcidType;       /* Which Acid type: 20 types */
        int     SeqId[20];      /* decide the connection of atom */
        char    Atomname[20];   /* atom name, N, O, C, ... */
        int     Start;          /* offset of residue, N=0,CA=1 */
} PeptideResidue;
PeptideResidue PeptideDrug[30];
int PepResCount;               /* numer of residue of peptide drug */

/*==== BEGIN: Data structure for saving original receptor protein ======*/
const char acidName[26][4] = {  /* 20 amino acid */
    "GLY", "ALA", "VAL", "LEU", "ILE",
    "SER", "CYS", "THR", "MET", "PRO",
    "PHE", "TYR", "TRP", "HIS", "LYS",
    "ARG", "ASP", "GLU", "ASN", "GLN",
    "MG ", "ZN ", "MN ", "FE ", "CA ", "CU ",
};

typedef struct
{
    char   atomName;        // Atom Name
    double charge;          // atom charge
    Point  position;        // 3D position
    int    EngType;         // parameters for energy functions
    int    HbType;      // for h-bond and disulfide bond
    //char   UniN[5];         // unique name
    int    acidIdx;         // the array index of acid in Protein
    int    AcidSeqId;
    int    AtomSeqId;
    int    acidID;
    char   flag;            // only for protein to ReceptorATOMS
#ifdef HAVE_PROSTR
    char   orgStr[82];      // store original pdb string
#endif
} ProteinAtom;

typedef struct
{
    int acidID;             // amino acid id
    int AtomNum;            // number atoms
    char flag;              // only for protein to ReceptorATOMS
}  ProteinAcid;


#define MAX_ACID  25000      // max number of residues in a protein
#define MAX_ACID_ATOM 14    // MAX atom number of an acid excluding h atom

typedef struct
{
    ProteinAtom AtomList[MAX_ACID*MAX_ACID_ATOM];   // GLOBAL ATOM LIST
    ProteinAcid AcidList[MAX_ACID];                 // GLOBAL ACID LIST
    int acidNum;                                    // total residues of a protein
} ProteinStr;

ProteinStr Protein;        /* Store original receptor protein */


int    atomIndex;          /* atom index in whole protein */
int    acidCount;          /* total residues of receptor protein  */
int    atomCount;          /* total atoms of a residue in receptor protein*/
int    lastAcidSer;
/*==== END: Data structure for saving original receptor protein ========*/

/* ===== Begin AMBER Model ========================================*/
typedef struct
{
    int num;                // atom label number
    int nA, nB, nC;         // connect atom number
    double length;          // bond length num - nA
    double angle;           // angle between nB-nA-num
    double diAngle;         // dihedral-angle nC-nB-nA-num
    double charge;          // atom charge
    char treeType;          // The topological type
    char atomName[4 + 1];   // Atom Unique Name
    char atomType[2 + 1];   // Atom Type
    Point position;       // 3D position
}  Model_Atom;

Model_Atom model[20][29];
int modelAtomNum[20];
/* ===== End == AMBER Model ================================*/

/* =====BEGIN == Matrix Structure ================================*/
typedef struct{
        double M[4][4] ;         /* 4 by 4 matrix */
} Matrix,*MatrixPtr  ;
Matrix tran_matrix,scale_matrix,rotate_matrix;

/* =====BEGIN == chromosome Structure ================================*/
typedef struct {
    double   para[NUMP] ;   /* parameter list */
    double   var[NUMP],cvar[NUMP],svar[NUMP];  /* step of parameter */
    double   value ;        /* function evaluation */
    double   singleE[2] ;   /* energy contributed by single bonds */
    double   Amean;
    double   rmsd;          /* rmsd between predicted and original */
    double   CE;            /* electron energy caused by charge energy*/
    double   RangeRatio;    /* testing only */
    int      idx;
} Chromosome;

Chromosome  pop[NUMPOP];    /* population */
Chromosome  BEST_IND;       /* best population up to now */

/* =====END == chromosome Structure ================================*/


#define func1(x1,x2,x3,x4)      ( (x4)=(x1)*(x1)+(x2)*(x2)+(x3)*(x3) )
double range[NUMP][3] ;     /* range of NO_VARS (0,1,2: x,y,z bound values) others: 0 to 2PI*/
double BOX[3][2];           /* the search area */

int    NO_VARS;             /* number of NO_VARS */
char   chrfile[256];        /* chromosome file */
int    NO_POPSIZE;          /* population size */
int    generation = 0;      /* number of generation */
long   CountE     = 0;      /* total number of evaluation */
FILE   *PerFP;              /* record performance */
FILE   *fc;                 /* file descriptor for config file */
FILE   *rlogfp, *dlogfp;    /* record screening results*/
double curBest;             /* current generation best value */
double curWorst;            /* current generation worst value */
double curAvg;              /* current generation average value */
double SBE[2];              /* single bond energy */


int    val[MAXATOM] ;
Point  Xaxis,Yaxis,Zaxis ;
int    hcount = 0;

//======================BEGIN: parameters of FCEA ==============================
time_t InitTime;
time_t LastTime;
int LastGen;
double DEC_Expect_Imv,ADA_Expect_Imv,CAU_Expect_Imv;
double  DEC_Times_Imv,ADA_Times_Imv,CAU_Times_Imv;
double ADAGlobalMean,DECGlobalMean,REPLACE_RATE;
double ES_G_CONST,ES_L_CONST;         /* constant of self-adaptive mutation */
#define DEC_REPLACE_RATE 0.995
#define INIT_ADA_STEPSIZE 0.05
#define DEC_CROSSOVER_RATE 0.5 //0.3
#define ADA_CROSSOVER_RATE 0.3 //0.2
#define CAU_CROSSOVER_RATE 0.2 //0.1
#define ES_DECRATE 0.95
#define ADA_DECRATE 0.95
#define VARIANCE_MIN 0.001

Chromosome nowInd,bestInd,temM_PInd;
//======================End: parameters of FCEA ==============================

#ifndef DEGREE
# define DEGREE(x)  ((x) * 180.0 / M_PI)
# define RADIAN(x)  ((x) * M_PI / 180.0)
#endif

/* ===========================================================
BEGIN: Read PDB (Receptor and Drug) FILE and Molecular Configuration
==========================================================*/
void ReadAMBERModel(void);
double AtomCharge(char * uniqName, int curAcidId);
double AtomFormCharge(char * uniqName, int curAcidId);
int HBond(char * uniqName, int curAcidId);
int CatchAcidID(char *string ,int startPos);
void RecordAtom(char * line);
void RecordHetAtm(char * line);
void RecordRefLig(char *line);
void RecordConn(char * line);
bool ReadSingleBond(char *line);    // new return type
int GetEndAtomNum(int idx,char AType);      /* count the number of atom Atype ("O","N") in Leaf*/
int GetCircleAtomSer(int idx);          /* Get serial of atom idx */
int GetCircleAtomType(int idx, char AType);

int GetDrugIdx(int serId);
//void GetChainID(char * line);
void DecideDrugConn(void);
int CheckConnRep(int *ptr, int id, int own, int v);
void ReadPDB(char *fileName);
void ReadLigand(char *fileName);
void GenPeptideConn(int SeqId, int ConnId, char SBType);
void GenPeptideResidueConn(int idx, int curAcidId);
void DecideBindingAtoms(void);
//-------- Modification Adding--------
void SetSBRange();
void AdjustReceptorPosi();
void WritePDBFile(void);
//void ReadSingleBond(char *line);
void DecideBoundArea(void);
void DecideDrugAtomHBondType(void);
void AdjustDrugPosi(void);
void DecideDrugCircle(void);
int TestDrugCircle(int idx);
int TestCircle(int idx, int count, int aNum[]);
/* ===========================================================
End: Read PDB (Receptor and Drug) FILE and Molecular Configuration
==========================================================*/

/* ===========================================================
BEGIN: FCEA FUNCTIONS: Mutations, recombinations, selections, ...
==========================================================*/
void ES_Mutation(char type);
void Cau_Mutation(char type);
void DEC_Mutation(char type);
void FIT_Mutation(char type);
void TorAngleOper(void);
void DE_LOOP(void);
void ROT_Mutation(void);

Chromosome ROT_Mutation_Oper(int idx);
Chromosome ES_Mutation_Oper(int idx,double Mul);
Chromosome Cau_Mutation_Oper(int idx,double Mul);
Chromosome DEC_Mutation_Oper(int idx);
Chromosome FIT_Mutation_Oper(int idx);
Chromosome DE_Oper(int idx);

Chromosome BLX05CrossOver(Chromosome *p1,Chromosome *p2,char type);
Chromosome IntervalCrossOver(Chromosome *p1,Chromosome *p2,char type);
Chromosome ESCrossOver(Chromosome *p1,Chromosome *p2,char type);
void ESSelection(void);

double Cauchy(double scale_val);
double Gauss(double mean, double stdev);
double FA04AS(double x);
void GA_ES(void);
void FCES(void);
/* ===========================================================
END: FCEA FUNCTIONS: Mutations, recombinations, selections, ...
==========================================================*/

double CalculateEnergy(Molecule *m1, Molecule *m2);
double ComputeElecStatic(double q1, double q2, double r);
double ComputeVDW(Atoms *atom1, Atoms *atom2, double r);
double ComputeVDW95(Atoms *atom1, Atoms *atom2, double r);
double ComputeVDW96(Atoms *atom1, Atoms *atom2, double r);
double ComputeVDW965(Atoms *atom1, Atoms *atom2, double r, bool intra, bool *HbF);

void ChangeConformation(Chromosome *ind, int idx);
double FunctionEvaluation(Chromosome *ind, int idx);
double CalculateIntraEnergy(Molecule *m, Chromosome *ind);
double ComputeTorEnergy(double dang, char type);
double GetLigBaseEnergy(void);

void InitGlobalVar(void);
void InitialPopulation(int size,int arg);
void InitRandom(void);
double RMSD(Molecule *m);

void WriteChromosome(void);

double RadToDeg(double rad);
double OutsideBox(Point p);
double frandom(double n);
double Power(double x,int n);
char* trimSpace(char *s);
char* trimAllSpace(char *s);

void PartialTransform(Molecule *m,int which,int index,double theta);
Matrix UpdateScale_M(double scale);
Matrix UpdateTran_M(double x,double y,double z);
Matrix UpdateRotate_M(Point axis,double theta);
Matrix  MatrixMultMatrix(Matrix m1,Matrix m2, Matrix m3) ;
Matrix  MatrixMult(Matrix m1,Matrix m2) ;
Point CrossProduct(Point v1,Point v2);
double InnerProduct(Point *N,Point *L);
Point Normalize(Point vec);
Point MatrixMultVec(Matrix mat,Point vec) ;
double VectorLength(Point p1,Point p2);

void Transformation(Molecule *p,Matrix mat);
void TraverseAndModify(Molecule *m,int start,Matrix mat);
double ThresholdContactPairs;
double GetThresholdContactPairs(int LNum);

void StrDot(char *s);
void Itoa(int val, char *RunId, int Bas);
/*-----------------------------------------------------------------------
 * V e c t o r L e n g t h
 *----------------------------------------------------------------------*/
double VectorLength(Point p1,Point p2)
{
        double   length ;
        length=(p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)
                +(p1.z-p2.z)*(p1.z-p2.z) ;
        return(length) ;
}
/*---------
Give three points : return the angle (p2)
-----------------*/
Point p1, p2, p3;
double angle(Point p1, Point p2, Point p3)
{
    //int i; 
    double  acc, d1, d2;
    acc = 0;
    acc=(p2.x-p1.x)*(p2.x-p3.x)+(p2.y-p1.y)*(p2.y-p3.y)+
        (p2.z-p1.z)*(p2.z-p3.z);

    d1 = sqrt(VectorLength(p1, p2));
    d2 = sqrt(VectorLength(p3, p2));
    if (d1 !=0 || d2 != 0)
        acc = acc / (d1 * d2);
    else acc=1;
    if (acc > 1.0)
        acc = 1.0;
    else if (acc < -1.0)
        acc = -1.0;
    acc = acos (acc);
    return (acc);
}

/*--------------------------------------------------------------------
 * N o r m a l i z e
 *--------------------------------------------------------------------*/
Point Normalize(Point vec)
{
    double   v_length ;

    v_length=(double)sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z) ;
    vec.x/=v_length ;
    vec.y/=v_length ;
    vec.z/=v_length ;
    return(vec) ;
}

void Itoa(int val, char *RunId, int Bas)
{
    int i=0,j;
    char buf[30];
    for (i=0; val >= 10 ;i++)
    {
        buf[i]='0'+ val%Bas;
        val=val/10;
    }
    buf[i  ] = '0' + val;
    buf[i+1] = '\0';
    for (j=0;i>=0;i--,j++)
    {
        RunId[j]=buf[i];
        //printf("\n %c %c",RunId[j], buf[i]);
    }
    RunId[j]='\0';
}


int random(int val)
{
    return(rand()%val);
}



//////////////////////////////
//
//  Screening Result Record Structure for sorting  dock.log
//
//////////////////////////////

typedef struct
{
    char ProName[30];
    char DrgName[30];
    int  AtomNum;
    int  HbondNum;
    int  ElectNum;
    double fitness;
    int    run;
} DockRecord;

vector<DockRecord> DResult;

//-------- Modification --------
void setDResult(const string& pn,const string& dn,
                int a,int h,int e,double f,int r)
{        
    DockRecord dr;          // recording and sorting screening results
    strcpy(dr.ProName,pn.c_str() );
    strcpy(dr.DrgName,dn.c_str());
    dr.AtomNum  = a;
    dr.HbondNum = h;
    dr.ElectNum = e;
    dr.fitness  = f;
    dr.run      = r;
    DResult.insert(DResult.end(),dr);

    
}

void writeDResult(FILE* fp)
{
    fprintf(fp,"#\n# Top best docking result\n#\n"
            "%6s %-14.14s %-20.20s RUN Atom Hbond Elect\n",
            "# RANK", "FitnessValue",
            ((DrugFormat > 1) ? "DrugName":"LigandName"));
    for(int i=0;i<(int)DResult.size();i++) {
        fprintf(dlogfp,"%6d %14f %-20.20s %3d %4d %5d %5d\n",
                i+1,  DResult[i].fitness, DResult[i].DrgName,
                      DResult[i].run,     DResult[i].AtomNum,
                      DResult[i].HbondNum,DResult[i].ElectNum);
    }

}


#ifdef FORCERECORD
typedef struct
{
    int     ResAcidSeqId;
    char*   ResAcidName;
    int     ResAtomSeqId;
    char    ResAtomName;
    char    ResHbondType;

    int     DrgAtomSeqId;
    char    DrgAtomName;
    char    DrgHbondType;

    double  Energy;
    double  Distance;
} ForceRecord;

//
// 0: BEST  1: CURRENT
//
#define BEST 0
#define CURT 1

ForceRecord DrugHBond[2][MAXDRUGATOM]; // record all h-bond */
ForceRecord DrugElect[2][MAXDRUGATOM];

static int BestHBondCount=0;   // number of H-bond
static int BestElectCount=0;
static int CurtHBondCount=0;
static int CurtElectCount=0;

#endif

//vector<HBondRecord> HBResult;


//////////////////////////////
//
//  Qsort for DResult
//
//////////////////////////////

int _part(int l,int h)
{
    int i=l,j=h+1;
    while((i+1)!=j)
    {
        while(DResult[i+1].fitness<DResult[l].fitness)
        {
            i++;
            if((i+1)==j)  break;
        }
        if((i+1)!=j)
            swap(DResult[i+1],DResult[--j]);
    }
    swap(DResult[i],DResult[l]);

    return i;
}

void _qsort(int l,int h)
{
    if(l<h)
    {
        int m=_part(l,h);
        _qsort(l,m-1);
        _qsort(m+1,h);
    }
}

void dr_qsort()
{
    _qsort(0,DResult.size()-1);
}


//////////////////////////////
//
//  cut dir path
//
//////////////////////////////

void cutPath(char *pathname)
{
    char *filename, *token;
    int i;
//-------- Modification --------
    strcpy(pathname,STRING::getNoPathPrefix(pathname).c_str());
    return;

    token = pathname;
    do
    {
        filename = token;
        token = strstr(token,"/");
    }
    while(token++ != NULL);
    token=filename;
    do
    {
        filename = token;
        token = strstr(token,"\\");
    }
    while(token++ != NULL);
    for (i=0;*(filename+i)!=0 && *(filename+i)!='.';i++)
        *(pathname+i)=*(filename+i);
    *(pathname+i)='\0';    
    //printf("\n path=%s|| fN=%s",pathname,filename); getchar(); 
}


char *substr(const char *str,int pos,int len)
{
    static char dst[90];

    strncpy(dst,str+pos,len);
    dst[len] = '\0';

    return dst;
}

int InitialPositionNum=0;
double InitialPostionData[51][100];
void InitialPosition(char *fileN){
  FILE *fp;
  char string[256];
  char *data;
  int data_num=0,data_para;
  
  InitialPositionNum=0;
  if( (fp = fopen(fileN,"r")) == NULL )
        return;
   while((fgets(string,256,fp))!=NULL) {
       data_para=-1;
       data=strtok(string," ");
       while(data!=NULL) {
             data_para++;
             if (data_para <=5)
                 InitialPostionData[data_num][data_para]=atof(data);
             data=strtok(NULL," ");
             // printf("\n %8.3f",InitialPostionData[data_num][data_para]);
       }
       if (data_para>=5) // must have x,y,x, and x,y,x rotation angles
              data_num++;
       if (data_num>=50)
              break;
    }  
    InitialPositionNum=data_num;
    //printf("data_num=%d",InitialPositionNum);getchar();
    fclose(fp);
}


// --------------------------------------------
// Modification by Yu-Ju Chen 
// --------------------------------------------

bool    isUseImptAtom = false;          // used for reading imptatom.txt
bool    isInPbsMode = false;
string  PreStrPath2= "./PrePDB/docked_Pose/";   // relocatable output directory
string  currentDrug;
bool    isExceedMAXDRUGATOM;
ostringstream oss_profile,oss_fitness;
const   string tab("\t");

// public functions
void ReadPDBCavity(const string& fileName);
bool ReadDrugCompound(const string& file);
void setDResult(const string&,const string&,int ,int ,int ,double ,int );
void writeDResult(FILE* fp);

// private functions ( do not use this function for read/write your file )
bool ReadLigand2(char *fileName);
void RecordHetAtm2(char*);
void ReadPDB2(const string& fileName);
void WritePDBFile2();


// All error message must print through UniversalErrorMessageLogger
void UniversalErrorMessageLogger(const string& msg);

