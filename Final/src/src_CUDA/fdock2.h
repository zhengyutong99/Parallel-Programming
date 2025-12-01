/*
 *  function name ended with 2 means new function used by new_main function
 *  this is for compatiable between new and old version
 *  for example: 
        WritePDB() --> WritePDB2()
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
#include "fiolib.h" 
using namespace std;

#ifndef FDOCK2_HEADER
#define FDOCK2_HEADER

/* BEGIN: SPECIAL PARAMETERS  */
#define  WST        /* if define WST run in LINUX environment else Windows*/
#undef   WEB         /* run in web environment */ 
#define  RefineDockedConf  0 /* 1: assign initial docked ligand conformation for refinement*/ 
#define  DelTerminalSB   1   /* 1: not consider the terminal SB (freadent.cpp)*/
extern bool  GenInterProfile;  /* true: not run docking: give docked ligand conf. to generate protein-ligand interaction profile for post-analysis*/
#undef  HAVE_PROSTR  //  output Receptor with ligand
#define SHOWFORCE    //  display h-bond informations (GA_ES)
#define FORCERECORD  0 // 1: record (output) h-bond and electrostatic interaction to ligand conformation file
/* END: SPECIAL PARAMETERS*/

extern int RUNTIMES;                   /* number of independent runs */

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

extern bool DOCKING;  // false==> run screening (no .rlt but generate run.log and dock.log, dock.log records best run for each ligand)
extern bool DrugList; // false: single drug True: Multiple drugs and ID lised in the file
extern int  DrugFormat;     // 1:PDB 2. MOL format 3. Mol2 format

extern char FIT_TYPE;     /* 5: 95 6: 965 A: Amber R:rmsd */
extern double FIT_CHARGE;     /* 0: no 1: including charge, < 0: penalty for charged ligand */
extern double  FIT_HYDRO;       /* 1: neutral binding site, < 0: penalty for polar ligand */
extern char    FIT_BOX;     /* 0: no 1: including BOX penalty (important) */
extern char    FIT_SB;     /* 0: no 1: include single bond */
extern double  InitSizeRate;  /* initial size of step size */
extern int     ADASearchLen;
extern int     DECSearchLen;
extern char    OPER_DEC;
extern char    OPER_CAU;
extern char    OPER_DE;
extern char    OPER_GAU;
extern char    SerialNo[20];

#define NUMPOP          4001    /* maximum number of population */
extern int MAXGEN;      /* maximun generation */

#define ConsiderHetATOM 1   /* 1: consider the structure water or metal as part of protein */
#define ConsiderWater   0   /* 1: consider the structure water, metal as part of protein */
#define ConsiderHetal   1   /* 1: consider metal atoms as part of protein */
extern int RefCrossDrug;         /* 1: reference drug is the docking drug */

#define BoundedLen      2.5 /* the boundary of binding area (1 or 3 )*/
extern char LigandName[256];        /* (3 char) which ligand binding to protein */
extern char PDBName[256];           /* which protein (4 char) */
extern char ListName[256];
extern char RunId[10];             /* which run for generate predicated structure only */
#define BIND_DIST     1200  /* distance of binding atoms between protein and ligand */

#define MAXSB           150              /* maximal number of single bond */
#define MAXGRID         1               /* table lookup */
#define MAXATOM         90000            /* maximum receptor atom number */
#define MAXDRUGATOM     200             /* maximum drug atom number */
#define NUMP            70              /* maximum number of parameter */
#define MINIMUM         -1.0E32
#define MAXIMUM         1.0E32
extern double  LigBaseEnergy;          /* base energy of x-ray ligand */
extern int UnRotSBNum;      /* number of ind. don't change single bond */
extern double  RangeRatio;     /* ratio=(number of distance < 8)/all number of distance between ligand and protein */
extern double  TChargeEnergy;  /* electron energy */
extern double  LigandCenter[3]; /* Original center of Ligand*/
extern double  ProteinCenter[3]; /* original center of the Protein */

// O, N , C, S, P, H, metal
// energy parameters of AutoDock
//double EngR[6]={ 1.60, 1.75, 2.00, 2.00, 2.10, 1.00};
extern double EngE[7];
extern double EngR[7];
extern double EngA[7];
extern double EngAE[7]; 

#define SP3SP3Angle     2.0*(M_PI/3.0)
#define SP2SP3Angle     2.0*(M_PI/6.0)
#define SP3SP3Start     -2.0*M_PI/3.0
#define SP2SP3Start     -M_PI
typedef struct{
        double   x,y,z ;         /* 3D coordinate */
} Point ;

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

extern int   KeyAtomCount;        /* number of key atoms defined in  WeiType */
extern int   PosiAtomCount;        /* number of position atoms defined in  WeiType */
extern int   HitKeyAtom;          /* number of ligand conformation hit Key atoms */
extern int   HitPosiAtom;         /* number of ligand conformation hit Position atoms */
extern int   RefLigAtomCount;
extern char  RefLigAtom[MAXDRUGATOM][83];  /* store original(HETATM) in cavity PDB as reference*/
extern Point RefDrugOrgPosi[100];          /* save posi of referecne Drug in Cavity PDB */

typedef struct{
    int type;       /* circle type (4, 5, 6,..) */
        int cid;        /* belong to which circle */
} CircleType;

extern CircleType  DrugCircle[MAXDRUGATOM];      /* denote circle 1: circle atom */


extern char  OrgDrugSTR[MAXDRUGATOM*2][83];  /* store original pdb string (HETATM) */
extern int   OrgDrugSTRIndex;
extern Point PredDrug3D[MAXDRUGATOM];      /* store drug original 3D coordinate */
extern Atoms DrugATOMS[MAXDRUGATOM];       /* store ligand atom list */
extern Atoms DummyDrugATOMS[MAXDRUGATOM];  /* store temp Drug */
extern Point DrugOrgPosi[100];             /* save posi of the Drug */

extern Atoms ReceptorATOMS[MAXATOM];       /* store receptor atom list */
extern Atoms DummyRecATOMS[MAXATOM];       /* store temp receptor atom list */
extern int   ReceptorAtomCount;            /* total binding atoms of receptor */
extern int   DrugAtomCount;                /* total atoms of ligand molecular */
extern int   LigHbondCount;       /* total atoms of forming h-bond */
extern int   LigChargeAtomCount;       /* total atoms of formal charge >0 */
extern int   HBindingCount;       /* total h-bond binding*/
extern int   NonHBindingCount;       /* total nonh-bond binding*/
//float  HBindingCount=0;           /* total h-bond binding*/
//float  NonHBindingCount=0;        /* total nonh-bond binding*/

extern int   DrugAtomConnNum[MAXDRUGATOM]; /* work with DrugAtomAdj[]: total atoms connecting to a atom within 3 links */
extern int   DrugAtomAdj[MAXDRUGATOM][200]; /* link id within 3 links connecting to a atom for computing intraEnergy */


//== record the information of hetatm of receptor (such as metal, wator)
extern int ProHETAtomCount;            /* total atoms(for metal and water) of receptor */
extern int ProMetalAtomCount;
extern int ChargeContactCount;

#define MAXPROTEINHET     1000  /* maximum atoms(for metal and water) number */
extern Atoms ProteinHET[MAXPROTEINHET];    /* store ligand atom list */
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

extern Molecule mole[2];       /* molecular data mole[0]:drug mole[1]:protein */
                        /* link to arrays DrugATOMS and ReceptorATOMS */
extern Molecule dummy[2];      /* temp store: DummyDrug and TempATOMS */

/*==================================================================
    BEGIN: Single Bond structure for binding drug and receptor
=============================================================== */
typedef struct{
        int     numSB ;         /* total number of single bond */
        int     atom[MAXSB][2] ;/* index of two atoms */
        char    type[MAXSB];    /* single bond sp2-sp2, sp3-sp2, sp3-sp3 */
        double  baseE[MAXSB] ;  /* base energy of the two group */
} SingleBond ;

extern SingleBond      singleB[2] ;    /* single bond 0: ligand 1: receptor*/

/* ==== Peptide drug =====*/
typedef struct {
        int     AcidSer;        /* serial residue number of peptide drug */
        int     AtomNum;        /* atom count of a residue */
        char    AcidType;       /* Which Acid type: 20 types */
        int     SeqId[20];      /* decide the connection of atom */
        char    Atomname[20];   /* atom name, N, O, C, ... */
        int     Start;          /* offset of residue, N=0,CA=1 */
} PeptideResidue;
extern PeptideResidue PeptideDrug[30];
extern int PepResCount;               /* numer of residue of peptide drug */

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

extern ProteinStr Protein;        /* Store original receptor protein */


extern int    atomIndex;          /* atom index in whole protein */
extern int    acidCount;          /* total residues of receptor protein  */
extern int    atomCount;          /* total atoms of a residue in receptor protein*/
extern int    lastAcidSer;
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

extern Model_Atom model[20][29];
extern int modelAtomNum[20];
/* ===== End == AMBER Model ================================*/

/* =====BEGIN == Matrix Structure ================================*/
typedef struct{
        double M[4][4] ;         /* 4 by 4 matrix */
} Matrix,*MatrixPtr  ;

extern Matrix tran_matrix,scale_matrix,rotate_matrix;

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

extern Chromosome  pop[NUMPOP];    /* population */
extern Chromosome  BEST_IND;       /* best population up to now */

/* =====END == chromosome Structure ================================*/


#define func1(x1,x2,x3,x4)      ( (x4)=(x1)*(x1)+(x2)*(x2)+(x3)*(x3) )
extern double range[NUMP][3] ;     /* range of NO_VARS (0,1,2: x,y,z bound values) others: 0 to 2PI*/
extern double BOX[3][2];           /* the search area */

extern int    NO_VARS;             /* number of NO_VARS */
extern char   chrfile[256];        /* chromosome file */
extern int    NO_POPSIZE;          /* population size */
extern int    generation;      /* number of generation */
extern long   CountE;      /* total number of evaluation */
extern FILE   *PerFP;              /* record performance */
extern FILE   *fc;                 /* file descriptor for config file */
extern FILE   *rlogfp, *dlogfp;    /* record screening results*/
extern double curBest;             /* current generation best value */
extern double curWorst;            /* current generation worst value */
extern double curAvg;              /* current generation average value */
extern double SBE[2];              /* single bond energy */


extern int    val[MAXATOM] ;
extern Point  Xaxis,Yaxis,Zaxis ;
extern int    hcount;

//======================BEGIN: parameters of FCEA ==============================
extern time_t InitTime;
extern time_t LastTime;
extern int LastGen;
extern double DEC_Expect_Imv,ADA_Expect_Imv,CAU_Expect_Imv;
extern double  DEC_Times_Imv,ADA_Times_Imv,CAU_Times_Imv;
extern double ADAGlobalMean,DECGlobalMean,REPLACE_RATE;
extern double ES_G_CONST, ES_L_CONST;         /* constant of self-adaptive mutation */
#define DEC_REPLACE_RATE 0.995
#define INIT_ADA_STEPSIZE 0.05
#define DEC_CROSSOVER_RATE 0.5 //0.3
#define ADA_CROSSOVER_RATE 0.3 //0.2
#define CAU_CROSSOVER_RATE 0.2 //0.1
#define ES_DECRATE 0.95
#define ADA_DECRATE 0.95
#define VARIANCE_MIN 0.001

extern Chromosome nowInd,bestInd,temM_PInd;
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
extern Point Normalize(Point vec);
Point MatrixMultVec(Matrix mat,Point vec) ;
extern double VectorLength(Point p1,Point p2);

void Transformation(Molecule *p,Matrix mat);
void TraverseAndModify(Molecule *m,int start,Matrix mat);
extern double ThresholdContactPairs;
double GetThresholdContactPairs(int LNum);

void StrDot(char *s);
extern void Itoa(int val, char *RunId, int Bas);

/*---------
Give three points : return the angle (p2)
-----------------*/
extern Point p1, p2, p3;

extern double angle(Point p1, Point p2, Point p3);



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

extern vector<DockRecord> DResult;



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

extern ForceRecord DrugHBond[2][MAXDRUGATOM]; // record all h-bond */
extern ForceRecord DrugElect[2][MAXDRUGATOM];

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

extern int _part(int l,int h);
extern void _qsort(int l,int h);
extern void dr_qsort();
//////////////////////////////
//
//  cut dir path
//
//////////////////////////////
extern void cutPath(char *pathname);
extern char *substr(const char *str,int pos,int len);

extern int InitialPositionNum;
extern double InitialPostionData[51][100];
extern void InitialPosition(char *fileN);


// --------------------------------------------
// Modification by Yu-Ju Chen 
// --------------------------------------------

extern bool isUseImptAtom;          // used for reading imptatom.txt
extern bool isInPbsMode;
extern string  PreStrPath2;   // relocatable output directory
extern string  currentDrug;
extern bool    isExceedMAXDRUGATOM;
extern ostringstream oss_profile,oss_fitness;
const   string tab("\t");

// public functions
void ReadPDBCavity(const string& fileName);
bool ReadDrugCompound(const string& file);
extern void setDResult(const string&,const string&,int ,int ,int ,double ,int );
extern void writeDResult(FILE* fp);

// private functions ( do not use this function for read/write your file )
bool ReadLigand2(char *fileName);
void RecordHetAtm2(char*);
void ReadPDB2(const string& fileName);
void WritePDBFile2();


// All error message must print through UniversalErrorMessageLogger
void UniversalErrorMessageLogger(const string& msg);
#endif 
