#include "fiolib.h" // tools for quick i/o and string processing
#include "fdock2.h" // global variables and functions
#include "fanaly.h" // modules for analyzing rlt file
#include "fuilib.h" // modules for interfacing with GUI
#include "freadent2.cpp"
#include "freadmol2.cpp"
#include "fscoring2.cpp"


/* 
FIOLIB
 */
namespace FAlgo {
    using namespace std;
    
    void FloydWarshallShortestPath( vector<vector<int> >& D ) 
    {
      size_t num_of_atom = D.size();
      int a, b , x, y ;
      for( size_t k=0 ; k<num_of_atom ; k++) 
        for( size_t i=0 ; i<num_of_atom ; i++) 
          for( size_t j=0 ; j<num_of_atom ; j++) 
          {
            if( j==i ) D[i][j] =0;
            else 
            {
                a = D[i][j] ;
                x = D[i][k] ;
                y = D[k][j] ;
                b = x+y ;
                if( a < b ) D[i][j] = a ;
                else        D[i][j] = b ;
            }           
          }
    }

}

namespace FioLib {
    using namespace std;

    void qdel(string job_names) 
    {
        string log_name=".qstat.log";
        string cmd = "qstat | grep " + job_names + " > " + log_name;
        system( cmd.c_str() );
        ifstream fin( log_name.c_str() );
        string s; vector<string> v;
        while(getline(fin,s))
        {
            s = s.substr(0,s.find_first_of('.'));
            v.push_back( s );
        }
        fin.close();
        if(v.size()==0) 
        {
            cout << "No " << job_names << " exist" << endl;
            return;
        }
        s = "qdel";
        cout << endl << s ;
        for( size_t i=0;i<v.size();i++)
        {
            cout << " " << v[i];
            s += ( " " + v[i]) ;
        }
        system( s.c_str() );
        cout << endl << endl;
    }
}

namespace STRING {
    using namespace std;

    string getPrefix(const string& s)
    {   
        int pos = s.find_last_of('.');
        if( pos > 0 )
            return s.substr(0,pos);
        else return s;
    }
    
    string getSuffix(const string& s)
    {   
        size_t pos = s.find_last_of('.');
        if( pos == string::npos )
            return s;
        return s.substr(pos);
    }
    
    string RemovePath(const string& s)
    {

        std::string::size_type pos = s.find('/');
        if( pos == std::string::npos )
        {  
            std::string::size_type pos2 = s.find('\\');
            if( pos2==std::string::npos ) {
                return s;
            } 
            return s.substr(s.find_last_of('\\')+1);
        }
        return s.substr(s.find_last_of('/')+1);
    
    }
    
    string getNoPathPrefix(const string& s)
    {
        return getPrefix(RemovePath(s));
    }
    
    template<typename T> string Num2Str(const T &number)
    {
        ostringstream oss;
        oss << number;
        return oss.str();
    }

    template<typename T> string itos(const T &number)
    {
        ostringstream oss;
        oss << setfill('0') << number;
        return oss.str();
    }
}



/* 
FDOCK2 
 */
ostringstream oss_profile,oss_fitness;
double ES_G_CONST, ES_L_CONST;  
bool isInPbsMode = false;
bool isUseImptAtom = false;
bool isExceedMAXDRUGATOM;
int RUNTIMES = 3;
bool GenInterProfile = false;
bool DOCKING = true; 
bool DrugList = false;
int  DrugFormat = 1;
char FIT_TYPE = '6';
double FIT_CHARGE = 1;
double FIT_HYDRO = 1; 
char FIT_BOX = '1';
char FIT_SB;
double InitSizeRate = 0.0125;  /* initial size of step size */
int ADASearchLen = 2;
int DECSearchLen = 2;
char OPER_DEC = '1';
char OPER_CAU = '1';
char OPER_DE = '0';
char OPER_GAU = '0';
char SerialNo[20] = "0001";
int MAXGEN = 70;   
int RefCrossDrug = 0;  
char LigandName[256];     
char PDBName[256];        
char ListName[256];
char RunId[10]; 
double LigBaseEnergy;          
int UnRotSBNum = 10;     
double RangeRatio;    
double TChargeEnergy; 
double LigandCenter[3]; 
double ProteinCenter[3];
double EngE[7]={ 0.20, 0.16, 0.15, 0.20, 0.20, 0.02, 0.20};
double EngR[7]={ 3.20, 3.50, 4.00, 4.00, 4.20, 2.00, 4.20};
double EngA[7]={ 2.85, 3.15, 3.55, 3.55, 3.75, 1.75, 3.75};
double EngAE[7]={ 0.20, 0.16, 0.15, 0.20, 0.20, 0.02, 0.20}; 
int KeyAtomCount=0;        
int PosiAtomCount=0;     
int HitKeyAtom;       
int HitPosiAtom;      
int RefLigAtomCount;
char RefLigAtom[MAXDRUGATOM][83];  
Point RefDrugOrgPosi[100];   
CircleType DrugCircle[MAXDRUGATOM];

char  OrgDrugSTR[MAXDRUGATOM*2][83]; 
int   OrgDrugSTRIndex = 0;
Point PredDrug3D[MAXDRUGATOM];     
Atoms DrugATOMS[MAXDRUGATOM];     
Atoms DummyDrugATOMS[MAXDRUGATOM]; 
Point DrugOrgPosi[100];      

Atoms ReceptorATOMS[MAXATOM];  
Atoms DummyRecATOMS[MAXATOM];   
int   ReceptorAtomCount;      
int   DrugAtomCount;         
int   LigHbondCount      = 0; 
int   LigChargeAtomCount = 0;    
int   HBindingCount      = 0;      
int   NonHBindingCount   = 0;  

int   DrugAtomConnNum[MAXDRUGATOM]; /* work with DrugAtomAdj[]: total atoms connecting to a atom within 3 links */
int   DrugAtomAdj[MAXDRUGATOM][200]; 

int ProHETAtomCount;            /* total atoms(for metal and water) of receptor */
int ProMetalAtomCount;
int ChargeContactCount;
Atoms ProteinHET[MAXPROTEINHET];
Molecule mole[2];
Molecule dummy[2];
SingleBond      singleB[2] ;
PeptideResidue PeptideDrug[30];
int PepResCount;   
ProteinStr Protein; 
int atomIndex;       
int acidCount;      
int atomCount;  
int lastAcidSer;
Model_Atom model[20][29];
int modelAtomNum[20];
Matrix tran_matrix, scale_matrix, rotate_matrix;
Chromosome  pop[NUMPOP];  
Chromosome  BEST_IND;  
double range[NUMP][3] ;  
double BOX[3][2];          
int    NO_VARS;         
char   chrfile[256];      
int    NO_POPSIZE;          
int    generation = 0;   
long   CountE     = 0;    
FILE   *PerFP;           
FILE   *fc;              
FILE   *rlogfp, *dlogfp;   
double curBest;            
double curWorst;          
double curAvg;            
double SBE[2];           
int    val[MAXATOM] ;
Point  Xaxis,Yaxis,Zaxis ;
int    hcount = 0; 

time_t InitTime;
time_t LastTime;
int LastGen;
double DEC_Expect_Imv,ADA_Expect_Imv,CAU_Expect_Imv;
double  DEC_Times_Imv,ADA_Times_Imv,CAU_Times_Imv;
double ADAGlobalMean,DECGlobalMean,REPLACE_RATE;
Chromosome nowInd,bestInd,temM_PInd;
Point p1, p2, p3;

ForceRecord DrugHBond[2][MAXDRUGATOM];
ForceRecord DrugElect[2][MAXDRUGATOM];

double ThresholdContactPairs;
int InitialPositionNum = 0;
double InitialPostionData[51][100];
vector<DockRecord> DResult;

string PreStrPath2= "./PrePDB/docked_Pose/";   // relocatable output directory
string  currentDrug;

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

Point Normalize(Point vec)
{
    double   v_length ;

    v_length=(double)sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z) ;
    vec.x/=v_length ;
    vec.y/=v_length ;
    vec.z/=v_length ;
    return(vec) ;
}

double VectorLength(Point p1,Point p2)
{
        double   length ;
        length=(p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)
                +(p1.z-p2.z)*(p1.z-p2.z) ;
        return(length) ;
}

int random(int val)
{
    return(rand()%val);
}

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


/* 
FIOLIB
 */

/*--------------------------------------------------------------------
 * I n i t i a l P o p u l a t i o n
     Randomly generate an initial population.
     size: population size
     arg : number of parameters (NO_VARS)
 *--------------------------------------------------------------------*/
void InitialPopulation(int size,int arg)
{
    int      i,j, idx=0 ;
    double   t=0.0,diff ;             /* random number */

    curAvg   = 0.0;
    curBest  = MAXIMUM;
    curWorst = MINIMUM;

    for(i=0;i<size;i++)
    {        /* generate each individual */
        t=0.0;
        //printf("\n");
        if (InitialPositionNum > 0 && i < size/2)
            idx = random(InitialPositionNum);    
        for(j=0;j<arg;j++) {/* generate each parameters */
            diff=fabs(range[j][1]-range[j][0]);
            if  (j<=5) { // axis position of drug
                if (InitialPositionNum > 0 && i < size/2) {
                     pop[i].para[j]=InitialPostionData[idx][j] ;
                     pop[i].var[j]=InitSizeRate*diff*0.1; // 0.04  //0.08
                    // printf("pop=%d idx=%d %d value=%f",i,idx,j,pop[i].para[j]); getchar();                   
                }
                else {
                     pop[i].para[j]=frandom(diff)+range[j][0] ;
                     pop[i].var[j]=InitSizeRate*diff; // 0.04  //0.08
                }
            }
            else {
                pop[i].var[j]=InitSizeRate*diff;  // 0.04 0.08
                pop[i].para[j]=frandom(diff)+range[j][0];
            }
            pop[i].cvar[j]=pop[i].var[j];
            pop[i].svar[j]=4.0*pop[i].var[j];
            t=t+pop[i].var[j];

        } /* for j */
        pop[i].Amean = t/arg;
        pop[i].value=FunctionEvaluation(&pop[i],i);
        /* find the best one and worst one */
        if(pop[i].value<curBest)
            curBest=pop[i].value ;
        if(pop[i].value>curWorst)
            curWorst=pop[i].value ;
        curAvg+=pop[i].value ;
    } /* for i : each individual */

    ES_G_CONST=1.0/sqrt(2.0*NO_VARS);
    ES_L_CONST=1.0/sqrt(2.0*NO_VARS*NO_VARS);
    ADAGlobalMean=t/arg;
    DECGlobalMean=4.0*pop[1].var[0];

    curAvg/=(double)size ;
    printf("\n\nInitial complete\n");
}

/* ===============================================
*  The operator of rotamer library
 =============================================== */
void TorAngleOper()
{
   int i,k;
  // for (j=0;j<ADASearchLen;j++) {
        for (i=0;i<NO_POPSIZE;i++) {
            for (k=6;k<NO_VARS;k++) {
                pop[i+NO_POPSIZE]=pop[i];
                if (singleB[0].type[k-6]==1) // sp3-sp3
                    pop[i+NO_POPSIZE].para[k]=SP3SP3Start+SP3SP3Angle*random(3);
                else if (singleB[0].type[k-6]==2)
                    pop[i+NO_POPSIZE].para[k]=SP2SP3Start+SP2SP3Angle*random(6);
            }
            pop[i+NO_POPSIZE].value=FunctionEvaluation(&pop[i+NO_POPSIZE],i);
        }
        for (i=0;i<NO_POPSIZE;i++)
           if (pop[i+NO_POPSIZE].value < pop[i].value)
                   pop[i]=pop[i+NO_POPSIZE];
  // }
}



/*----------------------------------------------------------------------
 * Initial globale variables
 *---------------------------------------------------------------------*/
void InitGlobalVar()
{
    double t[4][4] = {  {1.0,0.0,0.0,0.0},
                        {0.0,1.0,0.0,0.0},
                        {0.0,0.0,1.0,0.0},
                        {0.0,0.0,0.0,1.0} } ;

    CountE=0;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            tran_matrix.M[i][j]   = t[i][j];
            scale_matrix.M[i][j]  = t[i][j];
            rotate_matrix.M[i][j] = t[i][j];
        }
    }
    BEST_IND.value=MAXIMUM ;           /* no best population */
    Xaxis.x=Yaxis.y=Zaxis.z=1.0 ;
    Xaxis.y=Xaxis.z=Yaxis.x=Yaxis.z=Zaxis.x=Zaxis.y=0.0 ;
    SBE[0]=0;
    SBE[1]=0;

    REPLACE_RATE = 1.0;

    BestHBondCount = 0;
    BestElectCount = 0;
    CurtHBondCount = 0;
    CurtElectCount = 0;

}

/*------------------------------------------------------------------
    Generate the distribution of Cauchy
    Scale_value is the scale
--------------------------------------------------------------------*/
double Cauchy(double scale_val)
{
     double t1,val;
     t1 = frandom(1.0);
     val=tan((t1-0.5)*M_PI)*scale_val;
     return(val);
}

double FA04AS(double x)
{
        if (x > 0.0)
                return(frandom(1.0));
        else
                return(frandom(1.0) * 2.0 - 1.0);
}

/*
C***************
C
C   SUBROUTINE FA05A(ALPHA,BETA,Z) RETURNS A PSEUDO-RANDOM
C   NUMBER, Z, FROM THE NORMAL DISTRIBUTION WITH MEAN
C   ALPHA AND STANDARD DEVIATION BETA.
C
C   ALPHA,BETA AND Z ARE ALL SINGLE PRECISION.
C
C   THE H.S.L. ROUTINE FA04AS IS USED AS A SOURCE
C   OF UNIFORM PSEUDO-RANDOM NUMBERS ON (0,1).
C
C   THE METHOD USED IS DESCRIBED IN HARWELL
C   REPORT NUMBER CSS 89.
C
C***************
*/
double Gauss(double mean, double stdev)
{
    double PP[2]  = { 0.60776046, 0.87673085 };
    double T1[14] = {
            0.85776389,     -0.058610738,   0.43285607,
            0.0,            1.9827695,      1.0,
            0.0,            0.85776389,     1.0581054,
            -0.21549043,    4.4310323,      -4.4310323,
            -1.0,           -1.0 };
    double T2[5] = {
            0.85776389,     0.37424533,     0.46836965,
            0.77773058,     -0.095360829    };
    double V, U, BETAD, T, IV, ZERO = 0.0;
    int K;

    if ((V = FA04AS(-1)) >= 0.0)
        BETAD = fabs(stdev);
    else
        BETAD = -fabs(stdev);
/*      BETAD = fabs(stdev) * (V/fabs(V)); */
    V = fabs(V);
    if (V < PP[0]) goto x200;
    K = 1;
    if (V > PP[1])  K = 8;
x100:
    IV = 1;
    if (V < ZERO) IV = -1;
    V = FA04AS(IV);
    IV = 1;
    if (IV < ZERO) IV = -1;
    U = FA04AS(IV);
    if (V <= U) goto x110;
    T = V;
    V = U;
    U = T;
x110:
    V = T1[K-1]*V;
    U = T1[K]+T1[K+1]*U;
    if (V < T1[K+2]+T1[K+3]*U) goto x300;
    if (T1[K+4]*U <= T1[K+5]) goto x100;
    if (V*V > -4.0*U*U*log(U)) goto x100;
    goto x300;
x200:
    IV = 1;
    if (V < ZERO) IV = -1;
    V = T2[0]*FA04AS(IV);
    IV = 1;
    if (IV < ZERO) IV = -1;
    U = T2[1] + T2[2] * FA04AS(IV);
    if (V <= T2[3] + T2[4]*U) goto x300;
    if (V*V > -4.0*U*U*log(U)) goto x200;
x300:
    if ( U != 0)
        T = mean+BETAD*(V/U);
    else
        T = 1.0;
    if ( T > 4.0 ) T = 4.0;
    return(T);
}

double frandom(double n)
{
    n = fabs(n);
    return((double) n * (random(32767))/(double)(32767));
}


/*----------------------------------------------------------------------
 * C r o s s P r o d u c t
 *---------------------------------------------------------------------*/
Point CrossProduct(Point v1,Point v2)
{
    Point N ;

    N.x=v1.y*v2.z-v2.y*v1.z;
    N.y=v2.x*v1.z-v1.x*v2.z;
    N.z=v1.x*v2.y-v2.x*v1.y;
    return(N) ;
}

/*----------------------------------------------------------------------
 * I n n e r P r o d u c t
 *---------------------------------------------------------------------*/
double InnerProduct(Point *N,Point *L)
{
    double result;

    result=(N->x)*(L->x)+(N->y)*(L->y)+(N->z)*(L->z) ;
    return(result) ;
}

/*-------------------------------------------------------------------
 * M a t r i x M u l t V e c
 *-------------------------------------------------------------------*/
Point MatrixMultVec(Matrix mat,Point vec)
{
    Point ret;
    int i,j;
    double x[4],s[4] ;

    s[0]=vec.x;
    s[1]=vec.y;
    s[2]=vec.z;
    s[3]=1.0;

    for(i=0 ; i<4 ; i++)  {
        x[i]=0.0 ;
        for(j=0 ; j<4 ; j++)
            x[i]+=mat.M[i][j]*s[j] ;
    }
    ret.x=x[0]/x[3] ;
    ret.y=x[1]/x[3] ;
    ret.z=x[2]/x[3] ;
    return(ret) ;
}

/*-------------------------------------------------------------------
 * P o w e r
 *-----------------------------------------------------------------*/
double Power(double x, int n)
//double   x ;
//int     n ;
{
    int     i,j=1 ;
    double   value=1.0 ;

    i=n ;
    if(i) j=1 ; else j=0 ;
    while(i>1) {j<<=1 ; i>>=1 ;}
    while(j) {
        value*=value ;
        i=n/j ; n%=j  ;
        if(i) value*=x ;
        j>>=1 ;
    }
    return(value) ;
}

/*---------------------------------------------------------------------
 * U p d a t e S c a l e _ M
 *--------------------------------------------------------------------*/
Matrix UpdateScale_M(double scale)
//double   scale ;
{
    int i ;

    for(i=0 ;i<3 ;i++)
        scale_matrix.M[i][i]=scale ;
    return(scale_matrix) ;
}

/*---------------------------------------------------------------------
 * U p d a t e T r a n _ M
 *--------------------------------------------------------------------*/
Matrix UpdateTran_M(double x,double y,double z)
//double   x,y,z;
{
    tran_matrix.M[0][3]=x ;
    tran_matrix.M[1][3]=y ;
    tran_matrix.M[2][3]=z ;
    return(tran_matrix) ;
}

/*----------------------------------------------------------------------
 * U p d a t e R o t a t e _ M
 *****
 *      Arbitrarily rotation theta with respect to axis.
 *----------------------------------------------------------------------*/
Matrix UpdateRotate_M(Point axis,double theta)
//Point   axis ;
//double   theta ;
{
    double  c,s ;
    double   uxx,uyy,uzz,uxy,uyz,uxz ;
    double   uxS,uyS,uzS ;

    // sincos((double)theta,&s,&c) ;
    s=sin(theta);
    c=cos(theta);
    uxx=axis.x*axis.x ; uyy=axis.y*axis.y ;
    uzz=axis.z*axis.z ; uxy=axis.x*axis.y ;
    uyz=axis.y*axis.z ; uxz=axis.x*axis.z ;
    uxS=axis.x*s ; uyS=axis.y*s ; uzS=axis.z*s ;

    rotate_matrix.M[0][0]=uxx+c*(1.0-uxx) ;
    rotate_matrix.M[0][1]=uxy*(1.0-c)-uzS ;
    rotate_matrix.M[0][2]=uxz*(1.0-c)+uyS ;
    rotate_matrix.M[1][0]=uxy*(1.0-c)+uzS ;
    rotate_matrix.M[1][1]=uyy+c*(1.0-uyy) ;
    rotate_matrix.M[1][2]=uyz*(1.0-c)-uxS ;
    rotate_matrix.M[2][0]=uxz*(1.0-c)-uyS ;
    rotate_matrix.M[2][1]=uyz*(1.0-c)+uxS ;
    rotate_matrix.M[2][2]=uzz+c*(1.0-uzz) ;

    return(rotate_matrix) ;
}


Matrix MatrixMult(Matrix m1, Matrix m2)
{
    int     i,j,k;
    Matrix  res ;
    double t1;

    for (i=0;i<4;i++)
        for (j=0;j<4;j++) {
            t1=0.0;
            for (k=0;k<4;k++)
                t1=t1+m1.M[i][k]*m2.M[k][j];
            res.M[i][j]=t1;
        } // end of j
    return(res);
}

Matrix MatrixMultMatrix(Matrix m1, Matrix m2,Matrix m3)
{
    Matrix  res ;
    res=MatrixMult(m3,m2);
    res=MatrixMult(res,m1);
    return(res) ;
}


/*---------------------------------------------------------------------
 * T r a n s f o m a t i o n
 *--------------------------------------------------------------------*/
void Transformation(Molecule *p,Matrix mat)
//Molecule        *p ;
//Matrix          mat ;
{
    int i ;

    for(i=0 ;i<p->count ;i++)
        p->atom[i].pos=MatrixMultVec(mat,p->atom[i].pos) ;
}

/*--------------------------------------------------------------------
 * T r a v e r s e A n d M o d i f y
 *--------------------------------------------------------------------*/
void TraverseAndModify(Molecule *m,int start,Matrix mat)
//Molecule        *m ;
//int             start ;
//Matrix          mat ;
{
    int     t,index ;

    val[start]=1 ;    /* traversed */
    t=m->atom[start].numB ;
    while(t>0) {
        index=m->atom[start].adj[t-1] ;
        if(!val[index]) {
            m->atom[index].pos=MatrixMultVec(
                    mat,m->atom[index].pos) ;
            TraverseAndModify(m,index,mat) ;
        } /* end if : val */
        t-- ;
    } /* while : t */
}

/*----------------------------------------------------------------------
 * P a r t i a l T r a n s f o r m
 *---------------------------------------------------------------------*/
void PartialTransform(Molecule *m,int which,int index,double theta)
{
    int     i ;
    int     base,start ;
    Point   axis ;
    Matrix  mat ;
    //double   diff ;

    if(theta) {
        for(i=0 ; i<m->count ;i++)
            val[i]=0 ;              /* initial traversal array */
        base =singleB[which].atom[index][0] ;
        start=singleB[which].atom[index][1] ;
        axis.x=m->atom[start].pos.x-m->atom[base].pos.x ;
        axis.y=m->atom[start].pos.y-m->atom[base].pos.y ;
        axis.z=m->atom[start].pos.z-m->atom[base].pos.z ;
        axis=Normalize(axis) ;

        val[base]=1 ;   /* this way is forbidden */
        mat=MatrixMultMatrix(UpdateTran_M(-m->atom[base].pos.x,
                    -m->atom[base].pos.y,
                    -m->atom[base].pos.z),
                UpdateRotate_M(axis,theta),
                UpdateTran_M(m->atom[base].pos.x,
                    m->atom[base].pos.y,
                    m->atom[base].pos.z)
                ) ;
        TraverseAndModify(m,start,mat) ;
        /*diff=CalculateGroupE(m,base,start) ; */
    } /* end if :theta */
    /* return(diff) ;*/
}


/*--------------------------------------------------------------------------
 * I n i t R a n d o m
     Initial random number generator.
 *-------------------------------------------------------------------------*/
void InitRandom(void)
{
    time_t  tloc ;
    time(&tloc) ;
    srand(tloc) ;
}

double RadToDeg(double rad)
{
    double deg;
    deg=rad*180/M_PI;
    deg=deg-360*((int)(deg/360));
    if (deg < 0)
        deg=deg+360;
    return(deg);
}

/*--------------------------------------------------------------------
 * W r i t e C h r o m o s o m e
 *****
 *      Write out the best 5 chromosomes.
 *--------------------------------------------------------------------*/
void WriteChromosome()
//int     ending ;                        /* last one call or not */
{
    int j;
    /* write out the solution vector */
    //fprintf(PerFP, "\nEC=%7ld Var=%3d POP=%4d best=%17.10f time=%8ld\n",CountE,NO_VARS,NO_POPSIZE,BEST_IND.value,LastTime);

    //  if (_WEB_)
    //      return;

    /* write out the solution vector */
    for(j=0 ;j<NO_VARS ; j++) {
        if  ( j >= 3)
            fprintf(PerFP,"%f ",RadToDeg(BEST_IND.para[j])) ;
        else
            fprintf(PerFP,"%f ",BEST_IND.para[j]) ;

        /* write out funtion value */
        fprintf(PerFP,"%f %5d %f %f %f %4d %6ld\n",BEST_IND.CE,
                (int)(time(NULL)-InitTime),
                BEST_IND.value,
                BEST_IND.singleE[0],
                BEST_IND.rmsd,
                generation,CountE);
        fflush(PerFP);
    }
}

void ESSelection()
{
   int i,j,k;

  for (i=0;i<NO_POPSIZE;i++) {
       k=i;
       for (j=i+1;j<2*NO_POPSIZE;j++)
            if (pop[k].value>=pop[j].value) k=j;
        temM_PInd=pop[i];
        pop[i]=pop[k];
        pop[k]=temM_PInd;
   }
}

//--------------------------------------------------------------
// Interval Crossover
// interval: 1/2*fabs(p1-p2) + p1 - 1/2*fabs(p1-p2) + p2
//---------------------------------------------------------------
Chromosome BLX05CrossOver(Chromosome *p1,Chromosome *p2,char type)
{
   double MoveRate,diff,t;
   int i;

   temM_PInd=*p1;
   for (i=0;i<NO_VARS;i++) {
        MoveRate=frandom(1.5)-0.25;
        diff=MoveRate*(p2->para[i]-p1->para[i]);
        t=p1->para[i]+diff;
        if (t < range[i][1] && t > range[i][0])
               temM_PInd.para[i] = t;

       if (type == 'G')
               temM_PInd.var[i]=0.5*(p1->var[i]+p2->var[i]);
       if (type == 'C')
               temM_PInd.cvar[i]=0.5*(p1->cvar[i]+p2->cvar[i]);
       if  (type == 'D' )
           temM_PInd.svar[i]=0.5*(p1->svar[i]+p2->svar[i]);
   }
   return(temM_PInd);
}

//--------------------------------------------------------------
// ES recombination for varibales
// interval: 1/2*fabs(p1-p2) + p1 - 1/2*fabs(p1-p2) + p2
//---------------------------------------------------------------
Chromosome ESCrossOver(Chromosome *p1,Chromosome *p2, char type)
{
   int i;
   temM_PInd=*p1;
   for (i=0;i<NO_VARS;i++) {
           if (frandom(1.0)>=0.2)
                   temM_PInd.para[i]=p1->para[i];
           else
                   temM_PInd.para[i]=p2->para[i];
       if (type == 'G')
               temM_PInd.var[i]=0.5*(p1->var[i]+p2->var[i]);
       if (type == 'C')
               temM_PInd.cvar[i]=0.5*(p1->cvar[i]+p2->cvar[i]);
       if  (type == 'D' )
           temM_PInd.svar[i]=0.5*(p1->svar[i]+p2->svar[i]);
   }
   return(temM_PInd);
}

//--------------------------------------------------------------
// ES recombination for varibales
// interval: 1/2*fabs(p1-p2) + p1 - 1/2*fabs(p1-p2) + p2
//---------------------------------------------------------------
Chromosome IntervalCrossOver(Chromosome *p1,Chromosome *p2,char type)
{
   int i;
   //if (frandom(1.0)>=0.5) chvar=0;
   temM_PInd=*p1;
   for (i=0;i<NO_VARS;i++) {
        temM_PInd.para[i]=0.5*(p1->para[i]+p2->para[i]);

       if (type == 'G')
               temM_PInd.var[i]=0.5*(p1->var[i]+p2->var[i]);
       if (type == 'C')
               temM_PInd.cvar[i]=0.5*(p1->cvar[i]+p2->cvar[i]);
       if  (type == 'D' )
           temM_PInd.svar[i]=0.5*(p1->svar[i]+p2->svar[i]);
   }
   return(temM_PInd);
}

//------------------------------------------------------
// 1. original Chromosome is located
//          pop[NO_POPSIZE] - pop[NO_POPSIZE[2*NO_POPSIZE-1]
// 2. put best children into
//          pop[0] - pop[NO_POPSIZE[NO_POPSIZE-1]
//--------------------------------------------------------
void ES_Mutation(char type)
{
   int i,j;
   //double t1;

   ADA_Expect_Imv=0;
   ADA_Times_Imv=0;
   ADAGlobalMean=0;
   for (i=0;i<NO_POPSIZE;i++) {
       pop[i+NO_POPSIZE]=ES_Mutation_Oper(i,1.0);

       if (pop[i+NO_POPSIZE].value > pop[i].value) {
           pop[i].Amean=pop[i].Amean*ADA_DECRATE;
           for (j=0;j<NO_VARS;j++)
                  pop[i].var[j]=pop[i].var[j]*ADA_DECRATE;
       }

       if (type == 'R') {
           if (pop[i+NO_POPSIZE].value <= pop[i].value)
              pop[i]=pop[i+NO_POPSIZE];
       }
       ADAGlobalMean=ADAGlobalMean+pop[i].Amean;
  }
  ADAGlobalMean=ADAGlobalMean/(double)NO_POPSIZE;
  ADA_Expect_Imv=ADA_Expect_Imv/(double)(ADASearchLen*NO_POPSIZE);
  ADA_Times_Imv=ADA_Times_Imv/(double)(ADASearchLen*NO_POPSIZE);
}

//------------------------------------------------------
// 1. original Chromosome is located
//          pop[NO_POPSIZE] - pop[NO_POPSIZE[2*NO_POPSIZE-1]
// 2. put best children into
//          pop[0] - pop[NO_POPSIZE[NO_POPSIZE-1]
//--------------------------------------------------------
void Cau_Mutation(char type)
{
   int i,j;
   //double t1;

   CAU_Expect_Imv=0;
   CAU_Times_Imv=0;
   for (i=0;i<NO_POPSIZE;i++) {
       pop[i+NO_POPSIZE]=Cau_Mutation_Oper(i,1.0);
       if (pop[i+NO_POPSIZE].value > pop[i].value)
          for (j=0;j<NO_VARS;j++)
               pop[i].cvar[j]=pop[i].cvar[j]*ADA_DECRATE;

       if (type == 'R') {
            if (pop[i+NO_POPSIZE].value <= pop[i].value)
              pop[i]=pop[i+NO_POPSIZE];
       }
  }
  CAU_Expect_Imv=CAU_Expect_Imv/(double)(ADASearchLen*NO_POPSIZE);
  CAU_Times_Imv=CAU_Times_Imv/(double)(ADASearchLen*NO_POPSIZE);
}


//------------------------------------------------------
// 1. original Chromosome is located
//          pop[NO_POPSIZE] - pop[NO_POPSIZE[2*NO_POPSIZE-1]
// 2. put best children into
//          pop[0] - pop[NO_POPSIZE[NO_POPSIZE-1]
//--------------------------------------------------------
void DEC_Mutation(char type)
{
   int i,j;
   DEC_Expect_Imv=0;
   DEC_Times_Imv=0;
   DECGlobalMean=0;
   for (i=0;i<NO_POPSIZE;i++) {
       pop[i+NO_POPSIZE]=DEC_Mutation_Oper(i);
       if  (type == 'R')
            if (pop[i+NO_POPSIZE].value <= pop[i].value) {
               pop[i]=pop[i+NO_POPSIZE];
           }
    }
    for (i=0;i<NO_POPSIZE;i++) {  //
          for (j=0;j<NO_VARS;j++) {
               if (pop[i].svar[j] > VARIANCE_MIN ) {
                   pop[i].svar[j]=pop[i].svar[j]*ES_DECRATE;
                   DECGlobalMean=DECGlobalMean+pop[i].svar[j];
               }
         }
    }
    DECGlobalMean=DECGlobalMean/(double)(NO_POPSIZE*NO_VARS);
    DEC_Expect_Imv=DEC_Expect_Imv/(DECSearchLen*NO_POPSIZE);
    DEC_Times_Imv=DEC_Times_Imv/(DECSearchLen*NO_POPSIZE);
}


void FIT_Mutation(char type)
{
   int i;
   for (i=0;i<NO_POPSIZE;i++) {
       pop[i+NO_POPSIZE]=FIT_Mutation_Oper(i);
       if  (type == 'R')
            if (pop[i+NO_POPSIZE].value <= pop[i].value) {
               pop[i]=pop[i+NO_POPSIZE];
           }
    }
}

//------------------------------------------------------
// 1. original Chromosome is located
//          pop[NO_POPSIZE] - pop[NO_POPSIZE[2*NO_POPSIZE-1]
// 2. put best children into
//          pop[0] - pop[NO_POPSIZE[NO_POPSIZE-1]
//--------------------------------------------------------
void ROT_Mutation()
{
   int i;
   for (i=0;i<NO_POPSIZE;i++) {
       pop[i+NO_POPSIZE]=ROT_Mutation_Oper(i);
       if (pop[i+NO_POPSIZE].value <= pop[i].value)
            pop[i]=pop[i+NO_POPSIZE];
    }
}


//------------------------------------------------------
// ES-mutation for local search
//--------------------------------------------------------
Chromosome ES_Mutation_Oper(int idx, double Mul)
{
   double gv,v1,mean,imv,t;
   int i,j;

   bestInd.value=9e20;
   imv=0;
          //nowInd=pop[idx];
   for (i=0;i<ADASearchLen;i++) {
          mean=0;
          nowInd=pop[idx];
          if (frandom(1.0) < ADA_CROSSOVER_RATE) { // 0.1 for 24 , 0.2 for 28,30
              if (frandom(1.0)<=0.8) // 0.8 for 43
                 nowInd=ESCrossOver(&nowInd,&pop[random(NO_POPSIZE)],'G');
              else
                 nowInd=IntervalCrossOver(&nowInd,&pop[random(NO_POPSIZE)],'G');
          }
            gv = Gauss(0,1)*ES_G_CONST;
            for (j=0;j<NO_VARS;j++) {
                  v1=nowInd.var[j]*exp(gv+ES_L_CONST*Gauss(0,1));
                  t=nowInd.para[j]+v1*Mul*Gauss(0,1);
                  if (t > range[j][0] && t < range[j][1]) {
                          nowInd.para[j]=t;
                          nowInd.var[j]=v1;
                  }
                  mean+=nowInd.var[j];
           }
          nowInd.value=FunctionEvaluation(&nowInd,idx);
          if   (pop[idx].value-nowInd.value > 0)
               imv = imv+pop[idx].value-nowInd.value;
          mean=mean/NO_VARS;
          bestInd.Amean=mean;
          if (nowInd.value < pop[idx].value) {
              ADA_Expect_Imv = ADA_Expect_Imv + (pop[idx].value-nowInd.value);
              ADA_Times_Imv++;
              for (j=0;j<NO_VARS;j++)
                  if (nowInd.var[j] > 10.0 * nowInd.svar[j])
                      nowInd.svar[j]=0.1*nowInd.var[j];
          }
          if ( nowInd.value < bestInd.value)
                 bestInd=nowInd;
   }
   return(bestInd);
}

//------------------------------------------------------
// ES-mutation for local search
//--------------------------------------------------------
Chromosome Cau_Mutation_Oper(int idx, double Mul)
{
   double gv,v1,mean,imv,t;
   int i,j;

   bestInd.value=9e20;
   imv=0;
   for (i=0;i<ADASearchLen;i++) {
          mean=0;
          nowInd=pop[idx];
          if (frandom(1.0) < CAU_CROSSOVER_RATE) { // 0.1 for 24 , 0.2 for 28,30
              if (frandom(1.0)<=0.8) // 0.9 for 43
                 nowInd=ESCrossOver(&nowInd,&pop[random(NO_POPSIZE)],'C');
              else
                     nowInd=IntervalCrossOver(&nowInd,&pop[random(NO_POPSIZE)],'C');
          }
            gv = Gauss(0,1)*ES_G_CONST;
            for (j=0;j<NO_VARS;j++) {
                  v1=nowInd.cvar[j]*exp(gv+ES_L_CONST*Gauss(0,1));
                  t=nowInd.para[j]+v1*Cauchy(Mul*1.0);
                  if (t > range[j][0] && t < range[j][1]) {
                          nowInd.para[j]=t;
                          nowInd.cvar[j]=v1;
                  }
                  mean+=nowInd.var[j];
           }
          nowInd.value=FunctionEvaluation(&nowInd,idx);
          if   (pop[idx].value-nowInd.value > 0)
               imv = imv+pop[idx].value-nowInd.value;
          mean=mean/NO_VARS;
          if (nowInd.value < pop[idx].value) {
              CAU_Expect_Imv = CAU_Expect_Imv + (pop[idx].value-nowInd.value);
              CAU_Times_Imv++;
              for (j=0;j<NO_VARS;j++)
                  if (nowInd.cvar[j] > 100.0 * nowInd.svar[j])
                      nowInd.svar[j]=0.1*nowInd.cvar[j];
          }
          if ( nowInd.value < bestInd.value)
                 bestInd=nowInd;
   }
   return(bestInd);
}

//------------------------------------------------------
// ES-mutation for local search
//--------------------------------------------------------
Chromosome DEC_Mutation_Oper(int idx)
{
   double v1,imv=0,t;
   int i,j;

   bestInd.value=9e20;
   for (i=0;i<DECSearchLen;i++) {
       nowInd=pop[idx];
       if (frandom(1.0) < DEC_CROSSOVER_RATE)
       {  // 0.1 for 24  // 0.9 for 0.6 for 28,30
           if (frandom(1.0)<=0.8) // 0.9 for 43  // 0.1 for 6
               nowInd=ESCrossOver(&pop[idx],&pop[random(NO_POPSIZE)],'D');
           else
               nowInd=IntervalCrossOver(&pop[idx],&pop[random(NO_POPSIZE)],'D');
       }

       for (j=0;j<NO_VARS;j++) {
           v1=nowInd.svar[j]*Gauss(0,1);
           t=nowInd.para[j]+v1;
           if (t > range[j][0] && t < range[j][1])
               nowInd.para[j]=t;
       }
       nowInd.value=FunctionEvaluation(&nowInd,idx);
       if (pop[idx].value-nowInd.value > 0)
           imv = imv+pop[idx].value-nowInd.value;
       if (nowInd.value < pop[idx].value) {
           DEC_Expect_Imv = DEC_Expect_Imv + (pop[idx].value-nowInd.value);
           DEC_Times_Imv++;
       }
       if ( nowInd.value < bestInd.value)
           bestInd=nowInd;
   }
   return(bestInd);
}


//------------------------------------------------------
// ES-mutation for local search
//--------------------------------------------------------
Chromosome ROT_Mutation_Oper(int idx)
{
   double v1,t,t3=M_PI/3.0,t6=M_PI/6.0,gv;
   int i,j;

   bestInd.value=9e20;
   for (i=0;i<DECSearchLen;i++) {
          nowInd=pop[idx];
          gv = Gauss(0,1)*ES_G_CONST;
          for (j=0;j<NO_VARS;j++) {
              if (j < 6 || frandom(1.0) <= 0.5) {
                v1=nowInd.var[j]*exp(gv+ES_L_CONST*Gauss(0,1));
                t=nowInd.para[j]+v1*Gauss(0,1);
                if (t > range[j][0] && t < range[j][1]) {
                        nowInd.para[j]=t;
                        nowInd.var[j]=v1;
                }
              }
              else {
                    if (singleB[0].type[j-6]==1) // sp3-sp3
                        t=t3*random(3);
                    else if (singleB[0].type[j-6]==2) // sp3-sp2
                        t=t6*random(6);
                    else
                        t=0;
                    nowInd.para[j]=t;
              }
           }
          nowInd.value=FunctionEvaluation(&nowInd,idx);
          if ( nowInd.value < bestInd.value)
                 bestInd=nowInd;
   }
   return(bestInd);
}

void DE_LOOP()
{
   int i,j;
   for (j=0;j<ADASearchLen;j++) {
       for (i=0;i<NO_POPSIZE;i++) {
           pop[i+NO_POPSIZE]=DE_Oper(i);
       }
       for (i=0;i<NO_POPSIZE;i++)
           if (pop[i+NO_POPSIZE].value <= pop[i].value) {
                   pop[i]=pop[i+NO_POPSIZE];
           }
   }
}

//------------------------------------------------------
// Differential Evolution
//--------------------------------------------------------
#define CR 0.8
#define F  0.5
Chromosome DE_Oper(int idx)
{
   double t;
   //doublr t3=M_PI/3.0,t6=M_PI/6.0;
   int j,p1,p2,L,p;

//   F=frandom(0.5)+0.5;
//     CR=frandom(0.7)+0.3;
   p=random(NO_POPSIZE);
   while  (p==idx) {
       p=random(NO_POPSIZE);
   }
   p1=random(NO_POPSIZE);
   while  (p1==idx || p1==p)  {
       p1=random(NO_POPSIZE);
   }
   p2=random(NO_POPSIZE);
   while  (p2==idx || p2==p || p2==p1)  {
       p2=random(NO_POPSIZE);
   }
   nowInd=pop[idx];
   j= random(NO_VARS);
   L=0;
   do {
      //  if (L < 6 || frandom(1.0) <= 0.5) {
            t=nowInd.para[j]+F*(pop[p1].para[j]-pop[p2].para[j]);
            if (t > range[j][0] && t < range[j][1])
               nowInd.para[j] = t;
            j=(j+1) % NO_VARS;
         /*
        }
        else {
            if (singleB[0].type[L-6]==1) // sp3-sp3
                t=t3*random(3);
            else if (singleB[0].type[L-6]==2) // sp3-sp2
                t=t6*random(6);
            else
                t=0;
            nowInd.para[j]=t;
        }
        */
        L=L+1;
  } while ( frandom(1.0) < CR && L < NO_VARS);
    nowInd.value=FunctionEvaluation(&nowInd,idx);
    if (nowInd.value > pop[idx].value)
       nowInd=pop[idx];
    return(nowInd);
}

//------------------------------------------------------
// ES-mutation for local search
//--------------------------------------------------------
Chromosome FIT_Mutation_Oper(int idx)
{
   double v1,t;
   int i,j;

   bestInd.value=9e20;
   for (i=0;i<DECSearchLen;i++)
   {
       nowInd=pop[idx];
       if (frandom(1.0) < CAU_CROSSOVER_RATE)
       {  // 0.1 for 24  // 0.9 for 0.6 for 28,30

           if (frandom(1.0)<=0.8) // 0.9 for 43  // 0.1 for 6
               nowInd=ESCrossOver(&pop[idx],&pop[random(NO_POPSIZE)],'D');
           else if (frandom(1.0) <= 0.5) // 0.8 for 6,28, but worse for 30
               nowInd=IntervalCrossOver(&pop[idx],&pop[random(NO_POPSIZE)],'D');
           else
               nowInd=BLX05CrossOver(&pop[idx],&pop[random(NO_POPSIZE)],'D');
       }

       for (j=0;j<NO_VARS;j++)
       {
           v1=nowInd.svar[j]*Cauchy(1.0);
           t=nowInd.para[j]+v1;
           if (t > range[j][0] && t < range[j][1])
               nowInd.para[j]=t;
       }
       nowInd.value=FunctionEvaluation(&nowInd,idx);

       if ( nowInd.value < bestInd.value)
           bestInd=nowInd;
   }
   return(bestInd);
}

void StrDot(char *s)
{
    int i,j,L=0;
    for (i=0;s[i] && s[i]!= '.';i++);
    L=i;
    for (i=L-8,j=0; s[i] && j < 8; i++,j++) {
        s[j]=s[i];
    }
    s[8]='\0';
    //printf("\n %s", s); getch();
}

void showElecHbond()
{
    char buf[80];
    int i;
    for(i=0;i<BestHBondCount;i++)
    {
        sprintf(buf,"HBOND P%4d %3s %4d %c %c <-> D%4d %c %c %8.3f\n",
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
        printf("%s",buf);
    }
    for(i=0;i<BestElectCount;i++)
    {
        sprintf(buf,"ELECT P%4d %3s %4d %c <-> D%4d %c %12.6f\n",
                DrugElect[BEST][i].ResAcidSeqId,
                DrugElect[BEST][i].ResAcidName,
                DrugElect[BEST][i].ResAtomSeqId,
                DrugElect[BEST][i].ResAtomName,
                DrugElect[BEST][i].DrgAtomSeqId,
                DrugElect[BEST][i].DrgAtomName,
                DrugElect[BEST][i].Energy
               );
        printf("%s",buf);
    }
}

void GA_ES()
{
    //  if (REPLACE_RATE > 0.5)
    //       REPLACE_RATE = REPLACE_RATE * DEC_REPLACE_RATE;
if( isInPbsMode==false ) 
{
    if (DOCKING) {
        printf( "generation=%3d v=%2d EC=%6ld "
                //"Fit=%8.2f S0=%5.1f CE=%6.2f RM=%5.3f Ra=%4.3f t=%5ld\n",
                "Fit=%8.2f S0=%5.1f CE=%6.2f RM=%5.3f Ra=%4.3f t=%5d\n",
                generation, NO_VARS, CountE,
                BEST_IND.value,
                BEST_IND.singleE[0],
                BEST_IND.CE,
                BEST_IND.rmsd,
                BEST_IND.RangeRatio,
                (int)(time(NULL)-InitTime)
              );

    #if  0
        showElecHbond();
    #endif
    }
    //-------- Display Process --------
    else {
        printf("message: Current generation: %3d  complete\n",generation);
        //printf("Starting mutation and selection procedure...");
    }
    fflush(stdout);
}
    if (generation > 30 && frandom(1.0) > 0.5 )
    {
        if (OPER_DEC == '1' ) {
            DEC_Mutation('S');
            ESSelection();
        }
    }
    else
    {
        if (OPER_DEC == '1' ) DEC_Mutation('R');
    }

    if (OPER_CAU == '1') Cau_Mutation('R');
    if (OPER_DE  == '1') DE_LOOP();
    if (OPER_GAU == '1') ES_Mutation('R');

    //printf("complete\n");
}


void FCES()
{
    InitTime = time(NULL);
    int i;

    /* initial population P0'*/
    InitialPopulation(NO_POPSIZE,NO_VARS);
    
    //-------- Display Process --------
    printf("Applying FC_adaptive procedure\n\r");

    generation = 1;
#ifdef WEB
    for (i=0;i<MAXGEN;++i,++generation)
      Generating predicted   GA_ES();
#else
    if(!DOCKING) // screening
    {
        for(i=0;i<MAXGEN;++i,++generation)
            GA_ES();
    }
    else
    {
        WriteChromosome() ;
        for(int i=0;i<MAXGEN;++i,++generation)
        {
            if (generation%20==0)  WriteChromosome();
            GA_ES();
        }
        WriteChromosome() ;
    }
    // write docking ligand conformation
    
    ChangeConformation(&BEST_IND,BEST_IND.idx);
    for (i=0;i<DrugAtomCount;i++)
        PredDrug3D[i]=dummy[0].atom[i].pos;
    //-------- Modification : Use WritePDBFile2 --------
    WritePDBFile2();

#endif
}

void SetImptAtom()
{
    FILE *fp;
    char buf[100],WeiType=' ';


    //-------- Modification --------
    if( not isUseImptAtom ) return;
    if( (fp = fopen("imptatom.txt","r")) == NULL )
        return;
    KeyAtomCount=0;
    PosiAtomCount=0;
    while(fgets(buf, 80, fp) != NULL)
    {
        int AtomSeqId,i=0;
        float Weight;
        char *token;
        bool found = false;

        WeiType=' ';
        token = strtok(buf,  " \t");
        for (i=0;token == NULL && i<=20;i++) token = strtok(NULL, " \t");
        if(token == NULL || i>=20) break;
        AtomSeqId = atoi(token);

        token = strtok(NULL, " \t");
        for (i=0;token == NULL && i<=20;i++) token = strtok(NULL, " \t");
        if(token == NULL || i>=20) break;
        Weight    = atof(token);

        token = strtok(NULL, " \t");
        for (i=0;token == NULL && i<=20;i++) token = strtok(NULL, " \t");

        if(token != NULL) {
            WeiType   = token[0];
            if (WeiType < 'A' || WeiType > 'Z') WeiType = ' ';
        }

        for(int i=0;i<=ReceptorAtomCount;i++)
        {
            if(ReceptorATOMS[i].AtomSeqId == AtomSeqId)
            {
                ReceptorATOMS[i].Weight = Weight;
                ReceptorATOMS[i].WeiType = WeiType;
                if (WeiType == 'K') KeyAtomCount++;
                if (WeiType == 'P') PosiAtomCount++;
                found = true;
                break;
            }
        }
        #if 0
            if(!found)
                printf("Can't find the %d %f atom!\n",AtomSeqId,Weight);
            else
                printf("|%d |%f ||%c|| K=%d P=%d\n",AtomSeqId,Weight,WeiType,KeyAtomCount,PosiAtomCount);
            getchar();
        #endif
    }
    fclose(fp);
    return;
}


void ReadDrug()
{
    char path[100];
    FILE *fp;

    sprintf(path, "%s%s", DrugFormat > 1 ? DrgPath : LigPath, LigandName);
    std::cout << "LigandName is " << LigandName << std::endl;
    cutPath(LigandName);
    if( (fp = fopen(path, "r")) == NULL ) 
    {
        printf("Drug file %s open error!\n", path);
        exit(1);
    }
    
    switch(DrugFormat)
    {
        case 1:     ReadLigand(fp);  break;
        case 2:     ReadMol(fp);     break;
        case 3:     ReadMol2(fp);    break;
                    
        default:    ReadLigand(fp);  break;
    }
    fclose(fp);
    
    return;
}

bool ReadDrug2(const string& file)
{
    FILE *fp;
    char path[256];

    //DrugFormat = FioLib::getFormatID(file);
    //sprintf(path, "%s%s", DrugFormat > 1 ? DrgPath : LigPath, file.c_str());
    
    strcpy(LigandName,file.c_str());
    cutPath(LigandName);

    sprintf(path, "%s%s", DrgPath, file.c_str());
    if( (fp = fopen(path,"r")) == NULL ) {
        sprintf(path, "%s%s", LigPath, file.c_str());
        if( (fp = fopen(path,"r")) == NULL ) {
            sprintf(path, "%s", file.c_str());
            if( (fp = fopen(path,"r")) == NULL ) {
                //fprintf(stderr,"Error: Drug file %s : Open Error !\n", file.c_str());
                ostringstream oss;
                oss << "Error: Drug file " << file << " Open Error" << endl;
                UniversalErrorMessageLogger(oss.str());
                return false;
            }
        }
    }
    
    DrugFormat = FioLib::getFormatID(path);
    bool success;
    switch(DrugFormat)
    {
        case 1: success=ReadLigand2(fp);  break;
        case 2: success=ReadMol(fp);    break;
        case 3: success=ReadMol2(fp);   break;
        case -1: return false; 

        default: success=ReadLigand(fp); break;
    }
    fclose(fp);
    return success;
}

void SingleDocking() 
{
    ReadDrug();

#if 0
    for(int i=0;i<OrgDrugSTRIndex;i++)
        printf("%s",OrgDrugSTR[i]);
    return;
#endif
    
    DecideDrugCircle();  // decide circle-atom to judge the h-bond and charge

    DecideDrugAtomHBondType();
    DecideDrugConn();     //  within 3 links for computing intraEnergy
    SetSBRange();
    LigBaseEnergy = GetLigBaseEnergy();

#ifdef WEB  // for linking with PHP
    sprintf(chrfile,"%sID-%s.txt", PreStrPath, SerialNo);

    if ( (PerFP = fopen(chrfile,"w")) == NULL ) 
    {
        printf("File %s open error!\n",chrfile);
        exit(1);
    }
    fclose(PerFP);
#else
    sprintf(chrfile,"./PrePDB/log_Rlt/T-%s-%s-%d.rlt", PDBName, LigandName, NO_POPSIZE);

    if ( (PerFP = fopen(chrfile,"a")) == NULL ) 
    {
        printf("File %s open error!\n",chrfile);
        exit(1);
    }
#endif
    
    printf("\n\n"
            "Protein  : %-20.20s Atom Num : %4d  Receptor Atom Num : %4d \n"
            "%s   : %-20.20s Atom Num : %4d  Single Bond Num   : %4d \n"
            "X-ray Fit: %f Ratio: %8.5f CE: %f %f \n",
            PDBName,
            Protein.acidNum, ReceptorAtomCount,
            (DrugFormat > 1) ? "Drug  " : "Ligand",
            LigandName, DrugAtomCount, singleB[0].numSB,
            FunctionEvaluation(&nowInd, -1),
            RangeRatio, nowInd.CE, LigBaseEnergy);

    // ajust the centers of drug & Protein to original point
    AdjustReceptorPosi();       // adjust recpertor axis
    AdjustDrugPosi();
    InitRandom();

    for(int i=0;i<RUNTIMES;i++) /* test for same drug */
    {
            /* initial global variables */
#ifdef WEB
        sprintf(RunId, "%s", SerialNo);
        trimAllSpace(RunId);
        InitGlobalVar();
        FCES();

        unlink(chrfile);
        return 0;
#else
        Itoa(i, RunId, 10);

        trimAllSpace(RunId);
        InitGlobalVar();
        FCES();
#endif
    } // for
    
}


void DrugScreening()
{
    char buf[300];
    FILE *listfp;
    double bestfit = 9e20;  // for recording best run for each drug
    int    bestidx = 0;     // for recording best run for each drug
    DockRecord dr;          // recording and sorting screening results

    
    //char* tmp = substr(ListName,0,strlen(ListName)-4);

    if((listfp = fopen(ListName,"r")) == NULL ) 
    {
        printf("List file %s open error!\n", ListName);
        exit(1);
    }

    /* Log File:runtime log and docking result log */
    cutPath(ListName);
    sprintf(buf, "./PrePDB/screening_log/%s_run.log", ListName);
    if( (rlogfp = fopen(buf,"w")) == NULL )
    {
        printf("Runtime log file %s open error!\n", buf);
        exit(1);
    }

    sprintf(buf, "./PrePDB/screening_log/%s_dock.log", ListName);
    if( (dlogfp = fopen(buf,"w")) == NULL )
    {
        printf("Docking result log file %s open error!\n", buf);
        exit(1);
    }

    if (GenInterProfile) { //do not run docking: generate protein-ligand profile for statistic only
        fprintf(dlogfp,"LinandName ");
        for(int i=0 ; i<mole[1].count ;i++) {        /* protein */  
                sprintf(buf,"%4d%c-%3s",mole[1].atom[i].AtomSeqId, mole[1].atom[i].name,
                   (char*)( (mole[1].atom[i].AcidType == -1) ? mole[1].atom[i].OtherName :
                          //acidName[mole[1].atom[i].AcidType]));
                          acidName[(int)mole[1].atom[i].AcidType]));
                fprintf(dlogfp,"%s-V %s-H %s-E", buf,buf,buf);
                //printf("%s-V %s-H %s-E", buf,buf,buf);
         } 
         fprintf(dlogfp," AverConPair  ElecVal    TotalEnergy\n");  
         DrugFormat = 1;  
         // printf("atom=%d", mole[1].count);              
    }        
    /* ****  Read List File  ******/
    
    memset(buf ,0 , sizeof(buf));
    while (fgets(buf, 200, listfp) != NULL)
    {
        buf[strlen(buf)-1] = '\0';
        switch(buf[0])
        {
            case '#': // read comment
            {
                if(buf[1]=='!') // decide the drug fromat
                {
                    if     ( strncmp(buf+2,"PDB" ,3)==0 ) DrugFormat = 1;
                    else if( strncmp(buf+2,"MOL" ,3)==0 ) DrugFormat = 2;
                    else if( strncmp(buf+2,"MOL2",4)==0 ) DrugFormat = 3;
                    else 
                    {
                        printf("List format error!!\n");
                        exit(1);
                    }
                }
                continue;
            }
            case '\0':
            case ' ':
                continue;

            default  :
                break;
        } // switch
    
        /****  Open and Read Durg file or Ligand file   *****/
        sprintf(LigandName, "%s", buf); 
        ReadDrug();

        // decide circle-atom to judge the h-bond and charge
        DecideDrugCircle();  
        DecideDrugAtomHBondType();
        DecideDrugConn();         
        SetSBRange();
        
        if (GenInterProfile) { //do not run docking: generate protein-ligand profile for statistic only
            fprintf(dlogfp,"%s ",LigandName);
            FunctionEvaluation(&nowInd, -1);    
            continue;
        }        
        
#ifdef DEBUG
        printf("AFT: ||%s||_%s_||%d\n",buf,LigandName,DrugAtomCount);
#endif  
        //////////////////////////////
        //
        //  Start docking
        //
#ifndef DEBUG
        printf("\n\n"
               "Protein  : %-20.20s Atom Num : %4d  Receptor Atom Num : %4d \n"
               "%s   : %-20.20s Atom Num : %4d  Single Bond Num   : %4d \n",
               PDBName,
               Protein.acidNum, ReceptorAtomCount,
               (DrugFormat > 1) ? "Drug  " : "Ligand",
               LigandName, DrugAtomCount, singleB[0].numSB);

        // ajust the centers of drug & Protein to original point
        AdjustDrugPosi();
        InitRandom();
        
        bestfit = 9e20;
        for(int i=0;i<RUNTIMES;i++) /* test for same drug */
        {
            /* initial global variables */
            Itoa(i, RunId, 10);
            trimAllSpace(RunId);
            InitGlobalVar();
            FCES();
            //  Record and Write best result into runtime log
            if(bestfit > BEST_IND.value) 
            {
                bestfit = BEST_IND.value;
                bestidx = i;
            }
            sprintf(buf,"%2d %14f %s\n",i, BEST_IND.value, LigandName);
            fprintf(rlogfp,"%s",buf);

        } // for

        //////////////////////////////
        //
        //  Record best docking result for each drug
        //
        strcpy(dr.ProName,PDBName   );
        strcpy(dr.DrgName,LigandName);
        dr.AtomNum  = DrugAtomCount;
        dr.HbondNum = BestHBondCount;
        dr.ElectNum = BestElectCount;
        dr.fitness  = bestfit;
        dr.run      = bestidx;
        DResult.push_back(dr);
#endif

    } // end of while loop

    if (GenInterProfile) { // generate protein-ligand profile only
        fclose(rlogfp); 
        fclose(dlogfp);
        return;
    }        
        

    //////////////////////////////
    //
    //  Sorting and Write best result into docking result log
    //
    dr_qsort();

    sprintf(buf,"#\n# Top best docking result\n#\n"
            "%6s %-14.14s %-20.20s RUN Atom Hbond Elect\n",
            "# RANK", "FitnessValue",
            ((DrugFormat > 1) ? "DrugName":"LigandName"));
    fprintf(dlogfp,"%s",buf);

    //for(int i=0;i<DResult.size();i++) 
    for(int i=0;i<(int)DResult.size();i++) 
    {
        sprintf(buf,"%6d %14f %-20.20s %3d %4d %5d %5d\n",
                i+1,
                DResult[i].fitness,
                DResult[i].DrgName,
                DResult[i].run,
                DResult[i].AtomNum,
                DResult[i].HbondNum,
                DResult[i].ElectNum
               );
        fprintf(dlogfp,"%s",buf);
    }


    fclose(rlogfp); 
    fclose(dlogfp);

    return;
}


void usage(char *s)
{
    char p[] = "6 1 1 1 2 0 70 1 1 0 0";
    
    printf("Usage   : %s -[Options] POPSIZE PDBFILE DRUGFILE   \n"  ,s);
    printf("          %s -l [Options] POPSIZE PDBFILE LISTFILE \n\n",s);
    printf("Options : p -- DRUGFILE is PDB  format             \n"    );
    printf("          m -- DRUGFILE is MOL  format             \n"    );
    printf("          t -- DRUGFILE is MOL2 format             \n"    );
    printf("          l -- Using LISTFILE                      \n\n"  );
    printf("          s -- Drug Screening (only in list mode)  \n\n"  );
    printf("          h -- Show Usage                          \n\n"  );
    printf("Default : DRUGFILE is PDB format                   \n\n"  );
    printf("Example : %s 100 1hvr.pdb XK2                      \n"  ,s);
    printf("          %s -m 100 1hvr.pdb drug.mol              \n"  ,s);
    printf("Docking(PDB format):                               \n"    );
    printf("    fcdock 100 4dfr.pdb mtx.ent %s                 \n"  ,p);
    printf("Docking(MOL format):                               \n"    );
    printf("    fcdock -m 100 4dfr.pdb mtx.mol %s              \n"  ,p);
    printf("Screening:                                         \n"    );
    //printf("    fcdock -l -s 100 4dfr.pdb list.txt %s          \n"    );
    printf("    fcdock -l -s 100 4dfr.pdb list.txt             \n"    );
    exit(0);
}



void ParaParse(int argc, char *argv[])
{
    int  i;
    char Para[20];

    bool pstart = false;
    int  para   = 0;
    char *op, opt;

    for (i=1;i<argc;i++)
    {
        op = argv[i];

        //////////////////////////////
        //
        //  Control parameter (-m -l)
        //
        if (!pstart && *op == '-')
        {
            while( (opt = *(++op)) != '\0')
            {
                switch(toupper(opt))
                {
                    case 'P': DrugFormat  = 1;      break; // PDB  format
                    case 'M': DrugFormat  = 2;      break; // MOL  format
                    case 'T': DrugFormat  = 3;      break; // MOL2 format

                    case 'L': DrugList    = true;   break; // Multiple Drug
                    case 'S': DOCKING     = false;  break; // Screening
                    
                    case 'I': GenInterProfile=true; 
                        DOCKING   = false; DrugList    = true; break; // generate interaction profile(only for post-analysis)
                    
                    case 'A': isUseImptAtom=true;   break;  //-------- Modification (imptatom.txt) --------

                    case 'H':
                    case '?': usage(argv[0]);   break;
                    default : usage(argv[0]);   break;
                }
            }
            para++;
        } // if

        else if(!pstart)
        {
            NO_POPSIZE = atoi(argv[i]);         // population size
            std::cout << "NO_POPSIZE is: " << NO_POPSIZE << std::endl;
            sprintf(PDBName, "%s", argv[i+1]);  // receptor file

            if(DrugList) 
            {
                //  drug-ID List file for multipe-drug docking
                sprintf(ListName, "%s", argv[i+2]);
            }
            else// if !DrugList
            {
                // single drug file
                sprintf(LigandName, "%s", argv[i+2]);
            }

            pstart = true;
            memset(Para, 0, sizeof(Para));
            i += 2;
        }

        //////////////////////////////
        //
        //  Other parameter
        //
        else if(pstart)
        {
            switch(i+1-para)
            {
                case  5: FIT_TYPE   = op[0];    break; // get energy functype
                case  6: FIT_CHARGE = atof(op); break; // get the charge type
                case  7: FIT_HYDRO  = atof(op); break; // get hydrophbic type 1:(netral binding site) (<0 hydrophobic)
                case  8: FIT_SB     = op[0];    break; // consider the intro energy

                case  9: // get the search length
                {
                         ADASearchLen = atoi(op);
                         DECSearchLen = ADASearchLen;
                         break;
                }
                case 10: // get rotable rate of population
                {
                         int j = atoi(op);
                         j = ((j>=10) ? 10 : j);
                         UnRotSBNum=NO_POPSIZE*j/10;
                         break;
                }
                case 11: // set maximum enerations
                {
                         MAXGEN = atoi(op); // printf("MAXG=%d",MAXGEN);getchar();
                         if (MAXGEN < 20) MAXGEN=20;
                         break;
                }
                case 12: // set RUNTIMES, i.e., set the number of independent runs 
                         RUNTIMES   = atoi(op); 
                         if (RUNTIMES < 1) RUNTIMES=1; 
                         if (RUNTIMES > 20) RUNTIMES=20;   
                         break; 
                case 13: OPER_CAU   = op[0];    break; // Self-Ada Cauchy
                case 14: OPER_DE    = op[0];    break; // DE
                case 15: OPER_GAU   = op[0];    break; // Self-Gaussian
#ifdef WEB
                case 16:
                {
                        sprintf(SerialNo, "%s", op);
                        break;
                }
#else
                case 16: {
                        InitialPosition(op);
                        break;
                }                                  
#endif
            } // switch
            strcat(Para,op);
        } //else
       // printf("\n %d %s %d %d",i,op,para,i-para);
    } // for

#if 0
    printf("\n FIT_TYPE=%c FIT_CHARGE=%f FIT_SB=%c",FIT_TYPE,FIT_CHARGE,FIT_SB);
    printf("\n OPER_DEC=%c OPER_CAU=%c OPER_DE=%c OPER_GAU=%c",
          OPER_DEC, OPER_CAU, OPER_DE, OPER_GAU);
    printf("\n InitSizeRate=%f ADASearchLen=%d", InitSizeRate,ADASearchLen);
    getchar();
#endif
}

//int main(int argc,char *argv[])
int old_main(int argc,char *argv[])
{
    
    if(argc<4)              // parameters > 4
        usage(argv[0]);     // display correct parameter format
    
    ParaParse(argc,argv);   // parameter parsing

    ReadAMBERModel();
    ReadPDB(PDBName);       // using PDBPath
    DecideBindingAtoms();   // decide the binding atoms of receptor
    //printf("\n Bef: PDBN=%s",PDBName);
    cutPath(PDBName);
    //printf("\n After: PDBN=%s",PDBName);
    //ReadLigCharge();
    DecideBoundArea();      // decide the search ares
    SetImptAtom();          // set important atoms or intercations on receptor
    
    if ( !DOCKING)              // drug screening
    {
        if (! GenInterProfile) //do not run docking: generate protein-ligand profile for statistic only
            AdjustReceptorPosi();   // adjust recpertor axis
        DrugScreening();
    }
    else
        SingleDocking();
    
    return 0;
}


// --------------------------------------------
// Modification by Yu-Ju Chen 
// --------------------------------------------

void debug_arguments(GemdockInterface &g)
{
    printf("NO_POPSIZE=%d\n",NO_POPSIZE);
    printf("Cavity=%s\nLigand=%s\n",g.getCavityName().c_str(),g.getLigandName().c_str());
    printf("FIT_TYPE=%c\nFIT_CHARGE=%f\nFIT_HYDRO=%f\nFIT_SB=%c\n",FIT_TYPE,FIT_CHARGE,FIT_HYDRO,FIT_SB);
    printf("DECSearchLen=%d =ADASearchLen=%d\n\n", DECSearchLen,ADASearchLen);
    printf("UnRotSBNum=%d\n",UnRotSBNum);
    printf("MAXGEN=%d\nRUNTIMES=%d\n",MAXGEN,RUNTIMES);
    printf("OPER_CAU=%c OPER_DE=%c OPER_GAU=%c\n",
              OPER_CAU, OPER_DE, OPER_GAU);
}

void ReadPDBCavity(const string& fileName)
{
    sprintf(PDBName,"%s",STRING::getNoPathPrefix(fileName).c_str());
    ReadPDB2(fileName);
    DecideBindingAtoms();
    DecideBoundArea();
    SetImptAtom();
    if (! GenInterProfile)
        AdjustReceptorPosi();   // adjust recpertor axis
}

bool ReadDrugCompound(const string& file)
{
    currentDrug=file;
    bool ret=ReadDrug2(file);
    if( ret==false ) 
        return false;
    // decide circle-atom to judge the h-bond and charge
    DecideDrugCircle();  
    DecideDrugAtomHBondType();
    DecideDrugConn();         
    SetSBRange();
    return ret;
}

void UniversalErrorMessageLogger(const string& msg)
{
    cerr << msg ;
    cerr.flush();
}

int new_main(int argc,char *argv[])
{
    GemdockInterface gemdockInterface(argc,argv);
    gemdockInterface.runGemdock();

    /*
    OO gemdock usage:
        GemdockInterface g("cav1gwr_EST.pdb","MFCD00000222.mol2");
        g.setPopulation(100);
        g.setGeneration(70);
        g.setGARuntimes(3);
        g.setConfigFile("gemdock.gem");
        g.setOutputPath("outputPathByAPI");

        g.printParameters();
        //g.runFromConfigFile("gemdock.gem");
        g.runGemdock();

    */

    return 0;
}

bool contains(const string& str,const string& pattern) 
{
    string::size_type pos = str.rfind(pattern);
    return ( pos != string::npos ); // pos == string::npos means pattern not be found
}

int main(int argc,char *argv[])
{
    if( contains(argv[0],"bestdock") ) {
        cout << "Run gemdock in traditional mode" << endl;
        old_main(argc,argv);
    }
    else {
        cout << "Run gemdock in extended mode" << endl;
        new_main(argc,argv);
    }
}