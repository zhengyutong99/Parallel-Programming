#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stdarg.h>
#include <values.h>
#include <time.h>


//-------- interface ----------------
bool getRltAnalysis(const char *rlt_file, const char* output_file);






//-------- implementation details ----------------

int LineIdx;
FILE *outf, *inf;
void ReadRltFile(char *buf);
char PDBId[50];

char GetNum(char *line, char *tbuf);
void GetPDBId(char *buf);

/*

analyze line (string) and cut tbuf based on space such as
123 AA EE 12356:
1st: tbuf=123 type='N'
2nd: tubf=AA  type='S'
...
*/
char GetNum(char *line, char *tbuf)
{
    int i;
    char type='N',fnonsp=1;

    for (i=0; (fnonsp || line[LineIdx+i] != ' ') && line[LineIdx+i];i++) {
        tbuf[i]=line[LineIdx+i];
        if (tbuf[i] != ' ' && tbuf[i] != '.' && tbuf[i] != '-' )
            if (tbuf[i] > '9' || tbuf[i] < '0')
                type = 'C';
        if (fnonsp && tbuf[i] != ' ')
            fnonsp=0;
    }
    if (fnonsp)
        type='S';
    LineIdx=i+LineIdx+1;
    tbuf[i]='\0';
    return(type);
}

void GetPDBId(char *buf)
{
    int i,j,k=0;
    char tf[100];
    i=strlen(buf);
    for (j=0; j<i; j++) {
            k=0;
            j=j-1;
            while (buf[j] != '\\' && buf[j]!= '/' && buf[j] && j >= 0) {
                tf[k]=buf[j];
                k++; j--;
            }
            for (i=0,k=k-1;k>=0;i++,k--)
                PDBId[i]=tf[k];
            PDBId[i]='\0';
    }  
}

void ReadRltFile(char *buf)
//int     ending ;                        /* last one call or not */
{
    int Len, EndL, num=0, i, tidx, correct=0;
    // bRmsd, bFit are best value for the last generation
    //correct = counts of brmsd < 2.0
    double bRmsd=1e20, aRmsd, sRmsd=0, bFit=1e20,aFit,sFit=0,stime=0, rmsdct=0;
    double anybRmsd=1e20, anybFit=1e20; // best in any generation
    double t=0, tGen=0;
    double tRmsd[2], tFit[2], time[2], rmsd2bfit[2][3];//rmsd2bfit[][0]: bfit, rmsd2bfit[][1]: rmsd at bfit, rmsd2bfit[][2]: gen
    double tval[100];
    FILE *inFile;
    char line[2001], tbuf[50],type, MF=0;

    i=strlen(buf)-1;
    while (1) {
        if ( buf[i]=='-' ) { --i; continue; }
        if ( (buf[i] <= 'Z' && buf[i] >= 'A') ||
              (buf[i] <= 'z' && buf[i] >= 'a') )
            break;
        i--;
    }
    buf[i+1]='\0';

    if((inFile=fopen(buf,"r"))==NULL) {     /* open config file */
        fprintf(outf,"open result error: ==%s== =%d=\n", buf,i) ;
        return;
    }
    GetPDBId(buf);

    while(fgets(line, 2000, inFile)!=NULL) {
        LineIdx=0;
        EndL=0; // count the number of segment
        tidx=0; //
        Len=strlen(line);
        while (LineIdx < Len && EndL < 100 ) {
            type=GetNum(line,tbuf);
            t=atof(tbuf);
            if ( strlen(tbuf) >= 2 && type == 'N' ) {
                tval[tidx]=t;
                tidx++;
            }
            else
                EndL=110;
            EndL++;
            //printf("\n %s %f endL=%d %d %c", tbuf, t, EndL, LineIdx, type);
            //getch();
        }
        tidx=tidx-5;
        i=0;//index for last 5 segments, i=0    1   2           3     4
            //                             time fit single-bond rmsd, generation
        while (tidx > 5 && i < 5) {
           t=tval[tidx+i]; // get value for each segments
           MF=1;
           //printf("\n%f %d %d \n", t, tidx, i);
           switch (i) {
               case 0: // time
                       time[0]=t;
                       break;
               case 1: // fitness
                       tFit[0]=t;
                       if (t < anybFit)
                          {anybFit=t;}
                       break;
               case 2: break; // single bond
               case 3: // rmsd
                       tRmsd[0]=t;
                       if (t < anybRmsd) 
                            anybRmsd=t;
                      break;
               case 4: // generation
                       if (t < tGen) {
                       		rmsd2bfit[0][2]=t;
                          if (tFit[1] < bFit){
                          	bFit=tFit[1];
                          	rmsd2bfit[1][0]=rmsd2bfit[0][0];//get bfit
                            rmsd2bfit[1][1]=rmsd2bfit[0][1];//get rmsd at bfit
                            rmsd2bfit[1][2]=rmsd2bfit[0][2];//get gen
                           }
                          if (tRmsd[1] < bRmsd) {
                              bRmsd=tRmsd[1];
                          }
                          
                          sRmsd=sRmsd+tRmsd[1];
                          sFit=sFit+tFit[1];
                          stime=stime+time[1];
                          
                          if (tRmsd[1] < 2.0) { 
                            correct++;
                            num++;
                          }
                          //printf("\n %f %f %d", bFit, bRmsd,num);
                          //getch();
                       }

                       tFit[1]= rmsd2bfit[0][0]= tFit[0];
                       tRmsd[1]= rmsd2bfit[0][1]=tRmsd[0];
                       time[1]= time[0];
                       tGen=rmsd2bfit[0][2]=t;
                       //printf("%12.4f %10.5f %12.4f %10.5f %2d  \n", tFit[1], tRmsd[1], tFit[0], tRmsd[0],  tGen);
                       break;
           }  // switch
           i++;
        } // while
    }
    if (tFit[1] < bFit) {
                bFit=tFit[1];
                rmsd2bfit[1][0]=rmsd2bfit[0][0];//get bfit
                rmsd2bfit[1][1]=rmsd2bfit[0][1];//get rmsd at bfit
                rmsd2bfit[1][2]=rmsd2bfit[0][2];//get gen
                //printf("==> %12.4f %10.5f %3d\n",rmsd2bfit[1][0],rmsd2bfit[1][1],rmsd2bfit[1][2]);
    }
    if (tRmsd[1] < bRmsd) {bRmsd=tRmsd[1]; }
    if (tRmsd[1] < 2.0) correct++;
    num++;
    sRmsd=sRmsd+tRmsd[1];
    sFit=sFit+tFit[1];
    stime=stime+time[1];
    if (MF) {
        aFit= sFit/num;
        aRmsd= sRmsd/num;
        stime= stime/num;
        rmsdct= correct/num;
        printf("====> %12.4f %10.5f %3f\n",rmsd2bfit[1][0],rmsd2bfit[1][1],rmsd2bfit[1][2]);
        fprintf(outf, "%s %12.4f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %2d %2d %10.5f %3.2f \n", PDBId, bFit, rmsd2bfit[1][1], bRmsd, aFit, aRmsd, anybRmsd, anybFit, correct, num, stime, rmsdct);
        //printf("%s %12.4f %10.5f %10.5f %10.5f %10.5f %10.5f %2d %2d %10.1f\n", PDBId, bFit, bRmsd,aFit, aRmsd, anybRmsd, anybFit, correct, num, stime);
//        rmsd2bfit[0][0]=rmsd2bfit[0][1]=rmsd2bfit[0][2]=rmsd2bfit[1][0]=rmsd2bfit[1][1]=rmsd2bfit[1][2]=1e20;
    }
    else {
        fprintf(outf, "%s ==> No results  \n", PDBId );
        printf(" %s ==> No results  \n", PDBId );
    }
    //getch();
    fclose(inFile);
}

bool getRltAnalysis(const char *argv1, const char* argv2)
{
    char buf[100];
    if((inf=fopen(argv1,"r"))==NULL) {
        fprintf(stderr,"File %s open input file error!\n",argv1); 
        return false;
    }
    if((outf=fopen(argv2,"w"))==NULL) {
        fprintf(stderr,"File %s open output file error!\n",argv2);
        return false;
    }
    fprintf(outf, "PDBID                      bFit       rmsd2bFit    bRmsd      aFit        aRmsd    anybRmsd  anybFit  correct Run time correct-rate\n");
    while(fgets(buf, 90, inf)!=NULL) {
        ReadRltFile(buf);
    }
    return true;
}


/*
int main(int argc,char *argv[])
{
    int i;
    char buf[100];

    if(argc<3) {
            printf("Usage :: inputfile outfile\n",argv[0]) ;
            printf("Input list of .rlt files, and assign output name\n",argv[0]) ;
            exit(0) ;
    }

    if((inf=fopen(argv[1],"r"))==NULL) {
             printf("File %s open input file error!\n",argv[1]); exit(1) ;
    }

    if((outf=fopen(argv[2],"w"))==NULL) {
             printf("File %s open output file error!\n",argv[2]); exit(1) ;
    }

    fprintf(outf, "\nPDBID                      bFit       rmsd2bFit    bRmsd      aFit        aRmsd    anybRmsd  anybFit  correct Run time correct-rate\n");
    while(fgets(buf, 90, inf)!=NULL) {
        ReadRltFile(buf);
    }
    return EXIT_SUCCESS;
}

*/
