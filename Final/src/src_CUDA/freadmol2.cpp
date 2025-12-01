//////////////////////////////////////////////////
//
//  Read information from Mol or RDF format file
//  into structure that needed by docking
//
//  Structure: (defined in fdock.h)
//
//    DrugATOMS   : record information of all atoms
//    DrugOrgPosi : record orginal pos of all atoms
//    singleB     : record single bond infromation
//    OrgDrugSTR  : transform information to PDB format
//                  and will be outputed latter
//
//  last modified by Ian Chiu  2003.05.27
//  last modified by YuJu Chen 2007.09.03 (keyword: Modification)
//
//////////////////////////////////////////////////
#define MAXLINE 256

int GetDrugIdx(int serId);


bool ReadMol(FILE *fp)
{
    /*
      aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv

      Where::
      aaa  = number of atoms (current max 255)*  [Generic]
      bbb  = number of bonds (current max 255)*  [Generic]
      lll  = number of atom lists (max 30)*   [Query]
      fff  = (obsolete)
      ccc  = chiral flag:  0=not chiral, 1=chiral  [Generic]
      sss  = number of stext entries  [ISIS/Desktop]
      xxx  = (obsolete)
      rrr  = (obsolete)
      ppp  = (obsolete)
      iii  = (obsolete)
      mmm  = number of lines of additional  properties,
             including the M  END line
      No longer supported, the default is set to 999. [Generic]
      vvvvvv  = Ctab version: .V2000. or .V3000.  [Generic]
    */

    Atoms *curAtom;
    char line[MAXLINE], *token;
    int  info[10], AtomNum, ConNum, SBondNum;
    int  DrugStart, i, j, StrNum;

    memset(DrugATOMS  , 0, sizeof(DrugATOMS)   );
    memset(DrugOrgPosi, 0, sizeof(DrugOrgPosi) );
    memset(singleB    , 0, sizeof(singleB)     );
    memset(OrgDrugSTR , 0, sizeof(OrgDrugSTR)  );

    singleB[0].numSB = 0;
    DrugAtomCount=0;
    
    //-------- Modification (Critical) : use new mechanism to recognize MOL header --------
    int _atom_num,_bond_num,_rec_1,_rec_2,_rec_3;   

    while( fgets(line,MAXLINE,fp) != NULL)
    {

        //////////////////////////////////////////////////
        //
        //  Only process the lines including
        //  Counts line, Atom block, Bond block
        //
        
        //-------- Modification (Critical) : use new mechanism to recognize MOL header --------
        if( strlen(line) > 16 ) {
            _atom_num = atoi( substr(line,0,3)  );
            _bond_num = atoi( substr(line,3,3)  );
            _rec_1    = atoi( substr(line,6,3)  );
            _rec_2    = atoi( substr(line,9,3)  );
            _rec_3    = atoi( substr(line,12,3) );
        } else {
            _atom_num = _bond_num = _rec_1 = _rec_2 = _rec_3 = -1;
        }

        if( _atom_num > 1 && _bond_num > 1 &&  _rec_1>= 0 &&  _rec_2 >= 0 &&  _rec_3>=0)
        //if( strncmp(line+2,"-ISIS-",6) == 0 ||
        //    strncmp(line+34,"V2000",5) == 0 ||
        //    strncmp(line+34,"V3000",5) == 0   )
        {

            //////////////////////////////
            //
            //  Process Counts line
            //
            if( strncmp(line+2,"-ISIS-",6) == 0 ) 
            {
                if(fgets(line,MAXLINE,fp) == NULL) continue;
                if(fgets(line,MAXLINE,fp) == NULL) continue;
            }    
            
            for(i=0;i<10;i++)
            {
                token = substr(line,3*i,3);
                info[i] = atoi(token);
            }

            AtomNum   = info[0];
            ConNum    = info[1];
            SBondNum  = 0;
            curAtom   = DrugATOMS;
            DrugStart = 9000;
            StrNum    = AtomNum;

            //////////////////////////////
            //
            //  Process Atom block
            //
            for(i=0;i<AtomNum;i++)
            {
                if(fgets(line,MAXLINE,fp) != NULL)
                {
                    char AtomName = *substr(line,31, 1);
                    char AtomName2 = *substr(line,32, 1);//2�Ӧr��Atom type�ثe�bgemdock���������P���Ttype, ������C�B�z(�r�����S���ŦXOINSP).   
                    //  Dont record and do nothing when read H
                    if(AtomName == 'H') continue;
                        
                    j = DrugAtomCount++;
                    
                    curAtom[j].pos.x     = atof(substr(line, 0,10));
                    curAtom[j].pos.y     = atof(substr(line,11,10));
                    curAtom[j].pos.z     = atof(substr(line,21,10));
                    curAtom[j].name      = AtomName;
                    curAtom[j].SeqId     = j+1;
                    curAtom[j].AtomSeqId = i+1+DrugStart; 
                    // for GetDrugIdx() so use original AtomSeqId with H
                    curAtom[j].AcidSeqId = 0;
                    curAtom[j].charge    = 0.0;
                    curAtom[j].numB      = 0;
                    
                    
                    ///////////////////////////////////////////////////////////////////////////////////
                    //  ���ݵL�P�_ (Ca,Mg,Na,Zn,Fe ) 
                    //  Se,Sn,Si �|�P�_�� S (���n) 
                    //  Cu,Cr.oh,Cr.th,Co.oh �|�P�_�� C (NO Charge)
                    //  Mo,Mn,K,Li,Al ���w�qtype (NO Charge) 
                    ///////////////////////////////////////////////////////////////////////////////////
                    switch(curAtom[j].name)
                    {                
                        case 'O':
                        case 'I':   curAtom[j].EngType = 0; break;
                        case 'N':   curAtom[j].EngType = 1; break;
                        case 'C':   curAtom[j].EngType = 2; break;
                        case 'S':   curAtom[j].EngType = 3; break;
                        case 'P':   curAtom[j].EngType = 4; break;
                        default :   curAtom[j].EngType = 5;
                    }

                    //  Record the original position
                    DrugOrgPosi[j].x = curAtom[j].pos.x;
                    DrugOrgPosi[j].y = curAtom[j].pos.y;
                    DrugOrgPosi[j].z = curAtom[j].pos.z;
                    
                    //////////////////////////////
                    //
                    //  Create PDB format text and record it (HETATM)
                    //
                    if(AtomName2==' '){//2�Ӧr��Atom type�ثe�bgemdock���������P���Ttype, ������C�B�z(�r�����S���ŦXOINSP). 
                      sprintf(line,
                      //0         1         2         3    3    4    5    6
                      //123456    34     7890123  78901    9    7    5    1
                       "HETATM%5d  %-3c UNK  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                      ,curAtom[j].AtomSeqId           //  Atom serial number
                      ,curAtom[j].name                //  Atom name
                      ,curAtom[j].SeqId               //  Residue sequence number
                      ,curAtom[j].pos.x               //  coordinates for X
                      ,curAtom[j].pos.y               //  coordinates for Y
                      ,curAtom[j].pos.z               //  coordinates for Z
                      ,1.0                            //  Occupancy
                      ,0.0                            //  Temperature factor
                      );
                    }
                    else{
                      sprintf(line,
                      //0         1         2         3    3    4    5    6
                      //123456    34     7890123  78901    9    7    5    1
                       "HETATM%5d %c%-2c  UNK  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                      ,curAtom[j].AtomSeqId           //  Atom serial number
                      ,curAtom[j].name                //  Atom name
                      ,AtomName2                      //  2�Ӧr��Atom type�ثe�bgemdock���������P���Ttype, ������C�B�z(�r�����S���ŦXOINSP). 
                      ,curAtom[j].SeqId               //  Residue sequence number
                      ,curAtom[j].pos.x               //  coordinates for X
                      ,curAtom[j].pos.y               //  coordinates for Y
                      ,curAtom[j].pos.z               //  coordinates for Z
                      ,1.0                            //  Occupancy
                      ,0.0                            //  Temperature factor
                      );    
                    }
                    strcpy(OrgDrugSTR[j],line);
                    //printf("\n %s",line); getchar();
                }
            }

#if 0
            printf("DrugAtomCount %d\n",DrugAtomCount);
            printf("Atom    Count %d\n",AtomNum);
            printf("Connect Count %d\n",ConNum);       
#endif

            //////////////////////////////
            //
            //  Process Bond block
            //
            
            //  Initial 
            for(i=0;i<DrugAtomCount;i++)
                for(j=0;j<6;j++)    
                    curAtom[i].adj[j] = -1;
            
            //  Read Bond Info
            for(i=0;i<ConNum;i++)
            {
                if(fgets(line,MAXLINE,fp) != NULL)
                {
                    int con1, con2, connum;
                    
                    con1   = atoi(substr(line,0,3))-1;
                    con2   = atoi(substr(line,4,3))-1;
                    connum = atoi(substr(line,7,3));
                   
                    if( (con1 = GetDrugIdx(con1 + DrugStart+1)) < 0 || 
                        (con2 = GetDrugIdx(con2 + DrugStart+1)) < 0 )
                        continue;
                        
                    if(curAtom[con1].numB>6 || curAtom[con2].numB>6)
                    {
                        fprintf(stderr,"Too much connection at Atom %d\n",
                                con1>con2 ? con1 : con2);
                        fflush(stderr);
                        return false ;
                    }
                    
                    curAtom[con1].adj[curAtom[con1].numB++] = con2;
                    curAtom[con2].adj[curAtom[con2].numB++] = con1;

                    if(connum == 1)
                    {
                        //////////////////////////////
                        //
                        //  Record singlebond information
                        //
                        SBondNum = singleB[0].numSB;
                        singleB[0].atom[SBondNum][0] = con1;
                        singleB[0].atom[SBondNum][1] = con2;
                        singleB[0].type[SBondNum]    = 3;
                        singleB[0].numSB++;
                    }
                }
            }

            //////////////////////////////
            //
            //  Record connect information using PDB format
            //
            for(i=0;i<DrugAtomCount;i++)
            {
                int  seqid, connum=0;
                char buf[MAXLINE];
                memset(buf, 0, sizeof(buf));

                // Create PDB format text and record it (CONECT)
                
                sprintf(line, "CONECT%5d", curAtom[i].AtomSeqId);
                
                for(j=0;j<6;j++)
                {
                    seqid = curAtom[curAtom[i].adj[j]].AtomSeqId;
                    
                    if(curAtom[i].adj[j] != -1)
                    {
                        sprintf(buf, "%5d", seqid);
                        strcat(line, buf);
                        connum++;
                    }
                }
                strcat(line, "\n");
#if 0        
                printf("\n %s",line); getchar();
#endif

                if(connum)
                    strcpy(OrgDrugSTR[StrNum++],line);
            }

            //////////////////////////////
            //
            //  Record single bond information using PDB format
            //
            for(i=0;i<singleB[0].numSB;i++)
            {
                sprintf(line, "SINGLE%5d%5d 2 2\n",
                        curAtom[singleB[0].atom[i][0]].AtomSeqId,
                        curAtom[singleB[0].atom[i][1]].AtomSeqId);

                strcpy(OrgDrugSTR[StrNum++],line);
            }

            OrgDrugSTRIndex = StrNum;
            mole[0].count = DrugAtomCount;  // total atoms of drug
            mole[0].atom  = DrugATOMS;
            dummy[0].atom = DummyDrugATOMS;

            //  Ok, read a drug information
            return true;
        }
    }
    //  No read
    return false;
}



bool ReadMol2(FILE *fp)
{
    /*
      # ... comments
      
      @<TRIPOS>MOLECULE
      mol_name
      num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
      mol_type
      charge_type

      @<TRIPOS>ATOM
      atom_id atom_name x y z atom_type [subst_id 
      [subst_name [charge [status_bit]]]]

      where:
       atom_id          int     : ID # of the atom
       atom_name        string  : name of the atom
       x                real    : x coordinate of the atom
       y                real    : y coordinate of the atom
       z                real    : z coordinate of the atom
       atom_type        string  : SYBYL atom type for the atom
       subst_id         int     : ID # of the substructure containing the atom
       subst_name       string  : name of the substructure containing the atom
       charge           real    : charge associated with the atom
       status_bit       string  : internal SYBYL status bits associated 
                                  with the atom
      
      @<TRIPOS>BOND
      bond_id origin_atom_id target_atom_id bond_type [status_bits]

      where:
       bond_id          int     : ID # of the bond
       origin_atom_id   int     : ID # of the atom at one end of the bond
       target_atom_id   int     : ID # of the atom at the other end of the bond
       bond_type        string  : SYBYL bond type
                                  1  = single 
                                  2  = double
                                  3  = triple
                                  am = amide
                                  ar = aromatic
                                  du = dummy
                                  un = unknown
                                  nc = not connected
       status_bits      string  : internal SYBYL status bits associated 
                                  with the bond
    */

    Atoms *curAtom;
    char line[MAXLINE], *token;
    int  AtomNum, ConNum, SBondNum;
    int  DrugStart, i, j, StrNum;

    memset(DrugATOMS  , 0, sizeof(DrugATOMS)   );
    memset(DrugOrgPosi, 0, sizeof(DrugOrgPosi) );
    memset(singleB    , 0, sizeof(singleB)     );
    memset(OrgDrugSTR , 0, sizeof(OrgDrugSTR)  );
    singleB[0].numSB = 0;
    // Modification -- critical -- clean Drug Atom count
    DrugAtomCount=0;
    
    
    //cout << "Read Drug format mol2 " << endl;

    while( fgets(line,MAXLINE,fp) != NULL)
    {

        //////////////////////////////////////////////////
        //
        //  Only process the lines begin with
        //  @<TRIPOS>
        //

        //  Do nothing for comments
        if( strncmp(line,"#",1) == 0)
        {
            continue;
        }

        else if( strncmp(line,"@<TRIPOS>MOLECULE",17) == 0 )
        {
            //cout << "Not comment" << endl;
            //////////////////////////////
            //
            //  Process MOLECULE Record
            //
            //
            //  @<TRIPOS>MOLECULE
            //  mol_name
            //  num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
            //  mol_type
            //  charge_type
           
            if(fgets(line,MAXLINE,fp) == NULL)    
            {   
                fprintf(stderr,"No mol_name information\n");    
                fflush(stderr);
                return false;  
            }
            if(fgets(line,MAXLINE,fp) == NULL)    
            {   
                fprintf(stderr,"No num_atoms num_bonds information\n");     
                fflush(stderr);
                return false;
            }
            
            //token     = strtok(line,"\0");
            AtomNum   = atoi(strtok(line," \t")); //cout << AtomNum << endl;

            //-------- Modification (Critical) : fix bugs --------
            //ConNum    = atoi(strtok(line,"\0"));//cout << ConNum << endl;
            ConNum    = atoi(strtok(NULL," \t")); //cout << ConNum << endl;
            SBondNum  = 0;
            curAtom   = DrugATOMS;
            DrugStart = 9000;
            StrNum    = AtomNum;

            //cout << "Finish strtok" << endl;
            //////////////////////////////
            //
            //  Process Atom block
            //
            
            //  Find ATOM Record
            while(fgets(line,MAXLINE,fp) != NULL)
                if( strncmp(line,"@<TRIPOS>ATOM",13) == 0 ) break;
            
            for(i=0;i<AtomNum;i++)
            {
                if(fgets(line,MAXLINE,fp) != NULL)
                {
                    char *AtomInfo[10] = {nullptr};
                    char AtomName;
                    char AtomName2 = '*';  
                
                    j = 0;
                    token = strtok(line," \t\n");
                    while(j<10 && token != NULL )
                    {
                        AtomInfo[j++] = token;
                        token = strtok( NULL," \t\n");
                    }
                    
                    //  Dont record and do nothing when read H
                    if( (AtomName = AtomInfo[1][0]) == 'H') continue;
                    //lsr
                    ///////////////////////////////////////////////////////////////////////////////////
                    //  ���ݵL�P�_ (Ca,Mg,Na,Zn,Fe ) 
                    //  Se,Sn,Si �|�P�_�� S (���n) 
                    //  Cu,Cr.oh,Cr.th,Co.oh �|�P�_�� C (NO Charge)
                    //  Mo,Mn,K,Li,Al ���w�qtype (NO Charge) 
                    //  Du ���� C 
                    //  "."�ᤣ�Ҽ{(C.3 C.ar...) 
                    //  ANY,HEV HET HAL (QSPR ONLY) �R�� 
                    //  LP:long pair �R��
                    ///////////////////////////////////////////////////////////////////////////////////
                    if( !strcmp(AtomInfo[5],"LP")  || !strcmp(AtomInfo[5],"ANY") || 
                        !strcmp(AtomInfo[5],"HEV") || !strcmp(AtomInfo[5],"HET") ||
                        !strcmp(AtomInfo[5],"HAL") ) continue;
                    if( !strcmp(AtomInfo[5],"Du")  || !strcmp(AtomInfo[5],"Du.c")){
                        AtomName = 'C';
                    }
                    else if(AtomInfo[5][1] != '.'){
                        AtomName2 = AtomInfo[5][1];              
                    }  
                    //lsr 

                    j = DrugAtomCount++;
                        
                    curAtom[j].pos.x     = atof(AtomInfo[2]);
                    curAtom[j].pos.y     = atof(AtomInfo[3]);
                    curAtom[j].pos.z     = atof(AtomInfo[4]);
                    curAtom[j].name      = AtomName;
                    curAtom[j].SeqId     = j+1;           // info[0] 
                    curAtom[j].AtomSeqId = i+1+DrugStart; // info[0]+DrugStart
                    // for GetDrugIdx() so use original AtomSeqId with H
                    curAtom[j].AcidSeqId = 0;
                    curAtom[j].charge    = 0.0;
                    curAtom[j].numB      = 0;

                    switch(AtomName)
                    {                
                        case 'O':
                        case 'I':   curAtom[j].EngType = 0; break;
                        case 'N':   curAtom[j].EngType = 1; break;
                        case 'C':   curAtom[j].EngType = 2; break;
                        case 'S':   curAtom[j].EngType = 3; break;
                        case 'P':   curAtom[j].EngType = 4; break;
                        default :   curAtom[j].EngType = 5;
                    }
                    
                    //  Record the original position
                    DrugOrgPosi[j].x = curAtom[j].pos.x;
                    DrugOrgPosi[j].y = curAtom[j].pos.y;
                    DrugOrgPosi[j].z = curAtom[j].pos.z;
                    
                    //////////////////////////////
                    //
                    //  Create PDB format text and record it (HETATM)
                    //
                    if(AtomName2 == '*'){
                      sprintf(line,
                      //0         1         2         3    3    4    5    6
                      //123456    34     7890123  78901    9    7    5    1
                      "HETATM%5d  %-3c UNK  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                      ,curAtom[j].AtomSeqId           //  Atom serial number
                      ,curAtom[j].name                //  Atom name
                      ,curAtom[j].SeqId               //  Residue sequence number
                      ,curAtom[j].pos.x               //  coordinates for X
                      ,curAtom[j].pos.y               //  coordinates for Y
                      ,curAtom[j].pos.z               //  coordinates for Z
                      ,1.0                            //  Occupancy
                      ,0.0                            //  Temperature factor
                      );
                    } 
                    else{
                      sprintf(line,
                      //0         1         2         3    3    4    5    6
                      //123456    34     7890123  78901    9    7    5    1
                      "HETATM%5d %c%-2c  UNK  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                      ,curAtom[j].AtomSeqId           //  Atom serial number
                      ,curAtom[j].name                //  Atom name
                      ,AtomName2                      //  lsr
                      ,curAtom[j].SeqId               //  Residue sequence number
                      ,curAtom[j].pos.x               //  coordinates for X
                      ,curAtom[j].pos.y               //  coordinates for Y
                      ,curAtom[j].pos.z               //  coordinates for Z
                      ,1.0                            //  Occupancy
                      ,0.0                            //  Temperature factor
                      );
                    } 

                    strcpy(OrgDrugSTR[j],line);
                }
            }

#ifdef DEBUG
            printf("DrugAtomCount %d\n",DrugAtomCount);
            printf("Atom    Count %d\n",AtomNum);
            printf("Connect Count %d\n",ConNum);       
#endif

            //cout << "Atom process done" << endl;

            //////////////////////////////
            //
            //  Process Bond Record
            //
            //
            //  Find BOND Record
           //cout << "Find BOND Record" << endl;
            while(fgets(line,MAXLINE,fp) != NULL)
                if( strncmp(line,"@<TRIPOS>BOND",13) == 0 ) break;
            
            //  Initial Connect Struct
            for(i=0;i<DrugAtomCount;i++)
                for(j=0;j<6;j++)    
                    curAtom[i].adj[j] = -1;
            
            //  Read Bond Info
           //cout << "Read Bond Info" << endl;
           //cout << "ConNum = " << ConNum << endl;
            for(i=0;i<ConNum;i++)
            {
               //cout << "ConNum #" << i << endl;
                if(fgets(line,MAXLINE,fp) != NULL)
                {
                    char *BondInfo[5] = {nullptr};
                    int con1, con2, type;
                    
                    j = 0;
                    token = strtok(line," \t\n");
                    while(j<5 && token != NULL)
                    {
                        BondInfo[j++] = token;
                        token = strtok( NULL," \t\n");
                    }
                    
                    //cout << "After inner strtok " << endl;

                    con1 = atoi(BondInfo[1])-1;
                    con2 = atoi(BondInfo[2])-1;

                    switch(BondInfo[3][0])
                    {
                        case '1':   type = 1;   break;
                        case '2':   type = 2;   break;
                        case '3':   type = 3;   break;
                                    
                        case 'a':   //  am(amide) and ar(aromatic)
                        {
                            //if(strncmp(BondInfo[3],"am") == 0)  {}
                            //if(strncmp(BondInfo[3],"ar") == 0)  {}
                        }
                        case 'd':   //  du(dummy) 
                        case 'u':   //  un(unknown) 
                        case 'n':   //  nc(not connected)
                        default:    type = 0;   break;
                    }

                    //cout << "After BondInfo " << endl;
                   
                    if( (con1 = GetDrugIdx(con1 + DrugStart+1)) < 0 || 
                        (con2 = GetDrugIdx(con2 + DrugStart+1)) < 0 )
                        continue;
                    //cout << "After GetDrugIdx " << endl;
                    
                    if(curAtom[con1].numB>6 || curAtom[con2].numB>6)
                    {
                        fprintf(stderr,"Too much connect at Atom %d\n",
                                con1>con2 ? con1 : con2);
                        fflush(stderr);
                        return false;
                    }

                    //cout << "After curAtom > 6  " << endl;
                    //printf("%d %d %d %d\n",con1,con1idx,con2,con2idx);

                    curAtom[con1].adj[curAtom[con1].numB++] = con2;
                    curAtom[con2].adj[curAtom[con2].numB++] = con1;

                    if(type == 1)
                    {
                        //////////////////////////////
                        //
                        //  Record singlebond information
                        //
                        SBondNum = singleB[0].numSB;
                        singleB[0].atom[SBondNum][0] = con1;
                        singleB[0].atom[SBondNum][1] = con2;
                        singleB[0].type[SBondNum]    = 3;
                        singleB[0].numSB++;
                    }
                }
            }

           //cout << "Record Connection" << endl;
            //////////////////////////////
            //
            //  Record connect information using PDB format
            //
            for(i=0;i<DrugAtomCount;i++)
            {
                int  seqid, connum=0;
                char buf[MAXLINE];
                memset(buf, 0, sizeof(buf));

                // Create PDB format text and record it (CONECT)
                
                sprintf(line, "CONECT%5d", curAtom[i].AtomSeqId);
                
                for(j=0;j<6;j++)
                {
                    seqid = curAtom[curAtom[i].adj[j]].AtomSeqId;
                    
                    if(curAtom[i].adj[j] != -1)
                    {
                        sprintf(buf, "%5d", seqid);
                        strcat(line, buf);
                        connum++;
                    }
                }
                strcat(line, "\n");
#if 0        
                printf("\n %s",line); getchar();
#endif

                if(connum)
                    strcpy(OrgDrugSTR[StrNum++],line);
            }

            //////////////////////////////
            //
            //  Record single bond information using PDB format
            //
            for(i=0;i<singleB[0].numSB;i++)
            {
                sprintf(line, "SINGLE%5d%5d 2 2\n",
                        curAtom[singleB[0].atom[i][0]].AtomSeqId,
                        curAtom[singleB[0].atom[i][1]].AtomSeqId);

                strcpy(OrgDrugSTR[StrNum++],line);
            }

            OrgDrugSTRIndex = StrNum;
            mole[0].count = DrugAtomCount;  // total atoms of drug
            mole[0].atom  = DrugATOMS;
            dummy[0].atom = DummyDrugATOMS;

            //  Ok, read a drug information
            return true;
        }
        // end if @<TRIPOS>MOLECULE
    }
    //  No read
    return false;
}



