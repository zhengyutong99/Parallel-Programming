#ifndef FUILIB_Header   // �������
#define FUILIB_Header

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
using namespace std;

#ifdef WIN32
    char quote='\"';
#else 
    char quote='\'';
#endif

void ErrorExitMsg(const string&);

enum FileType   { Pdb, Mol, Mol2, List, Bad };
FileType Get_FileType(const string& file)
{
    if(      FioLib::isPDB( file ) )    return Pdb;
    else if( FioLib::isMOL( file ) )    return Mol;
    else if( FioLib::isMO2( file ) )    return Mol2;
    else if( FioLib::isLST( file ) )    return List;
    else {                                
        ErrorExitMsg("The Cavity "+file+" is not in a legal format");
    }
    return Bad; // exit only, will not return
}

struct Scoring //
{
    int     Rank;
    double  Fitness;
    double  HB_VDW;
    double  Elect;
    double  Intra;
    string  Name;
    string  FullName;
    int     Run;  
    int     Atoms;
    int     SBond;

    void print();
    string str() const;
};
bool fitness_sorter(const Scoring &a,const Scoring &b);


//-------- This class depends on fdock.h --------
class GemdockInterface 
{
public:
//-------- interfaces --------
//  call by external program/functions

    GemdockInterface(int argc,char* argv[], bool debugMode=false);
    GemdockInterface(const string& cav, const string& lig, bool debugMode=false);
    void readFromCommandLine(int argc,char* argv[]);
    void runGemdock();
    void resumeFromConfigFile(const string& conf);
    void runFromConfigFile(const string& conf);
    void printParameters();
    void promptUsage();

    // setting varibles
    void setResume(bool resume) { m_resume_mode=resume; }
    void setOutputPath(const string& path) { m_Output_path=path; }
    void setConfigFile(const string& file) { m_configure_file=file; }
    void setGenInterProfile(bool gen=false) { m_genInterProfile=gen; }
    void setPbsJobNumber(int job=0) { m_Pbs_job_number=job; }
    void setPopulation(int size=300){ m_NO_POPSIZE=size; }
    void setGeneration(int gen=70)  { m_MAXGEN=gen; }
    void setGARuntimes(int run=3)   { m_RUNTIMES=run; }
    void setFitType(char fit='6')   { m_FIT_TYPE=fit; }
    void setCharge(int charge=1)    { m_FIT_CHARGE=charge; }
    void setHydro(int hydro=1)      { m_FIT_HYDRO=hydro; }
    void setSB(char sb)             { m_FIT_SB=sb; }
    void setOperatorCau(char op='1'){ m_OPER_CAU=op; }
    void setOperatorDE(char op='0') { m_OPER_DE=op; }
    void setOperatorGau(char op='0'){ m_OPER_GAU=op; }

    // getting varibles
    string getCavityName() { return m_Cavity_file; }
    string getLigandName() { return m_Ligand_file; }
    int getPopulation() { return m_NO_POPSIZE; }
    int getGeneration() { return m_MAXGEN; }
    int getGARuntimes() { return m_RUNTIMES; }


private: 
//-------- implementations --------
//   from here, indented function means 
//   subgroup of it's root function

    enum DockType   { Single, Screen, OneOneBat, Profile };
//
// major internal processes
//
    // initialization
    void mf_initialize_variables();
    void mf_parse_command_line();

    // process parameters
    void mf_parameters_check_and_set();
        void mf_detect_set_cav_lig_path();      
        void mf_detect_docking_mode();
        void mf_detect_set_output_directory();
        void mf_parameters_cout();

    // maintain configurations
    void mf_read_configuration(const string&);
    void mf_write_configuration();
        string mf_get_configurations();

    // execute by configurations 
    void mf_config_mode(const string&);
    void mf_resume_mode(const string&);
    void mf_ranking_mode(const string&);

    // execute Gemdock
    void mf_determine_dock_mode_and_run();

    // four major program execution types
    void mf_single_docking(const string&, const string& );
    void mf_oneToOneBatch_docking(const vector<string>&, const vector<string>&);
    void mf_virtual_screening(const string&, const vector<string>&);
    void mf_generate_interaction_profile(const string&, vector<string>&);

    // details, call functions in gemdock.cpp
    void mf_protein_ligand_docking(const string& cav,const string& lig);
    void mf_protein_ligand_screening(const string& cav, const vector<string>& vec_lig);
    
//
// implementation details
//
        // converts input files into memory
        void mf_cavity_ligand_to_vectors(vector<string>&,vector<string>&);
            void mf_read_file_to_vector(const string&,vector<string>&);
            string mf_gemdock_ligand_renamer(const string& lig);

        // handle resume problems
        void mf_set_unfinished_to_list(vector<string>&, vector<string>& );
            void mf_build_success_pose_map(map<string,int>&,const string&);
            void mf_build_ideal_pose_list(const vector<string>&, 
                                          const vector<string>&, vector<string>&);
        // use external PBS binary to run Gemdock
        void mf_write_pbs_script(const string&,const string&);

        // collect and scoring results
        void mf_write_scorings_to_file(const string& dir, const string& out);
            void mf_get_best_pose(map<string,int>,const string& path, vector<Scoring>& );
                void mf_get_Scoring_from_file(Scoring&, const string& file);
            string mf_get_ranking_table(const vector<Scoring>&);

        // get/set gemdock parameters
        void mf_feed_parameters_to_gemdock();
        string mf_get_docking_command(const string&, const string&);        
        // saving best pose
        void copyBestPose(const string& path, const string& cav_prefix, const string& lig_prefix, int run);

    // debug info
    void mf_debug_io_info();

//
//-------- data members --------
//

    bool    m_debug_mode;
    bool    m_preprocess_mode;
    int     m_argc;
    char**  m_argv;
    string  m_program_name;
    DockType  m_dock_mode;   // determined by m_cavity_type and m_ligand_type
    bool    m_config_mode;
    bool    m_resume_mode;
    bool    m_genInterProfile;
    bool    m_conf_save_as;
    string  m_configure_file;
    bool    m_bootstrap;    // guard writing configure file and dispatching jobs
    string  m_Output_path;
    string  m_pbs_output_path;
    string  m_rlt_output_path;
    string  m_pose_output_path;
    string  m_best_pose_output_path;
    int     m_Pbs_job_number;   // 0 means StandardMode

    // gemdock parameters
    string  m_Cavity_file;  string m_Cavity_file_prefix;
    string  m_Ligand_file;  string m_Ligand_file_prefix;
    int     m_NO_POPSIZE;
    int     m_MAXGEN;
    char    m_FIT_TYPE;
    double  m_FIT_CHARGE;
    double  m_FIT_HYDRO;
    char    m_FIT_SB;
    int     m_ADASearchLen;
    int     m_DECSearchLen;
    int     m_UnRotSBNum;
    int     m_RUNTIMES;
    char    m_OPER_CAU;
    char    m_OPER_DE;
    char    m_OPER_GAU;

    // scoring
    int     m_CavAcidNum;       // get from poses
    int     m_ReceptorAtomNum;  // get from poses
};





//-----------------------------------------------------------------------------
//-------- interfaces --------
//-----------------------------------------------------------------------------

GemdockInterface::GemdockInterface(int argc,char* argv[], bool debugMode)
    : m_debug_mode(debugMode)
{
    mf_initialize_variables();      // initialization
    readFromCommandLine(argc,argv); // read from command line interface
}

GemdockInterface::GemdockInterface(const string& cav, 
                                   const string& lig, bool debugMode)
    : m_debug_mode(debugMode), m_Cavity_file(cav), m_Ligand_file(lig)
{
    mf_initialize_variables();      // initialization
}

void GemdockInterface::readFromCommandLine(int argc,char** argv)
{
    m_argc=argc;
    m_argv=argv;
    mf_parse_command_line();        // parse command line 
}

void GemdockInterface::runFromConfigFile(const string& conf)
{
    mf_config_mode(conf);           // directly execute program by configurations
    runGemdock();
}

void GemdockInterface::resumeFromConfigFile(const string& conf)
{
    mf_resume_mode(conf);           // resume and execute by configurations 
    runGemdock();
}

void GemdockInterface::runGemdock()
{
    mf_parameters_check_and_set();      // process parameters
    mf_write_configuration();           // maintain configurations   
    mf_determine_dock_mode_and_run();   // execute Gemdock
}

void GemdockInterface::printParameters()
{
    mf_detect_set_cav_lig_path();
    mf_detect_docking_mode();
    mf_parameters_cout();
}


void GemdockInterface::promptUsage()
{
    char p[] = "6 1 1 1 2 0 70 1 1 0 0";
    char* s = m_argv[0];
    printf("\n");
    printf("Usage   : %s -[Options] POPSIZE PDBFILE DRUGFILE   \n\n",s);
    printf("          j JOB  -- Use PBS to run                 \n"    );
    printf("          d DIR  -- Specify output directory as DIR\n"    );
    printf("          f CONF -- Use config file CONF to run    \n"    );
    printf("          n CONF -- Obtain scoring data from CONF  \n"    );
    printf("          r CONF -- Use config file CONF to resume \n"    );
    printf("          w CONF -- Save config file as CONF       \n"    );
    printf("          i      -- Generate interaction profile   \n"    );
    printf("          h      -- Show Usage                     \n\n"  );
    printf("Example : %s 100 1hvr.pdb XK2                      \n"  ,s);
    printf("Docking(PDB format):                               \n"    );
    printf("    fcdock 100 4dfr.pdb mtx.ent %s                 \n"  ,p);
    printf("Docking(MOL format):                               \n"    );
    printf("    fcdock 100 4dfr.pdb mtx.mol %s                 \n"  ,p);
    printf("Screening:                                         \n"    );
    printf("    fcdock 100 4dfr.pdb list.txt                   \n\n"  );
    exit(0);
}




//-----------------------------------------------------------------------------
// initialization
//-----------------------------------------------------------------------------

void GemdockInterface::
mf_initialize_variables()
{
    m_preprocess_mode=false;
#ifndef WIN32
    m_program_name="./gemdock";
    m_Output_path="./PrePDB/";
    // m_Output_path="/data/weichengchang/cpp/PP/Final_Project/gemdock2_src/PrePDB/";
#else
    m_program_name="gemdock";
    m_Output_path="PrePDB/";
#endif
    m_rlt_output_path=m_Output_path+"log_Rlt/";
    m_pbs_output_path=m_Output_path+"log_Pbs/";
    m_pose_output_path=m_Output_path+"docked_Pose/";
    m_best_pose_output_path=m_Output_path+"best_Pose/";
    m_Pbs_job_number=0;
    m_resume_mode =false;
    m_config_mode =false;
    m_conf_save_as=false;
    m_configure_file="gemdock.gem";
    m_bootstrap=false;
    m_genInterProfile=false;
    m_NO_POPSIZE    =200;
    m_FIT_TYPE      ='6';
    m_FIT_CHARGE    =1;
    m_FIT_HYDRO     =1;
    m_FIT_SB        ='1';
    m_ADASearchLen  =2;
    m_DECSearchLen  =2;
    m_UnRotSBNum    =0;
    m_MAXGEN        =70;
    m_RUNTIMES      =3;
    m_OPER_CAU      ='1';
    m_OPER_DE       ='0';
    m_OPER_GAU      ='0';
}

void GemdockInterface::
mf_parse_command_line()
{
    m_program_name=m_argv[0];
    int opt;
    while( (opt=getopt(m_argc,m_argv,
                       "hH?d:D:f:F:iIj:J:n:N:r:Rw:W:BXYZ")) != -1 ) 
    {
        switch(opt) {
            case 'h':
            case 'H': 
            case '?': promptUsage();
            case 'd':   
            case 'D': m_Output_path = optarg; break;            
            case 'f':
            case 'F': mf_config_mode(optarg); return;
            case 'i': 
            case 'I': m_genInterProfile=GenInterProfile=true;break;//generate interaction profile(only for post-analysis)
            case 'j':
            case 'J': m_Pbs_job_number = atoi(optarg); break;
            case 'r': mf_resume_mode(optarg); return;
            case 'n': 
            case 'N': mf_ranking_mode(optarg); return;  // from dir
            case 'R': m_resume_mode=true; break;
            case 'w':
            case 'W': m_configure_file = optarg; m_conf_save_as=true; break;

            case 'B': m_bootstrap=true;break;
            case 'X': m_preprocess_mode = true; break;
            case 'Z': m_debug_mode=true; break;
            default: break;
        }
    }

    if( optind+2 < m_argc ) {
        m_NO_POPSIZE = atoi( m_argv[optind++] );
        m_Cavity_file= m_argv[optind++];
        m_Ligand_file= m_argv[optind++];
    }
    else 
        promptUsage();   

    int base_index = optind;
    while( optind < m_argc ) {
        char* op = m_argv[ optind ];
        int current_idx = optind - base_index + 5;
        switch( current_idx )
        {
            case  5: m_FIT_TYPE   = op[0];    break; // get energy functype
            case  6: m_FIT_CHARGE = atof(op); break; // get the charge type
            case  7: m_FIT_HYDRO  = atof(op); break; // get hydrophbic type 1:(netral binding site) (<0 hydrophobic)
            case  8: m_FIT_SB     = op[0];    break; // consider the intro energy
            case  9: m_DECSearchLen = m_ADASearchLen = atoi(op); break; // get the search length 
            case 10:  {// get rotable rate of population
                     int j = atoi(op); j = ((j>=10) ? 10 : j); 
                     m_UnRotSBNum=m_NO_POPSIZE*j/10;
                     break;
            }
            case 11: // set maximum enerations
                     m_MAXGEN = atoi(op); // printf("MAXG=%d",MAXGEN);getchar();
                     if (m_MAXGEN < 20) m_MAXGEN=20;
                     break;
            case 12: // set RUNTIMES, i.e., set the number of independent runs 
                     m_RUNTIMES   = atoi(op); 
                     if (m_RUNTIMES < 1) m_RUNTIMES=1; 
                     if (m_RUNTIMES > 1000) m_RUNTIMES=1000;
                     break; 
            case 13: m_OPER_CAU   = op[0];    break; // Self-Ada Cauchy
            case 14: m_OPER_DE    = op[0];    break; // DE
            case 15: m_OPER_GAU   = op[0];    break; // Self-Gaussia
        }    

        ++optind;
    }

    if( m_debug_mode ) mf_debug_io_info();
}





//-----------------------------------------------------------------------------
// process parameters
//-----------------------------------------------------------------------------

void GemdockInterface::mf_parameters_check_and_set()
{
    mf_detect_set_cav_lig_path();       // detecting real path and set variables to real path
    mf_detect_docking_mode();           // determine DockType by FileType and set m_dock_mode
    mf_detect_set_output_directory();   // directory creation ( create directory )
    mf_parameters_cout();               // display setting
}

void GemdockInterface::mf_detect_set_cav_lig_path()
{
    if( not FioLib::isFILE(m_Cavity_file) ) {
        string old_path("../WCavPDB/"+m_Cavity_file);
        // string old_path("/data/weichengchang/cpp/PP/Final_Project/WCavPDB/"+m_Cavity_file);
        if( not FioLib::isFILE(old_path) )
            ErrorExitMsg("Cannot open "+m_Cavity_file);
        else m_Cavity_file = old_path;
    }

    m_Ligand_file = mf_gemdock_ligand_renamer( m_Ligand_file );

// set Prefix of Cavity and Ligand
    m_Cavity_file_prefix = STRING::getNoPathPrefix(m_Cavity_file);
    m_Ligand_file_prefix = STRING::getNoPathPrefix(m_Ligand_file);
}

void GemdockInterface::mf_detect_docking_mode()
{
    FileType m_cavity_type = Get_FileType( m_Cavity_file );
    FileType m_ligand_type = Get_FileType( m_Ligand_file ); 

    if( m_cavity_type==Pdb && m_ligand_type!=List )
        m_dock_mode = Single;
    else if( m_cavity_type==Pdb && m_ligand_type==List ) {
        if( m_genInterProfile ) m_dock_mode = Profile;
        else m_dock_mode = Screen;
    }
    else if( m_cavity_type==List && m_ligand_type==List )
        m_dock_mode = OneOneBat;
    else {
        if( m_cavity_type==Bad ) ErrorExitMsg("The cavity file is not legal");
        if( m_ligand_type==Bad ) ErrorExitMsg("The ligand file is not legal");
    }
    if( m_debug_mode )  mf_parameters_cout();
}

void GemdockInterface::mf_detect_set_output_directory()
{
    m_Output_path = FioLib::createAndCheckDir( m_Output_path );

    if( m_dock_mode == Profile ) return;
    
    
    m_pose_output_path  =m_Output_path+"docked_Pose/";
    m_best_pose_output_path  =m_Output_path+"best_Pose/";
    m_pose_output_path = FioLib::createAndCheckDir( m_pose_output_path );
    m_best_pose_output_path = FioLib::createAndCheckDir( m_best_pose_output_path );

    m_rlt_output_path   =m_Output_path+"log_Rlt/";
    if( m_dock_mode==Single or m_dock_mode==OneOneBat )
        m_rlt_output_path = FioLib::createAndCheckDir( m_rlt_output_path );

    m_pbs_output_path   =m_Output_path+"log_Pbs/";
    if( m_Pbs_job_number > 0 )
        m_pbs_output_path = FioLib::createAndCheckDir( m_pbs_output_path );
    
    if( not m_config_mode and not m_conf_save_as )
        m_configure_file = m_Output_path+"gemdock_out/"+m_configure_file;
}

void GemdockInterface::mf_parameters_cout()
{
    cout << mf_get_configurations();
    cout << "Docking Mode    = ";
    switch(m_dock_mode) {
        case 0: cout << "Single" << endl;  break;
        case 1: cout << "Screen" << endl;  break;
        case 2: cout << "OneOneBat" << endl; break;
        case 3: cout << "Profile" << endl; break;
    }
    if( m_dock_mode==Single )
    cout << "RLT log dir     = " << m_rlt_output_path << endl;
    if( m_Pbs_job_number > 0 )
    cout << "PBS log dir     = " << m_pbs_output_path << endl;
    //cout << "Configuration   = " << m_configure_file << endl; 
    cout << "Docked Poses    = " << m_pbs_output_path << endl;
    cout << endl;
    cout.flush();
}





//-----------------------------------------------------------------------------
// maintain configurations
//-----------------------------------------------------------------------------

void GemdockInterface::mf_read_configuration(const string& file)
{
    ifstream fin( file.c_str() );
    if( not fin.good() or FioLib::isDIR(file) ) 
        ErrorExitMsg("Error: Cannot open configure file - "+file);

    string s; 
    string key, op, value;
    for(int i=1; getline(fin,s); i++ )   
    {        
        if( s.size()==0 or s[0]=='#' or s[0]=='[')  continue;   // ignore comment or empty line
        
        string::size_type pos = s.find("="); 
        if( pos == string::npos ) continue;   // ignore lines with '=' mark        
        op = s[pos];
        
        key = s.substr( 0, pos ); istringstream skey(key); skey >> key;    
        value = s.substr( pos+1 ); istringstream svalue(value); //svalue >> value;

        if( key.size()==0 or op.size()==0 or value.size()==0 or op==value ) 
            continue;                                // ignore lost key, no operator, or no value
        // begin reading
             if( key == "PBSJobNumber" ) m_Pbs_job_number = atoi(value.c_str()) ;
        else if( key == "Cavity"       ) m_Cavity_file    = value ;
        else if( key == "Ligand"       ) m_Ligand_file    = value ;
        else if( key == "OutputPath"   ) m_Output_path    = value ;
        else if( key == "NO_POPSIZE"   ) m_NO_POPSIZE     = atoi(value.c_str()) ;
        else if( key == "FIT_TYPE"     ) m_FIT_TYPE       = value[0] ;
        else if( key == "FIT_CHARGE"   ) m_FIT_CHARGE     = atof(value.c_str()) ;
        else if( key == "FIT_HYDRO"    ) m_FIT_HYDRO      = atof(value.c_str()) ;
        else if( key == "FIT_SB"       ) m_FIT_SB         = value[0] ;
        else if( key == "ADASearchLen" ) m_ADASearchLen   = atoi(value.c_str()) ;
        else if( key == "DECSearchLen" ) m_DECSearchLen   = atoi(value.c_str()) ;
        else if( key == "UnRotSBNum"   ) m_UnRotSBNum     = atoi(value.c_str()) ;
        else if( key == "MAXGEN"       ) m_MAXGEN         = atoi(value.c_str()) ;
        else if( key == "RUNTIMES"     ) m_RUNTIMES       = atoi(value.c_str()) ;
        else if( key == "OPER_CAU"     ) m_OPER_CAU       = value[0] ;
        else if( key == "OPER_DE"      ) m_OPER_DE        = value[0] ;
        else if( key == "OPER_GAU"     ) m_OPER_GAU       = value[0] ;
        else {
            //ostringstream oss;
            //oss << "Warning! In file: "+ file +"  Line: " << i << " -- Not a valid option: " << key << " =  " << value  << endl;
            //UniversalErrorMessageLogger(oss.str());
        }
    }
    fin.close();
    fin.clear();

    if(m_debug_mode) mf_debug_io_info();
}

void GemdockInterface::mf_write_configuration()
{
    if( m_bootstrap||m_genInterProfile ) return;
    ofstream fout(m_configure_file.c_str());
    if( not fout.good() )
        ErrorExitMsg("Error: Cannot write configure file - " + m_configure_file);
    fout << mf_get_configurations();
    if(m_debug_mode)  { cout << mf_get_configurations(); getchar(); }
    fout.close();
}

string GemdockInterface::mf_get_configurations()
{
    ostringstream oss; 
    oss << "#########################################"<< endl
        << "# GEMDOCK CONFIGURATION SETTINGS"         << endl
        << "# Format:   key = value"                  << endl
        << "# Comment:  leading with #"               << endl
        << "#########################################"<< endl
        << "                                         "<< endl
        <<  "Cavity          = " << m_Cavity_file     << endl
        <<  "Ligand          = " << m_Ligand_file     << endl
        <<  "OutputPath      = " << m_Output_path     << endl
        <<  "NO_POPSIZE      = " << m_NO_POPSIZE      << endl
        <<  "FIT_TYPE        = " << m_FIT_TYPE        << endl 
        <<  "FIT_CHARGE      = " << m_FIT_CHARGE      << endl 
        <<  "FIT_HYDRO       = " << m_FIT_HYDRO       << endl 
        <<  "FIT_SB          = " << m_FIT_SB          << endl 
        <<  "ADASearchLen    = " << m_ADASearchLen    << endl 
        <<  "DECSearchLen    = " << m_DECSearchLen    << endl
        <<  "UnRotSBNum      = " << m_UnRotSBNum      << endl 
        <<  "MAXGEN          = " << m_MAXGEN          << endl
        <<  "RUNTIMES        = " << m_RUNTIMES        << endl
        <<  "OPER_CAU        = " << m_OPER_CAU        << endl
        <<  "OPER_DE         = " << m_OPER_DE         << endl
        <<  "OPER_GAU        = " << m_OPER_GAU        << endl
        <<  "PBSJobNumber    = " << m_Pbs_job_number  << endl
    ;    
    return oss.str();
}





//-----------------------------------------------------------------------------
// execute by configurations 
//-----------------------------------------------------------------------------

void GemdockInterface::mf_config_mode(const string& file)
{
    mf_read_configuration(file);
    m_config_mode=true;
}

void GemdockInterface::mf_resume_mode(const string& file)
{
    mf_read_configuration(file);
    m_resume_mode=true;
}

void GemdockInterface::mf_ranking_mode(const string& file)
{
    mf_read_configuration(file);   
    mf_parameters_check_and_set();
    mf_write_scorings_to_file(m_pose_output_path,m_Output_path);
    exit(0);
}





//-----------------------------------------------------------------------------
// execute Gemdock
//-----------------------------------------------------------------------------

void GemdockInterface::mf_determine_dock_mode_and_run()
{
    vector<string> vec_cav,vec_lig; // save list into memory
    mf_cavity_ligand_to_vectors(vec_cav,vec_lig);

    if( m_resume_mode && !m_bootstrap ) {   // only for Screen and OneOneBat
        mf_set_unfinished_to_list(vec_cav,vec_lig);
        if( vec_lig.size() <= 0 )
            mf_write_scorings_to_file(m_pose_output_path,m_Output_path);
        else
            cout << endl <<"Start in resume mode" << endl;
    }

    mf_feed_parameters_to_gemdock();// convert internal info(except cav,lig) to GEMDOCK global variables
    switch( m_dock_mode ) 
    {
        case Single:    mf_single_docking(vec_cav[0],vec_lig[0]); return;
        case OneOneBat: mf_oneToOneBatch_docking(vec_cav,vec_lig); return;
        case Screen:    mf_virtual_screening(vec_cav[0],vec_lig); return;
        case Profile:   mf_generate_interaction_profile(vec_cav[0],vec_lig); return;
    }
}





//-----------------------------------------------------------------------------
// four major program execution types
//-----------------------------------------------------------------------------

void GemdockInterface::mf_single_docking(const string& cav, const string& lig)
{
    if( m_Pbs_job_number <= 0 || ! FioLib::isPBSexist() ) {   // non-pbs mode
        if( m_preprocess_mode ) 
            return;
        mf_protein_ligand_docking( cav,lig ); 
    }
    else {  // pbs mode
        string lig_prefix( STRING::getNoPathPrefix( lig ) );
        string pbs_script_location( m_pbs_output_path 
                                    + lig_prefix +"-pbs.sh" );
        string dock_cmd( mf_get_docking_command(cav,lig) );
        mf_write_pbs_script( dock_cmd, pbs_script_location );
    }
}

void GemdockInterface::mf_oneToOneBatch_docking(const vector<string>& vec_cav, 
                         const vector<string>& vec_lig)
{
    if( vec_cav.size() != vec_lig.size() ) // non-pbs mode
        ErrorExitMsg("Error: The input list must in pair when running batch docking mode");
    
    if( m_Pbs_job_number <= 0 || ! FioLib::isPBSexist() ) {   // non-pbs mode
        if( m_preprocess_mode ) return;
        for(size_t i=0;i<vec_cav.size();++i) 
            mf_protein_ligand_docking(vec_cav[i],vec_lig[i]);
    }
    else {  // pbs mode
        m_bootstrap=true;

        int totalNumber = (int)vec_cav.size();
        int jobInEachBatch = abs(totalNumber/m_Pbs_job_number)+1;

        int currentBathId = 1;
        ostringstream oss_command;
        for(int i=0;i<totalNumber; ++i) {
            oss_command << mf_get_docking_command(vec_cav[i],vec_lig[i]) << endl;
            if( (i+1)%jobInEachBatch==0 or (i+1)==totalNumber ) 
            {
                ostringstream id; id << setw(2) << setfill('0') << currentBathId;
                string pbs_script(m_pbs_output_path+
                                  m_Ligand_file_prefix+"-pbs"+id.str()+".sh");
                
                string pbs_cmd( oss_command.str() );
                mf_write_pbs_script(pbs_cmd,pbs_script);
                oss_command.str("");
                currentBathId++;
            }
        }
        m_bootstrap=false;
    }
}

void GemdockInterface::mf_virtual_screening(const string& cav, const vector<string>& vec_lig)
{
    if( m_Pbs_job_number <= 0 || ! FioLib::isPBSexist() ) {   // non-pbs mode
        if( m_preprocess_mode ) return;
        mf_protein_ligand_screening(cav,vec_lig);   
    }
    else {
        m_bootstrap=true;
        int totalNumber = (int)vec_lig.size();
        int jobInEachBatch = abs(totalNumber/m_Pbs_job_number)+1;

        int currentBathId = 1;
        ostringstream oss_ligand_file;
        for(int i=0;i<totalNumber; ++i) {
            oss_ligand_file << quote <<  vec_lig[i] << quote << endl;
            if( (i+1)%jobInEachBatch==0 or (i+1)==totalNumber ) 
            {
                ostringstream id; id << setw(2) << setfill('0') << currentBathId;
                string pbs_script(m_pbs_output_path+
                                  m_Ligand_file_prefix+"-pbs"+id.str()+".sh");
                string ligand_list(m_pbs_output_path+
                                   m_Ligand_file_prefix+"-pbs"+id.str()+".lst");
                
                ofstream splitted_ligand_list( ligand_list.c_str() );
                splitted_ligand_list << oss_ligand_file.str();

                string pbs_cmd( mf_get_docking_command(cav,ligand_list) );
                mf_write_pbs_script(pbs_cmd,pbs_script);
                oss_ligand_file.str("");
                currentBathId++;
            }
        }
        m_bootstrap=false;
    }
}

void GemdockInterface::mf_generate_interaction_profile(const string& cav, 
                                                       vector<string>& lig)
{
// set up cavity
    ReadPDBCavity(cav);

// write to old log file
    string result_log(m_Output_path+"gemdock_out/"+m_Ligand_file_prefix+"_dock.log");
    if( (dlogfp = fopen(result_log.c_str(),"w")) == NULL )
        ErrorExitMsg("Error: Docking result log file "+result_log+" open error!!");
    char buf[512];
    fprintf(dlogfp,"#%s",m_Cavity_file_prefix.c_str());
    for(int i=0 ; i<mole[1].count ;i++) {
        sprintf(buf, "%d%s-%d%c",mole[1].atom[i].AcidSeqId,(char*)( (mole[1].atom[i].AcidType == -1) ? mole[1].atom[i].OtherName : acidName[(int)mole[1].atom[i].AcidType]),mole[1].atom[i].AtomSeqId,mole[1].atom[i].name);
        fprintf(dlogfp,"\t%s-V\t%s-H\t%s-E", buf,buf,buf);
    } 
    fprintf(dlogfp,"\tVDW\tHBond\tElec\tAverConPair\tTotalEnergy\n");

    // reverse header
    fprintf(dlogfp,"#%s",m_Cavity_file_prefix.c_str());
    for(int i=0 ; i<mole[1].count ;i++) {        /* protein */  
        sprintf(buf, "%d%s-%d%c",mole[1].atom[i].AcidSeqId,(char*)( (mole[1].atom[i].AcidType == -1) ? mole[1].atom[i].OtherName : acidName[(int)mole[1].atom[i].AcidType]),mole[1].atom[i].AtomSeqId,mole[1].atom[i].name);
        fprintf(dlogfp,"\tV-%s\tH-%s\tE-%s", buf,buf,buf);
    }
    fprintf(dlogfp,"\tVDW\tHBond\tElec\tAverConPair\tTotalEnergy\n");

// write to new log file
    oss_fitness.str("");    // clean oss_fitness
    oss_fitness << "#"<< m_Cavity_file << " with " << m_Ligand_file << endl;
    oss_fitness <<"#Ligand\tTotalEnergy\tVDW\tHBond\tElec\tAverConPair\t" << endl;//lsr

    oss_profile.str("");    // clean oss_profile
    oss_profile << "#"<< m_Cavity_file << " with " << m_Ligand_file << endl;
    oss_profile << "#" << m_Cavity_file_prefix;
    for(int i=0 ; i<mole[1].count ;i++) {
        sprintf(buf, "%d%s-%d%c",mole[1].atom[i].AcidSeqId,(char*)( (mole[1].atom[i].AcidType == -1) ? mole[1].atom[i].OtherName : acidName[(int)mole[1].atom[i].AcidType]),mole[1].atom[i].AtomSeqId,mole[1].atom[i].name);
        oss_profile << tab << buf << "-V" << tab << buf << "-H" << tab << buf << "-E";
    }
    oss_profile << endl;
    // reverse header
    oss_profile << "#" << m_Cavity_file_prefix;
    for(int i=0 ; i<mole[1].count ;i++) {
        sprintf(buf, "%d%s-%d%c",mole[1].atom[i].AcidSeqId,(char*)( (mole[1].atom[i].AcidType == -1) ? mole[1].atom[i].OtherName : acidName[(int)mole[1].atom[i].AcidType]),mole[1].atom[i].AtomSeqId,mole[1].atom[i].name);
        oss_profile << tab << "V-" << buf  << tab << "H-" << buf  << tab << "E-" << buf ;
    }
    oss_profile << endl;

// begin generating profile
    string ligand_prefix;
    for(size_t idx=0; idx<lig.size();++ idx)
    {
        if( not ReadDrugCompound(lig[idx]) ) 
            continue;

        ligand_prefix=STRING::RemovePath(lig[idx]);
        oss_fitness << ligand_prefix;
        oss_profile << ligand_prefix;
        fprintf(dlogfp,"%s",ligand_prefix.c_str());
        FunctionEvaluation(&nowInd, -1);
    }
    fclose(dlogfp);

    ofstream fout_fitness, fout_profile;
    string   fitness( m_Output_path+"gemdock_out/"+"fitness.txt");
    string   profile( m_Output_path+"gemdock_out/"+"profile.txt");
    fout_fitness.open(fitness.c_str()); 
    fout_profile.open(profile.c_str());
    if( fout_profile.good() && fout_fitness.good() ) {
        cout << "write fitness file to " << fitness << endl;
        cout << "write interaction file to " << profile << endl; cout.flush();
    }
    fout_fitness << oss_fitness.str();
    fout_profile << oss_profile.str();
}





//-----------------------------------------------------------------------------
// details, call functions in gemdock.cpp
//-----------------------------------------------------------------------------

void GemdockInterface::mf_protein_ligand_docking(const string& cav, const string& lig)
{
// read cavity
    ReadPDBCavity(cav);

// read ligand
    if( not ReadDrugCompound(lig) ) ErrorExitMsg("Error: Read Drug Compound Error");

    LigBaseEnergy = GetLigBaseEnergy();

// open runtime log
    string chromosome_log(m_rlt_output_path
                  +"T-"+STRING::getNoPathPrefix(cav)
                  +"-"+STRING::getNoPathPrefix(lig)
                  +"-"+STRING::Num2Str(m_NO_POPSIZE)
                  +".txt");
    if ( (PerFP=fopen(chromosome_log.c_str(),"w")) == NULL ) 
        ErrorExitMsg("Warning: RLT File "+chromosome_log+"  open error!");

// prompt info
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

    //AdjustReceptorPosi();       // adjust recpertor axis
    AdjustDrugPosi();
    InitRandom();

// start GA
    for(int i=0;i<RUNTIMES;i++) {
        Itoa(i, RunId, 10);
        trimAllSpace(RunId);
        InitGlobalVar();
        FCES();
    } 
}

void GemdockInterface::mf_protein_ligand_screening(const string& cav, 
                                                   const vector<string>& lig)
{
    // open runtime log and result log
    string runtime_log(m_Output_path+"gemdock_out/"+m_Ligand_file_prefix+"_run.log");
    if( (rlogfp = fopen(runtime_log.c_str(),"w")) == NULL ) 
        ErrorExitMsg("Error: Runtime log file " + runtime_log + " open error!!");
    string result_log(m_Output_path+"gemdock_out/"+m_Ligand_file_prefix+"_dock_log");
    if( (dlogfp = fopen(result_log.c_str(),"w")) == NULL )
        ErrorExitMsg("Error: Docking result log file "+result_log+" open error!!");

// set up cavity
    ReadPDBCavity(cav);

// begin screening
    double bestfit = 9e20;  // for recording best run for each drug
    int    bestidx = 0;     // for recording best run for each drug
    string ligand_prefix;   // ligand prefix

    DResult.clear();
    for(size_t idx=0; idx<lig.size();++idx)
    {
        ligand_prefix=STRING::getNoPathPrefix(lig[idx]);
        if( not ReadDrugCompound(lig[idx]) ) 
            continue;
        

        printf("\r\n\r\n");
        fflush(stdout);
        printf("Protein  : %-20.20s Atom Num : %4d  Receptor Atom Num : %4d \n"
               "%s   : %-20.20s Atom Num : %4d  Single Bond Num   : %4d \n",
               PDBName, Protein.acidNum, ReceptorAtomCount,
               (DrugFormat > 1) ? "Drug  " : "Ligand",
               LigandName, DrugAtomCount, singleB[0].numSB);

        // ajust the centers of drug & Protein to original point
        AdjustDrugPosi();
        InitRandom();

        bestfit = 9e20;
        for(int i=0;i<RUNTIMES;i++) {
            Itoa(i, RunId, 10);
            trimAllSpace(RunId);
            InitGlobalVar();
            FCES();
            //  Record and Write best result into runtime log
            if(bestfit > BEST_IND.value) {
                bestfit = BEST_IND.value;
                bestidx = i;
            }
            fprintf(rlogfp,"%2d %14f %s\n",i, BEST_IND.value, LigandName);
        } // for
        cout << "Docking complete (ligand: " << lig[idx] << ")" << endl;
        cout << "Docking fitness = " << BEST_IND.value << endl;

        setDResult( m_Cavity_file_prefix, ligand_prefix,
                    DrugAtomCount, BestHBondCount, BestElectCount, 
                    bestfit, bestidx );     
        copyBestPose(m_Output_path,m_Cavity_file_prefix, ligand_prefix, bestidx);
    }
    dr_qsort();

    writeDResult(dlogfp);

    fclose(rlogfp); 
    fclose(dlogfp);
}




//-----------------------------------------------------------------------------
// converts input files into memory
//-----------------------------------------------------------------------------

void GemdockInterface::mf_cavity_ligand_to_vectors(vector<string>& vec_cav,
                                                   vector<string>& vec_lig)
{
    for(size_t i=0;i<vec_lig.size();++i) 
        vec_lig[i]=mf_gemdock_ligand_renamer(vec_lig[i]);

    switch( m_dock_mode ) 
    {
        case Single:    vec_cav.push_back(m_Cavity_file);
                        vec_lig.push_back(m_Ligand_file);
                        break;
        case Screen:    
        case Profile:   
                        vec_cav.push_back(m_Cavity_file);
                        mf_read_file_to_vector(m_Ligand_file,vec_lig);
                        break;
        case OneOneBat: 
                        mf_read_file_to_vector(m_Cavity_file, vec_cav);
                        mf_read_file_to_vector(m_Ligand_file,vec_lig);
                        break;

    }

    if( vec_lig.size() <= 0 )
        ErrorExitMsg("Error: There are no ligand for docking");
    if( vec_cav.size() <= 0 )
        ErrorExitMsg("Error: There are no protein for docking");
}

void GemdockInterface::mf_read_file_to_vector(const string& file, vector<string>& vec)
{
    ifstream fin(file.c_str());
    if( not fin.good() ) ErrorExitMsg("Error: Cannot open file "+ file);
    
    vec.clear();
    string read_line;
    while(getline(fin,read_line))
        if( read_line.size()> 0 && read_line[0]!='#')
            vec.push_back(read_line);               
    fin.close();
}

string GemdockInterface::mf_gemdock_ligand_renamer(const string& lig)
{
    vector<string> search_path;
    search_path.push_back( "../Drug/" + lig );
    search_path.push_back( "../Ligand/" + lig );
    // search_path.push_back( "/data/weichengchang/cpp/PP/Final_Project/Ligand/" + lig );

    if( not FioLib::isFILE( lig ) ) 
    {
        for(size_t i=0;i<search_path.size();++i)
        {
            if( FioLib::isFILE( search_path[i] ) ) {
                //cout << "Correct path = " << search_path[i]; getchar();
                return search_path[i];
            }
        }
        ErrorExitMsg("Cannot open ligand file"+lig);
    }
    return lig;
}





//-----------------------------------------------------------------------------
// handle resume problems
//-----------------------------------------------------------------------------

void GemdockInterface::mf_set_unfinished_to_list(vector<string>& vec_cav,vector<string>& vec_lig)
{
    vector<string> overall_output_list;     // �Ӳ��ͪ�poses
    mf_build_ideal_pose_list(vec_cav,vec_lig,overall_output_list);

    map<string,int> current_pose_status;    // �έp�ثe�������Tposes
    mf_build_success_pose_map( current_pose_status, m_pose_output_path );
    
    vector<string> left_cav, left_lig;
    int boundary = getGARuntimes();         // �n���D�h��run�~�����ǬO������poses
    for(size_t i=0;i<overall_output_list.size();++i) {
        int current_status = current_pose_status[ overall_output_list[i] ];
        if( current_status != boundary ) 
        {
            if( m_dock_mode==OneOneBat )
                left_cav.push_back( vec_cav[i] );
            left_lig.push_back( vec_lig[i] );
        }
    }
    if( m_dock_mode==OneOneBat ) vec_cav = left_cav;
    vec_lig = left_lig;

    if(m_debug_mode) {
        cout << "unfinished cavity" << endl; size_t i; for(i=0;i<vec_cav.size();++i)   cout << vec_cav[i] << endl;
        cout << "unfinished ligand" << endl; for(i=0;i<vec_lig.size();++i)   cout << vec_lig[i] << endl;
    }
}

void GemdockInterface::mf_build_ideal_pose_list( const vector<string>& cav, 
                          const vector<string>& lig, vector<string>& list)
{
    if( cav.size()<=0 || lig.size()<=0 ) return;

    list.clear();

    string ideal_name;
    string cavity;
    string ligand;
    for(size_t i=0;i<lig.size();++i)
    {
        if(m_dock_mode==Screen) cavity=cav[0];  // less efficiency, but clear
        else cavity=cav[i];

        ideal_name = //m_pose_output_path + 
                     STRING::getNoPathPrefix( cavity ) + '-' +
                     STRING::getNoPathPrefix( lig[i] );
        list.push_back( ideal_name );
        //cout << "ideal_name = " << ideal_name << endl;
    }
}

void GemdockInterface::mf_build_success_pose_map(map<string,int>& status, const string& path)
{
// need to scan output directory
    //cout << "Parse poses from " << path << endl;
    DIR* dir=FioLib::OpenDIR(path);
    struct dirent* dp;
    while( (dp=readdir(dir))!=NULL )
    if( dp->d_name[0] != '.' ) 
    {
        string item( STRING::getNoPathPrefix( dp->d_name ) );
        string::size_type pos = item.find_last_of('-');
        if( pos != string::npos )
            item = item.substr(0,pos);
        ++status[ item ];
    }
    closedir(dir);
}






//-----------------------------------------------------------------------------
// use external PBS binary to run Gemdock
//-----------------------------------------------------------------------------

void GemdockInterface::mf_write_pbs_script(const string& instruction,
                                           const string& script_file)
{
    ostringstream oss;
    oss << "#!/bin/sh" << endl;
    oss << "### Generated by GEMDOCK ###" << endl;
    oss << "#PBS -l nodes=1,walltime=9900:00:00" << endl;
    oss << "#PBS -N " << STRING::getNoPathPrefix(script_file) << endl;           // job name
    oss << "#PBS -o " << m_pbs_output_path
                      << STRING::getNoPathPrefix(script_file) 
                      << ".log" << endl; // output log name
    oss << "#PBS -j oe" << endl;
    oss << "cd $PBS_O_WORKDIR"  << endl;
    oss << instruction << endl;
    
    ofstream fout( script_file.c_str() );
    fout << oss.str();
    fout.close();

    if( m_preprocess_mode ) return;

    string command = "qsub " + script_file ;
    
    if (!m_debug_mode) {
        cout << "submiting job : "; cout.flush();
    }

    system( command.c_str() );
}






//-----------------------------------------------------------------------------
// get/set gemdock parameters
//-----------------------------------------------------------------------------

string GemdockInterface::mf_get_docking_command(const string& cavity, 
                                                const string& ligand)
{
    string enda(" ");
    ostringstream oss;
    oss << m_program_name << enda;

    if( m_genInterProfile ) 
        oss << "-i" << enda;
    if( m_bootstrap ) 
        oss << "-B" << enda;
    oss << "-d " << quote << m_Output_path << quote << enda;

    oss << m_NO_POPSIZE << enda;
    oss << quote << cavity << quote << enda;
    oss << quote << ligand << quote << enda;
    oss << m_FIT_TYPE << enda;    
    oss << m_FIT_CHARGE << enda;    
    oss << m_FIT_HYDRO << enda;    
    oss << m_FIT_SB << enda;    
    oss << m_ADASearchLen << enda;
    oss << (10*m_UnRotSBNum/m_NO_POPSIZE) << enda;
    oss << m_MAXGEN << enda;
    oss << m_RUNTIMES << enda;
    oss << m_OPER_CAU << enda;
    oss << m_OPER_DE << enda;
    oss << m_OPER_GAU << enda;

    return oss.str();
}

void GemdockInterface::mf_feed_parameters_to_gemdock()
{
    //GemdockOutputPath=m_Output_path;
    PreStrPath2         =   m_pose_output_path;
    GenInterProfile     =   m_genInterProfile;
    NO_POPSIZE          =   m_NO_POPSIZE;
    MAXGEN              =   m_MAXGEN;
    FIT_TYPE            =   m_FIT_TYPE;
    FIT_CHARGE          =   m_FIT_CHARGE;
    FIT_HYDRO           =   m_FIT_HYDRO;
    FIT_SB              =   m_FIT_SB;
    ADASearchLen        =   m_ADASearchLen;
    DECSearchLen        =   m_DECSearchLen;
    UnRotSBNum          =   m_UnRotSBNum;
    RUNTIMES            =   m_RUNTIMES;
    OPER_CAU            =   m_OPER_CAU;
    OPER_DE             =   m_OPER_DE;
    OPER_GAU            =   m_OPER_GAU;    
    if( m_dock_mode==Screen ) DOCKING=0;
    if( m_bootstrap) isInPbsMode=true; 
}








//-----------------------------------------------------------------------------
// collect and scoring results
//-----------------------------------------------------------------------------

void GemdockInterface::mf_write_scorings_to_file(const string& path, const string& output)
{
    string path_name( FioLib::checkDirRear(path) );

    map<string,int> current_pose_status;
    mf_build_success_pose_map( current_pose_status, path_name );


    if( m_dock_mode==Screen )
    {
        vector<Scoring> best_poses;
        mf_get_best_pose( current_pose_status, path_name, best_poses );

        sort( best_poses.begin(), best_poses.end(), fitness_sorter);

        if( best_poses.size() <= 0 ) {
            ostringstream oss;
            oss << "Warning: There's no poses generated in the given settings" << endl; 
            UniversalErrorMessageLogger(oss.str());
        }

        string ranking_file( output + "ranking.txt" );

        ofstream fout( ranking_file.c_str() );
        fout << mf_get_ranking_table( best_poses );
        //cout << mf_get_ranking_table( best_poses );
        fout.close();
        ostringstream oss;
        oss <<"Message: The scoring data has been dumped to " << ranking_file << endl;
        UniversalErrorMessageLogger(oss.str());
    }

    if( m_dock_mode==OneOneBat )
    {
        ostringstream rlt_list;
        DIR* dir=FioLib::OpenDIR(m_rlt_output_path);
        struct dirent* dp;
        while( (dp=readdir(dir))!=NULL )
        if( dp->d_name[0] == 'T' ) 
        {
            rlt_list << m_rlt_output_path << dp->d_name << endl;
        }
        closedir(dir);
        
        string analysis_file(m_Output_path+"rlt_list.txt");
        string analy_rlt_file(m_Output_path+"rlt_analysis.txt");
        ofstream fout( analysis_file.c_str() );
        fout << rlt_list.str() ;
        fout.close();
        getRltAnalysis( analysis_file.c_str(), analy_rlt_file.c_str() );

        ostringstream oss;
        UniversalErrorMessageLogger(oss.str());
        oss <<"Message: The analysis data has been dumped to " << analy_rlt_file << endl;
    }
}

string GemdockInterface::mf_get_ranking_table(const vector<Scoring>& table)
{
    ostringstream info;
    info  << "# Top best docking result" << endl
          << "#" << endl
          << "# Cavity: " << m_Cavity_file_prefix
          << "\tResNum: " << m_CavAcidNum
          << "\tAtomNum: "<< m_ReceptorAtomNum << endl
          << "# PoseLocation: " << m_pose_output_path << endl
          << "# RANK\tDrugName\tFitnessValue\tHB&VDW\tElect\tIntraEnergy\tAtoms\tSingle\tBestRun\tPoseName" << endl
          ;

    for(size_t i=0;i<table.size();i++) {
        info << i+1 << "\t" << table[i].str();  
        //info << m_pose_output_path << table[i].FullName << endl;
    }
    return info.str();
}

void GemdockInterface::mf_get_best_pose(map<string,int> all_poses,const string& path,
    vector<Scoring>& best_pose)
{
    map<string,int>::iterator current;
    string path_name( FioLib::checkDirRear(path) );
    best_pose.clear();

    int run = getGARuntimes();
    for(current=all_poses.begin(); current!=all_poses.end();++current)
    if( current->second == run )
    {
        double minFit=9999;
        Scoring minPose;
        for(int i=0;i<run; ++i) 
        {
            string pose_name(path_name+current->first
                             +'-'+STRING::itos(i)+".pdb");
            ///cout << pose_name << endl;
            Scoring sb;
            mf_get_Scoring_from_file(sb,pose_name);
            if(sb.Fitness < minFit )
            {
                minFit = sb.Fitness;
                minPose = sb;   // bad efficiency, but read easily
            }
        }
        best_pose.push_back(minPose);
        //cout << "==> " << minPose.Fitness << " " << minPose.Name << endl;
    }
}

void GemdockInterface::mf_get_Scoring_from_file(Scoring& sb, const string& file)
{

    ifstream fin( file.c_str() );
    if( ! fin.good() ) {
        ostringstream oss;
        oss << "Warning: Cannot open pose file - " << file << endl;
        UniversalErrorMessageLogger(oss.str());
        return;   
    }
    bool readCavHeader = false;
    string line;
    sb.FullName = STRING::RemovePath(file);
    while(getline(fin,line)) 
    if( line.size() > 50 ) if( line.substr(0,6) == "REMARK" )
    {
        
        if( not readCavHeader && line.substr(0,17) == "REMARK Protein : ")
        {
            istringstream sin;
            sin.str(line.substr(50-1));   sin >> m_CavAcidNum ;
            sin.str(line.substr(75-1));   sin >> m_ReceptorAtomNum;
            readCavHeader = true;
        }
        else if( line.substr(0,17) == "REMARK Drug    : ")
        {
            istringstream sin;
            sin.str(line.substr(17-1));   sin >> sb.Name ;
            sin.str(line.substr(50-1));   sin >> sb.Atoms ;
            sin.str(line.substr(75-1));   sin >> sb.SBond ;
        }
        else if( line.substr(0,19) == "REMARK FitnessValue" )
        {
            getline(fin,line);
            istringstream sin;
            sin.str(line.substr(8-1));   sin >> sb.Fitness;
            sin.str(line.substr(23-1));  sin >> sb.HB_VDW;
            sin.str(line.substr(38-1));  sin >> sb.Elect;
            sin.str(line.substr(53-1));  sin >> sb.Intra;

            sb.Run = getGARuntimes();
            break;
        }
    }
}

void GemdockInterface::copyBestPose(const string& path, const string& cav_prefix, const string& lig_prefix, int run)
{
    //sprintf(buf,"%s%s-%s-%s.pdb",PreStrPath2.c_str(),PDBName,LigandName,RunId);
    string best_pose_name(m_pose_output_path+cav_prefix+'-'+lig_prefix+='-'+STRING::itos(run)+".pdb");
    string copy_pose_name(m_best_pose_output_path+cav_prefix+'-'+lig_prefix+='-'+STRING::itos(run)+".pdb");
    cout << endl;    cout.flush();
    cout << "Best Pose: " << copy_pose_name << "ZZZKKKZZZ" << endl;
    cout << endl << endl; cout.flush(); fflush(stdout);
    FioLib::copyFile(best_pose_name,copy_pose_name);
}





//-----------------------------------------------------------------------------
// debug info
//-----------------------------------------------------------------------------

void GemdockInterface::mf_debug_io_info()
{
    printf("NO_POPSIZE=%d\n",m_NO_POPSIZE);
    printf("FIT_TYPE=%c\nFIT_CHARGE=%f\nFIT_HYDRO=%f\nFIT_SB=%c\n",
            m_FIT_TYPE,m_FIT_CHARGE,m_FIT_HYDRO,m_FIT_SB);
    printf("DECSearchLen=%d =ADASearchLen=%d\n\n", 
            m_DECSearchLen, m_ADASearchLen);
    printf("UnRotSBNum=%d\n",m_UnRotSBNum);
    printf("MAXGEN=%d\nRUNTIMES=%d\n",m_MAXGEN,m_RUNTIMES);
    printf("OPER_CAU=%c OPER_DE=%c OPER_GAU=%c\n", 
            m_OPER_CAU, m_OPER_DE, m_OPER_GAU);
}

void ErrorExitMsg(const string& s)
{
    ostringstream oss;
    oss << s << endl;
    UniversalErrorMessageLogger(oss.str());
    exit(0);
}





//-----------------------------------------------------------------------------
// Scoring
//-----------------------------------------------------------------------------

void Scoring::print() { 
    cout << " Name "    << Name 
         << " Fitness " << Fitness                    
         << " HB_VDW "  << HB_VDW                      
         << " Elect "   << Elect                        
         << " Intra "   << Intra                        
         << " Atoms "   << Atoms                        
         << " SingleBond " << SBond                   
         << endl; 
}

string Scoring::str() const { 
    ostringstream so; 
    so << setw(10) << Name      << "\t" 
       << setw(15) << Fitness   << "\t"
       << setw(10) << HB_VDW    << "\t"
       << setw(5)  << Elect     << "\t"
       << setw(10) << Intra     << "\t"
       << setw(10) << Atoms     << "\t"
       << setw(5)  << SBond     << "\t"
       << setw(5)  << Run       << "\t"
       << setw(20) << FullName  
    ;
    return so.str(); 
}

bool fitness_sorter(const Scoring &a,const Scoring &b)
{
    return a.Fitness < b.Fitness;
}

#endif
