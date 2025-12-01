// author:      Chen Yu-Ju
// compiler:    g++ 3.3.3 (or later)
//
#ifndef FIOLIB_HEADER
#define FIOLIB_HEADER

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cctype>
#include <ctime>
#include <cerrno>

#include <dirent.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/stat.h>




//-----------------------------------------------------------------------------
namespace STRING {
    using namespace std;
    // file name/path/prefix/suffix
    string           getPrefix(const string& s);
    string           getSuffix(const string& s);
    string           RemovePath(const string& s);
    string           getNoPathPrefix(const string& s);
    template<typename T> string Num2Str(const T &);
    template<typename T> string itos(const T &number);

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
        {   // no '/' contained
            std::string::size_type pos2 = s.find('\\');
            if( pos2==std::string::npos ) {
                // no '\\' contained
                return s;
                //return s.substr(s.find_last_of('/')+1);
            } 
            return s.substr(s.find_last_of('\\')+1);
        }
        return s.substr(s.find_last_of('/')+1);
        
        //string::size_type pos = s.find('\\');
        //if( pos == string::npos )
        //    return s.substr(s.find_last_of(Slash)+1);
        //return s.substr(s.find_last_of('\\')+1);
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
//-----------------------------------------------------------------------------



namespace FioLib {
    using namespace std;
    enum Msg { show, hide } ; 
    static const char Slash='/';

    //file and dir issue
    inline bool             isFILE(const string&,Msg msg=hide);
    inline bool             isDIR(const string& dir,Msg msg=hide);
    inline void             createDIR(const string&);
    inline DIR*             createAndOpenDIR(const string&);
    inline DIR*             OpenDIR(const string&);
    inline string           getCurrentDir();
    inline string           checkDirRear(const string&);
    //inline char*            checkDirRear(char*);
    inline string           isLegalFilePath(const string& s);
    
    //issues of protein/drug format 
    inline bool             isLegalDrug(const string&, Msg msg=hide);
    inline bool             isLegalCav(const string& );
    inline const string     getFormat(const string&);
    inline int              getFormatID(const string&);
    inline bool             isPDB(const string&);
    inline bool             isMO2(const string&);
    inline bool             isMOL(const string&);
    inline bool             isLST(const string&,Msg msg=hide);


    // Time problem
    inline time_t           getFileDate(string file);
    inline string           getTimeFormatted(time_t);
    inline string           printTime(time_t time_predict);
    inline bool             removeFile(const string& files);
    
    
    // Calculating molecular weight of compound
    typedef map<string,double> Mass;
    inline void             setAtomMass();
    inline double           getMolecularWeight(string);

    // Others
    inline bool             isPBSexist(bool show_msg=false);
    inline string           isLegalAlpha(const string);
    inline int              getDigit(int);
    inline void             justSleep(unsigned int mseconds);

 // end of namespace definition


//file and dir issue
    inline bool isFILE(const string& s,Msg msg)
    {
        bool return_value = true;
        ifstream fin(s.c_str());
        if( not fin.good() ) 
        {
            if(msg==show)
                cerr << string("Cannot open file: "+ s) << endl;
            fin.clear();
            return false;
        }
        fin.close();
    
        DIR* dir=opendir( s.c_str() );
        if( dir == NULL ) 
            return_value = true;
        else {
            return_value = false;
            closedir(dir);
        }
    
        return return_value;
    }
    
    inline bool isDIR(const string& s,Msg msg)
    {
        if (s.size() <= 0 ) return false;
        bool return_value = true;
        DIR* dir=opendir( s.c_str() );
        if( dir == NULL ) 
        {
            if(msg==show)
                cerr << string("Cannot open file: "+ s) << endl;
            return_value = false;
        }
        else {
            return_value = true;
            closedir(dir);
        }
    
        return return_value;
    }


    inline void createDIR(const string& s)
    {
        if(s.size() <= 0 ) return ;
        DIR* dir=opendir( s.c_str() ); 
        if( dir != NULL ) {
            closedir(dir);
            return ;
        }
        
        string ss = "\"" + s + "\"";
        switch(errno)
        {
            case ENOENT:
                //cout << "creating dir " + ss  << endl;
                if( system( string("mkdir "+ss).c_str()) != 0 ) {
                //if( mkdir( s.c_str() ) != 0 ) {
                    cerr << string("mkdir "+s+" failed!") << endl;
                    exit(1);
                }
            break;
    
            case EACCES: 
            case ENOTDIR: 
            case ENAMETOOLONG:
                cerr<<(string("directory open failed: " + s))<<endl;
                cerr.flush();
                exit(1);
        }
    }
    
    inline DIR* createAndOpenDIR(const string& s)
    {
        if(s.size() <= 0 ) return NULL;
        static DIR* dir=opendir( s.c_str() ); 
        if( dir != NULL ) 
            return dir;
        
        string ss = '\"' + s + '\"';
        switch(errno)
        {
            case ENOENT:
                //cout << "creating dir " + ss  << endl;
                if( system( string("mkdir "+ss).c_str()) != 0 ) {
                    cerr << string("Fetal Error: mkdir "+ss+" failed!") << endl;
                    exit(1);
                }
                dir = opendir(s.c_str());
            break;
    
            case EACCES: 
            case ENOTDIR: 
            case ENAMETOOLONG:
                cerr<<(string("Directory: " + s + "open failed"))<<endl;
                cerr.flush();
                exit(1);
        }
        return dir;
    }

    inline DIR* OpenDIR(const string& s)
    {
        if(s.size() <= 0 ) return NULL;

        DIR* dir=opendir( s.c_str() ); 
        if( dir == NULL ) {
            cerr << "Error: Cannot open directory " + s << endl;
            exit(-1);
        }
        return dir;
    }

    //inline char* checkDirRear(char* s)
    //{
    //    int size = strlen(s);
    //    if( size < 256-1 ) {
    //        if( s[size-1] != Slash ) {
    //            s[size++] = Slash;
    //            s[size++] = '\0';
    //            //#ifdef WIN32
    //            //if( size < 256-1 ) {
    //            //    s[size++] = Slash;
    //            //    s[size] = '\0';
    //            //}
    //            //#endif
    //        }
    //    }
    //    return s;
    //}

    inline string checkDirRear(const string& ss)
    {
        string s(ss);
        int size=(int)s.size();
        if( size <= 0 ) return s;
        if( s[size-1] == '/' or s[size-1] == '\\' ) return s;
        if( s[size-2] == '/' or s[size-2] == '\\' ) return s;
        if( s[size-1] != '/' or s[size-1] != '\\')  s += '/';
        return s;
    }

    inline string createAndCheckDir(const string& s)
    {
        createDIR(s);
        return checkDirRear(s);
    }


    inline string getCurrentDir()
    {
        const size_t len=1024;
        char bf[len];
        string s = getcwd(bf,len);
        for(size_t i=0;i<s.size();i++)
            if(s[i]=='\\') s[i]=Slash;
        return checkDirRear(s);
        //return (s);
    }
    
    inline string isLegalFilePath(const string& s)
    {
        string ret(s);
        int trimed_size=s.size();
        for(int i=ret.size()-1;i>0;i--)
        {
            if(isspace(ret[i]) or iscntrl(ret[i]) ) 
            {
                --trimed_size;
            }
            else break;
        }
        return ret.substr(0,trimed_size);
    }

    

//issues of protein/drug format 
    inline bool isLegalDrug(const string& file, Msg msg)
    {
        if( isMOL(file) or isMO2(file) ) return true;
        if( isPDB(file) )
        {
            if(msg==show) {
                cerr << "Warning: Drug file is in PDB format" << endl;
                cerr.flush();
            }
            return true;
        }
        return false;
    }

    inline bool isLegalCav(const string& file)
    {
        if( isPDB(file) ) return true;
        return false;
    }
    
    inline const string getFormat(const string& file)
    {
        static int FormatID = getFormatID(file);
        switch(FormatID)
        {
            case 1: return "PDB";
            case 2: return "MOL";
            case 3: return "MOL2";
            case 9: return "List";
            default : break;
        }
        return "";
    }

    inline int getFormatID(const string& file) 
    {
        int FormatID = -1;         //unrecognized format
             if( isPDB(file) ) FormatID = 1;
        else if( isMOL(file) ) FormatID = 2;
        else if( isMO2(file) ) FormatID = 3;
        else if(      isLST(file) ) FormatID = 9;
        return FormatID;
    }
    
    inline bool isLST(const string& file,Msg msg)
    {
        //cout << "Check whether the file is List or not" << endl;

        //if(isLegalDrug(file)) return false;
        ifstream fin(file.c_str());
        if( not fin.good() )
        {
            if( msg )
                cerr << "Cannot open file/list: " << file << endl;
            return false;
        }
        int  line = 0;
        vector<string> errfile;
        int  error_count = 0;
        bool flag = true;
        string buffer;
        while( getline(fin,buffer) && ++line )
        {
            if( buffer.size()>0 ) 
            if( buffer[0] != '#' ) {
                ifstream fin2( buffer.c_str() );

                if( not fin2.good() ) {                    
                    string buffer2 = string("..") + Slash + "Drug" + Slash + buffer;
                    fin2.clear();
                    fin2.open( buffer2.c_str() );
                    if( not fin2.good() ) {
                        buffer2 = string("..") + Slash + "Ligand" + Slash + buffer;
                        
                        fin2.clear();
                        fin2.open(buffer2.c_str());
                        if( not fin2.good() ) {
                            if( msg==show ) {
                                errfile.push_back( buffer );
                            }
                            flag = false;
                            if(++error_count>2) {
                                cerr << "Error: " << file << " is not a correct list file" << endl;
                                cerr << "An accetpable list file MUST contain drug location in each line in it" << endl;
                                cerr.flush();
                                return false;
                            }
                        }
                    }
                }
                fin2.close();
            }
        }
        fin.close(); 
        fin.clear(); 
        
        return flag;
    }
    
    
    inline bool isPDB(const string& file) 
    {
        ifstream fin(file.c_str());
        if( not fin.good() )
        {
            return false;
            //cerr << "Cannot open pdb file: " << file << endl;
            //fin.clear();
            //string traditional_WCav="../WCavPDB/" + file;
            //fin.open( traditional_WCav.c_str() );
            //if( not fin.good() ) 
            //{
            //    fin.clear();
            //    string traditional_Ligand="../Ligand/" + file;
            //    fin.open( traditional_Ligand.c_str() );
            //    if( not fin.good() )
            //        return false;
            //}
        }
    
        string buffer;
        int num_HETATM=0;
        int num_CONNECT=0;
        //int num_SINGLE=0;
        int num_ATOM=0;
        while(getline(fin,buffer)) 
        {
            if(buffer.size()>=54 && num_HETATM<3) 
            {
                if(buffer.substr(0,6)=="HETATM")
                {
                    bool have_atom=false;
                    string atom=buffer.substr(12,4);
                    for(size_t i=0;i<atom.size();i++) 
                    {
                        if( isalpha( atom[i] ) ) 
                        { 
                            have_atom = true;
                            //break;
                        }
                    }
                    
                    double x = atof( buffer.substr(30,8).c_str() );
                    double y = atof( buffer.substr(38,8).c_str() );
                    double z = atof( buffer.substr(46,8).c_str() );
                    if( //have_atom and
                        x != 0.0 and
                        y != 0.0 and
                        z != 0.0    ) ++num_HETATM;
                }
                if(buffer.substr(0,4)=="ATOM")
                {
                    if( atof( buffer.substr(30,8).c_str() )!= 0.0 and
                        atof( buffer.substr(38,8).c_str() )!= 0.0 and
                        atof( buffer.substr(46,8).c_str() )!= 0.0    ) ++num_ATOM;
                }
            }
            getline(fin,buffer);
            if( buffer.size()>16 ) {
                if(buffer.substr(0,6)=="CONECT")
                {   
                    if( atoi( buffer.substr(6,5).c_str() )!= 0 and
                        atoi( buffer.substr(11,5).c_str() )!= 0 ) ++num_CONNECT;
                    if( num_CONNECT > 1 ) break;
                }
                if( (num_HETATM>2 && num_CONNECT>1 ) || num_ATOM>2 ) break;
            }
            
            //if( buffer.size()>=16 ) {
            //    if(buffer.substr(0,6)=="SINGLE")
            //    {   
            //        if( atoi( buffer.substr(6,5).c_str() )!= 0.0 and
            //            atoi( buffer.substr(11,5).c_str() )!= 0.0 ) ++num_SINGLE;
            //        if( num_SINGLE > 1 ) break;
            //    }
            //}
            
        }
        fin.close();
        if( (num_HETATM>1 or num_CONNECT>=1 ) or num_ATOM>2 ) 
            return true;
        else return false;
    }
    
    inline bool isMOL(const string& file)
    {
        //cout << "Check whether the file is Mol or not" << endl;
        ifstream fin(file.c_str());
        if( not fin.good() )
        {
            return false;
            //fin.clear();
            //string traditional_Drug="../Drug/" + file;
            //fin.open( traditional_Drug.c_str() );
            //if( not fin.good() ) 
            //    return false;
            //cerr << "Cannot open mol file: " << file << endl;
        }
        string buffer;
        bool have_header=false;
        bool correct_atom_num=false;
        while(getline(fin,buffer)) 
        {
            if(buffer.size()>=6) 
            {
                int num_atoms=atoi(buffer.substr(0,3).c_str());
                int num_bonds=atoi(buffer.substr(3,3).c_str());

                if(0) {
                    cout << "num_bonds = " << num_bonds << " - " << buffer.substr(3,3) << endl;
                    cout << "num_atoms = " << num_atoms << " - " << buffer.substr(0,3) << endl;
                    getchar();
                }

                if( num_atoms>0 && num_bonds>0 ) 
                { 
                    have_header=true;
                    getline(fin,buffer);
                    if(buffer.size()>33) 
                    {
                        if( atof(buffer.substr(0,10).c_str())!=0.0  and
                            atof(buffer.substr(10,10).c_str())!=0.0 and
                            atof(buffer.substr(20,10).c_str())!=0.0  ) 
                        {
                                correct_atom_num=true;
                                break;
                        }
                        
                    }
                }
            }
        }
        fin.close();
        if( have_header && correct_atom_num) 
            return true;
        else return false;
    }
    
    
    inline bool isMO2(const string& file)
    {
        //cout << "Check whether the file is Mol2 or not" << endl;
        bool have_atom_header=false,have_bond_header=false;
        ifstream fin( file.c_str() );
        if( not fin.good() )
        {
            return false;
            //fin.clear();
            //string traditional_Drug="../Drug/" + file;
            //fin.open( traditional_Drug.c_str() );
            //if( not fin.good() ) 
            //    return false;
            //cerr << "Cannot open mol2 file: " << file << endl;
        }
        string buffer;
        while(getline(fin,buffer)) 
        {
            if(buffer.size()>=13)
            {
                if(buffer.substr(0,13)=="@<TRIPOS>ATOM") have_atom_header=true;
                if(buffer.substr(0,13)=="@<TRIPOS>BOND") have_bond_header=true;
                if( have_atom_header && have_bond_header ) break;
            }
        }
        fin.close();
        if( have_atom_header && have_bond_header ) 
            return true;
        else return false;
    }

// implementation for MW
    inline void setAtomMass(Mass& atomic_mass)
    {
        if( atomic_mass.size() ) return ;
        atomic_mass["H"] = 1.0079 ; 
        atomic_mass["He"] = 4.0026 ;
        atomic_mass["Li"] = 6.9410 ;
        atomic_mass["Be"] = 9.0122 ;
        atomic_mass["B"] = 10.8110;
        atomic_mass["C"] = 12.0107;
        atomic_mass["N"] = 14.0067;
        atomic_mass["O"] = 15.9994;
        atomic_mass["F"] = 18.9984;
        atomic_mass["Ne"] = 20.1797;
        atomic_mass["Na"] = 22.9898;
        atomic_mass["Mg"] = 24.3050;
        atomic_mass["Al"] = 26.9815;
        atomic_mass["Si"] = 28.0855;
        atomic_mass["P"] = 30.9738;
        atomic_mass["S"] = 32.0660;
        atomic_mass["Cl"] = 35.4527;
        atomic_mass["Ar"] = 39.9480;
        atomic_mass["K"] = 39.0983;
        atomic_mass["Ca"] = 40.0780;
        atomic_mass["Sc"] = 44.9559;
        atomic_mass["Ti"] = 47.9670;
        atomic_mass["V"] = 50.9415;
        atomic_mass["Cr"] = 51.9961;
        atomic_mass["Mn"] = 54.9380;
        atomic_mass["Fe"] = 55.8450;
        atomic_mass["Co"] = 58.9332;
        atomic_mass["Ni"] = 58.6934;
        atomic_mass["Cu"] = 63.5460;
        atomic_mass["Zn"] = 65.3900;
        atomic_mass["Ga"] = 69.6230;
        atomic_mass["Ge"] = 72.6100;
        atomic_mass["As"] = 74.9216;
        atomic_mass["Se"] = 78.9600;
        atomic_mass["Br"] = 79.9040;
        atomic_mass["Kr"] = 83.8000;
        atomic_mass["Rb"] = 85.4678;
        atomic_mass["Sr"] = 87.6200;
        atomic_mass["Ag"] = 107.8682;
        atomic_mass["I"] = 126.9045;
        return ;
    }
    
    inline double getMolecularWeight(string file)
    {
        static Mass mw_table;
        setAtomMass(mw_table);

        double molecular_weight=0;
        ifstream fin( file.c_str() );
        string parsing_buffer;
        if(not fin.good())
        {
            cerr << "Error: Cannot open file " << file << endl;
            return -1;
        }
    
        if(isMOL(file))
        {
            while( getline(fin,parsing_buffer) )
            if(parsing_buffer.size()>40)
            if(isalpha(parsing_buffer[31]))
            {
                istringstream sin(parsing_buffer.substr(31,3));
                string atom;  sin >> atom;
                molecular_weight += mw_table[ atom ];
            }
        }
        else if(isMO2(file))
        {
            while( getline(fin,parsing_buffer) )
            {
                if(parsing_buffer.size()>=13)
                if(parsing_buffer.substr(0,13)=="@<TRIPOS>ATOM")
                break;
            }
            while( getline(fin,parsing_buffer) )
            if(parsing_buffer.size()>50)
            {
                istringstream sin(parsing_buffer);
                string atom;  sin >> atom >> atom; 
                size_t pos=1;
                for(size_t i=0;i<atom.size();i++)
                {
                    if(atom[i]==' ' or isdigit(atom[i]))
                    {
                        pos = i;
                        break;
                    }
                }
                atom = atom.substr(0,pos);
                molecular_weight += mw_table[ atom ];
            }
        }
        else if(isPDB(file))
        {
            while( getline(fin,parsing_buffer) )
            if(parsing_buffer.size()>18)
            if(parsing_buffer.substr(0,6)=="HETATM")
            {
                istringstream sin(parsing_buffer.substr(12,5));
                string atom;  sin >> atom; 
                size_t pos=1;
                for(size_t i=0;i<atom.size();i++)
                {
                    if(atom[i]==' ' or isdigit(atom[i]))
                    {
                        pos = i;
                        break;
                    }
                }
                atom = atom.substr(0,pos);
                molecular_weight += mw_table[ atom ];
            }
        }
        else 
        {
            cerr << "Error: " << file << " - Unrecognized compound format (only for pdb,mol,mol2)" << endl;
            molecular_weight = -1;
        }
        fin.close();
        return molecular_weight;
    }


//Others
    inline int SystemCall(const string& cmd)
    {
        #ifdef WIN32
        return system(cmd.c_str());
        #else
        string cmd_sh("sh " + cmd );
        return system(cmd_sh.c_str());
        #endif
    }

    inline bool copyFile(const string& src, const string& des)
    {
        if(src==des) return true;
        ifstream fin(src.c_str());
        ofstream fout(des.c_str());
        if( not fin.good() ) return false;
        if( not fout.good() ) return false;
        string buffer;
        while(getline(fin,buffer))
            fout << buffer << endl;
        fin.close();
        fout.close();
        return true;
    }

    inline bool removeFile(const string& files)
    {
        return (remove(files.c_str())==0);
    }

    inline bool isPBSexist(bool show_msg)
    {
        if(show_msg) cout << "Checking for qstat of PBS... " ;
        int qstat= system("qstat > .check_pbs_status 2>&1");
    
        removeFile(".check_pbs_status");
        
        if( show_msg && qstat )
        {
            cout << "no!" << endl<<endl;
            cerr << "Your system does NOT support PBS" << endl; 
            cerr << "Please check whether PBS be installed or not" << endl; 
            cerr << "If you have installed PBS, please set PBS binaries in your $PATH" << endl;
            cerr << "Try to run the program again and do NOT use -j option" << endl;
            cerr << "( force job number to be 1 ) " << endl;
            exit(0);    
        }
        if(show_msg) cout << "available!" << endl;
        if(qstat) return false;
        return true;
    }

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
    
    inline string isLegalAlpha(const string& s)
    {
        string ret;
        string legal_str = ",.`~!@#$%^&()_+-='[]{}";
        bool special_legal_char;
        for(size_t i=0;i<s.size();i++)
        {
            special_legal_char=false;
            for(size_t j=0;j<legal_str.size();j++)
                if(s[i]==legal_str[j]) 
                    special_legal_char=true;
    
            if( isalnum(s[i]) or s[i]==Slash or special_legal_char)
            {
    
                ret.push_back(s[i]);
            }
            else break;
        }
        return ret;
    }
    
    inline int getDigit(int num)
    {
        int Digit = 0;
        int temp = num;
        while( temp > 0 )
        {
            temp /= 10;
            Digit++;
        }
        return Digit;
    }

    inline void justSleep(unsigned int mseconds)
    {
        clock_t goal = mseconds + clock();
        while( goal > clock() );
    }

    inline void justSleepSecond(unsigned int mseconds)
    {
        clock_t goal = 1000*mseconds + clock();
        while( goal > clock() );
    }


// time and date
    inline time_t getFileDate(string file)
    {
        struct stat st;
        stat(file.c_str(),&st);
        return st.st_mtime;
    }
    
    inline string getTimeFormatted(time_t rawtime)
    {
        struct tm * timeinfo;
        //time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        ostringstream os;
        os  << setw(2) << setfill('0') << timeinfo->tm_hour << ":"
            << setw(2) << setfill('0') << timeinfo->tm_min  << ":"
            << setw(2) << setfill('0') << timeinfo->tm_sec  << ", "
            << setw(2) << setfill('0') << timeinfo->tm_mon+1<< "/"
            << setw(2) << setfill('0') << timeinfo->tm_mday << ", "
            << setw(2) << setfill('0') << timeinfo->tm_year+1900;
        return string(os.str());
    }
    
    inline string printTime(time_t time_predict)
    {
        time_t min_left=0;
        time_t hr_left=0;
        time_t day_left=0;
        if( time_predict > 60 ) 
        {
            min_left = time_predict / 60;
            if( min_left > 60 )
            {
                hr_left = min_left / 60;
                min_left = min_left % 60;
                if( hr_left > 24 ) 
                {
                    day_left = hr_left / 24;
                    hr_left = hr_left % 24;
                }
            }
        }
        ostringstream oss;
        oss << day_left << "day " << hr_left << "hr " << min_left << "min" ;
        return oss.str();
    }
}



namespace FStat {
    using namespace std;

    const double O_O=0.000001;

    template<typename T>
    inline T sum_of_square(const vector<T> &v) 
    {
        T sum = 0;
        for(size_t i=0; i<v.size(); ++i) 
            sum +=  (v[i]*v[i]);
        return sum;
    }

    template<typename T> 
    inline T Summation(const vector<T> &v) 
    {
        T sum = 0;
        for(size_t i=0; i<v.size(); ++i) 
        {
            sum += v[i];
        }
        return sum;
    }


    template<typename T> 
    inline T Average(const vector<T> &v) 
    {
        T sum = 0;
        for(size_t i=0; i<v.size(); ++i) 
        {
            sum += v[i];
        }
        return sum / v.size();
    }

    template<typename T>
    inline T Maximum(const vector<T>& v) 
    {
        if( v.size() == 0 ) return 0;
        T max = v[0];
        for(size_t i=0; i<v.size(); ++i) {
        if( max > v[i] ) {
            max = v[i];
        }
        }
        return max;
    }


    template<typename T>
    inline T Minimum(const vector<T>& v) 
    {
        if( v.size() == 0 ) return 0;
        T min = v[0];
        for(size_t i=0; i<v.size(); ++i) {
            if( min < v[i] ) {
                min = v[i];
            }
        }
        return min;
    }

    template<typename T> 
    inline T Stdev(const vector<T>& v,const int freedom=1) 
    {
        if( v.size() <=1 ) return O_O;
        T avg = Average( v );
        T ss = 0;
        for(size_t i=0;i<v.size();++i)
            ss += ((v[i]-avg)*(v[i]-avg));
        return sqrt( ss / (v.size()-freedom) );
    }

    template<typename T>
    inline size_t check_width(const vector<T>& v1,const vector<T>& v2) 
    {
        return ( v1.size()==v2.size() );
    }

    template<typename T>
    inline T Pearson(const vector<T>& v1,const vector<T>& v2) 
    {
        if( ! check_width( v1,v2) ) 
        { 
            cerr<< "data is not consistent" << endl; 
            exit(0); 
        }
        T N = v1.size();
        T SumX =  Summation(v1);
        T SumY =  Summation(v2);
        T SsqX =  sum_of_square( v1 );
        T SsqY =  sum_of_square( v2 );
        double XiYi = 0.0;
        for(size_t i=0;i<N;++i) 
            XiYi += (v1[i])*(v2[i]);
        double ret = (N*XiYi) - (SumX*SumY);
        double dom = sqrt( (N*SsqX-SumX*SumX) *(N*SsqY-SumY*SumY) );
        ret /= dom;

        return 1-ret;
    }

    inline double Pearson(const vector<double>& v1,const vector<double>& v2, 
            const double& SumX, const double& SumY, 
            const double& SsqX, const double& SsqY) 
    {
        /*
        cout << "N =\t" << v1.size() << endl;
        cout << "SumX=\t" << SumX << endl;
        cout << "SumY=\t" << SumY << endl;
        cout << "SsqX=\t" << SsqX << endl;
        cout << "SsqY=\t" << SsqY << endl;
        cout << "Vec1=\t" << v1 << endl;
        cout << "Vec2=\t" << v2 << endl;
        */
        if( ! check_width( v1,v2) ) 
        { 
            cerr<< "data is not consistent" << endl; 
            exit(0); 
        }
        double N = v1.size();
        double XiYi = 0.0;
        for(size_t i=0;i<N;++i) 
            XiYi += (v1[i])*(v2[i]);
        
        double ret = (N*XiYi) - (SumX*SumY);
        double dom = sqrt( (N*SsqX-SumX*SumX) *(N*SsqY-SumY*SumY) );
        //cout << "XiYi=\t" << XiYi << endl;
        //cout << "ret =\t"  << ret << endl;
        //cout << "dom =\t"  << dom << endl;
        ret /= dom;
        //cout << "Pearsion CC =\t"  << ret << endl;
        return 1-ret;
    }


    template<typename T>
    inline T sqEuclidean(const vector<T>& v1,const vector<T>& v2) 
    {
        if( ! check_width( v1,v2) ) 
        { 
            cerr<< "Error: data is not consistent" << endl; 
            exit(0); 
        }   
        T sum = 0;
        for( size_t i=0; i<v1.size(); i++) 
        {
            sum += (v1[i]-v2[i])*(v1[i]-v2[i]);
        }
        return ( sum );
    }

    template<typename T>
    inline T Euclidean(const vector<T>& v1,const vector<T>& v2) 
    {
        if( ! check_width( v1,v2) ) 
        { 
            cerr<< "Error: data is not consistent" << endl; 
            exit(0); 
        }   
        T sum = 0;
        for( size_t i=0; i<v1.size(); i++) 
        {
            sum += (v1[i]-v2[i])*(v1[i]-v2[i]);
        }
        return sqrt( sum );
    }

    template<typename T>
    inline T Hamming(const vector<T>& v1,const vector<T>& v2) 
    {
        if( ! check_width( v1,v2) ) 
        { 
            cerr<< "Error: data is not consistent" << endl; 
            exit(0); 
        }   
        T sum = 0;
        for( size_t i=0; i<v1.size(); i++) 
        {
            if(v1[i]!=v2[i]) sum++;
        }
        //cout << "v1 = " << v1 << endl;
        //cout << "v2 = " << v2 << endl;
        //cout << "dist = " << sum << endl;
        //getchar();
        return ( sum );
    }


    template<typename T> 
    double Tanimoto(const vector<T> &s, const vector<T> &t) {
        if( s.size() != t.size() ) {
            cerr << "Error: no match type, cannot evaluate distance " << s.size() << ","  <<  t.size() << endl;
            exit(1);
        }
        double a=0.0, b=0.0, c=0.0;
        for(size_t i=0;i<s.size();i++) {
            
            if(s[i]) a++;
            if(t[i]) b++;
            if(s[i]*t[i]) c++;
        }
        double tanimoto = c / (a+b-c);
        return 1 - tanimoto;
    }


    template<typename T> 
    double TanimotoCont(vector<T> &s,vector<T> &t) {
        if( s.size() != t.size() ) {
            cerr << "Error: no match type, cannot evaluate distance " << s.size() << ","  <<  t.size() << endl;
            exit(1);
        }

        double a=0,b=0,c=0;
        for(size_t i=0;i<s.size();i++) {
            a += s[i];
            b += t[i];
            c += s[i]*t[i];
        }
        double up = c;
        double down = a+b-c;

        return up/down;
    }

} // end of namespace Statistics

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

#endif
