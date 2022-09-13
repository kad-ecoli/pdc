/* parse PDB file into data structure similar to Bio.PDB in biopython
 * (model - chain - residue - atom). */
#ifndef PDBParser_HPP
#define PDBParser_HPP 1

#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <cstdint>
#include "StringTools.hpp"
#include "pstream.h"

using namespace std;

struct AtomUnit    // struct for each atom entry
{
    string name;       // atom name
    vector<int32_t> xyz; // coordinate
    int16_t bfactor;     // temperature factor
};

struct ResidueUnit // struct for each residue
{
    int resi;               // residue sequence number
    char icode;             // insertion code
    string resn;            // residue name
    vector<AtomUnit> atoms; // list of atoms
};

struct ChainUnit  // struct for each chain
{
    string chainID_full;          // chain ID, might be more than 1 char
    char chainID;                 // short chain ID, must be 1 char
    string sequence;              // sequence converted from CA coordinate
    vector<ResidueUnit> residues; // list of residues
};

struct ModelUnit  // struct for each model in mult-model PDB
{
    vector<ChainUnit> chains; // list of chains
};

int32_t XYZtoint32(string line)
{
    float  xf=atof(line.c_str());
    int32_t x=atoi(line.substr(0,4).c_str())*1000;
    if (line[5]!='0') if (xf>0) x+=(int)(line[5]-'0')*100;
                      else      x-=(int)(line[5]-'0')*100;
    if (line[6]!='0') if (xf>0) x+=(int)(line[6]-'0')*10;
                      else      x-=(int)(line[6]-'0')*10;
    if (line[7]!='0') if (xf>0) x+=(int)(line[7]-'0');
                      else      x-=(int)(line[7]-'0');
    return x;
}

int32_t Btoint16(string line)
{
    float  xf=atof(line.c_str());
    int16_t x=atoi(line.substr(0,3).c_str())*100;
    if (line[4]!='0') if (xf>0) x+=(int)(line[4]-'0')*10;
                      else      x-=(int)(line[4]-'0')*10;
    if (line[5]!='0') if (xf>0) x+=(int)(line[5]-'0');
                      else      x-=(int)(line[5]-'0');
    return x;
}

/* parse one line in PDB file, append the data to pep. 
 * used by read_pdb_structure
 * allowX: 0 - ATOM, 1 - ATOM and MSE, converting MSE to MET
 * return 0 if line not parsed, 1 if line is parsed, 
 * 2 if line should be parsed by other subrountine
 */
int parse_pdb_line(const string line,ModelUnit &pep, ChainUnit &chain,
    ResidueUnit &residue, AtomUnit &atom,
    const int atomic_detail=2,const int allowX=1)
{
    if (StartsWith(line,"HEADER") || StartsWith(line,"TITLE ") ||
        StartsWith(line,"COMPND") || StartsWith(line,"SOURCE") ||
        StartsWith(line,"DBREF")) return 2;
    string record_name=line.substr(0,6);
    char altLoc=line[16];

    atom.name=line.substr(12,4);
    residue.resn=line.substr(17,3);

    if ((allowX==0 && record_name!="ATOM  ")||
        (allowX==1 && record_name!="ATOM  " &&  
            !(record_name=="HETATM" && residue.resn=="MSE"))||
        (allowX>=2 && record_name!="ATOM  " && record_name!="HETATM"))
        return 0;
    
    // ignore alternatively locating residues
    if (altLoc!=' ' && altLoc!='A') return 0;

    if ((atomic_detail==0 && atom.name!=" CA " && atom.name!=" C3'")||
        (atomic_detail==1 && atom.name!=" CA " && atom.name!=" C3'" &&
         atom.name!=" N  "&& atom.name!=" C  " && atom.name!=" O  ")) return 0;
    
    if (residue.resn=="MSE" && allowX<3)
    {
        record_name="ATOM  ";
        residue.resn="MET";
    }
    if (record_name=="HETATM") return 0;
    
    if      (atom.name=="SE  ") atom.name=" SD ";
    else if (atom.name==" O1P") atom.name=" OP1";
    else if (atom.name==" O2P") atom.name=" OP2";
    else if (atom.name[3]=='*') atom.name=atom.name.substr(0,3)+"'";
    chain.chainID=line[21];
    if (chain.chainID==' ') chain.chainID='_';
    residue.resi=atoi(line.substr(22,4).c_str());
    residue.icode=line[26];
    //atom.xyz[0]=1000*atof(line.substr(30,8).c_str());
    //atom.xyz[1]=1000*atof(line.substr(38,8).c_str());
    //atom.xyz[2]=1000*atof(line.substr(46,8).c_str());
    atom.xyz[0]=XYZtoint32(line.substr(30,8));
    atom.xyz[1]=XYZtoint32(line.substr(38,8));
    atom.xyz[2]=XYZtoint32(line.substr(46,8));
    atom.bfactor=0;
    //if (line.size()>=66) atom.bfactor=100*atof(line.substr(60,6).c_str());
    if (line.size()>=66) atom.bfactor=Btoint16(line.substr(60,6).c_str());

    int chain_index=-1;
    for (int c=0;c<pep.chains.size();c++)
        if (pep.chains[c].chainID==chain.chainID) chain_index=c;
    if (chain_index==-1)
    {
        pep.chains.push_back(chain);
        chain_index=pep.chains.size()-1;
    }

    if (pep.chains[chain_index].residues.size()==0||
        pep.chains[chain_index].residues.back().resi !=residue.resi||
        pep.chains[chain_index].residues.back().icode!=residue.icode)
        pep.chains[chain_index].residues.push_back(residue);
    
    pep.chains[chain_index].residues.back().atoms.push_back(atom);
    return 1;
}

/* atomic_detail: 0 - CA only, 1 - backbone heavy atoms (CA C N O), 2 - all atom
 * allowX: 0 - ATOM, 1 - ATOM and MSE, converting MSE to MET, 
 *         2 - all, converting MSE to MET, 3 - all, no conversion
 * filename: full filename path, stdin if filename=="-"
 */
ModelUnit read_pdb_structure(const char *filename,string &header,
    const int atomic_detail=2,const int allowX=1)
{
    ModelUnit pep;

    string line="";
    string record_name="ATOM  ";
    char altLoc=' ';

    AtomUnit atom;
    atom.xyz.assign(3,0);

    ResidueUnit residue;

    ChainUnit chain;

    string filename_str=(string) filename;
    
    int use_stdin=(filename_str=="-");
    int use_pstream=0; // input is compressed

    ifstream fp;
    redi::ipstream fp_gz; // if file is compressed
    if (filename_str.length()>=3 && 
        filename_str.substr(filename_str.length()-3,3)==".gz")
    {
        // gzip pdb
        fp_gz.open("zcat "+filename_str);
        use_pstream=1;
    }
    else
    {
        fp.open(filename,ios::in); //ifstream fp(filename,ios::in);
    }

    while(use_stdin?cin.good():(use_pstream?fp_gz.good():fp.good()))
    {
        if (use_stdin)
            getline(cin,line);
        else if (use_pstream)
            getline(fp_gz,line);
        else
            getline(fp,line);

        if (line.substr(0,3)=="END") break;
        if (line.length()<53) continue;
        
        if (parse_pdb_line(line,pep,chain,residue,atom,atomic_detail,allowX)==2)
            header+=rstrip(line)+'\n';
    }
    if (!use_stdin)
    {
        if (use_pstream==0)
            fp.close();
        else
            fp_gz.close();
    }
    
    chain.residues.clear();
    residue.atoms.clear();
    return pep;
}

/* i - first atom index */
string write_pdb_structure(ChainUnit &chain,int &i)
{
    stringstream buf;
    int r,a;
    char chainID=' ';
    int32_t x,y,z;

    for (r=0;r<chain.residues.size();r++)
    {
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            chainID=chain.chainID;
            if (chainID=='.' || chainID=='_') chainID=' ';
            x=chain.residues[r].atoms[a].xyz[0];
            y=chain.residues[r].atoms[a].xyz[1];
            z=chain.residues[r].atoms[a].xyz[2];
            buf<<"ATOM  "
               <<resetiosflags(ios::left)<<setw(5)<<i++<<' '
               <<chain.residues[r].atoms[a].name<<' '
               <<chain.residues[r].resn<<' '<<chain.chainID<<setw(4)
               <<chain.residues[r].resi<<chain.residues[r].icode<<"   "
               <<setiosflags(ios::fixed)<<setprecision(3)
               <<setw(8)<<0.001*x<<setw(8)<<0.001*y<<setw(8)<<0.001*z
               <<"  1.00"<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)
               <<0.01*chain.residues[r].atoms[a].bfactor
               <<"           "
               <<Trim(chain.residues[r].atoms[a].name)[0]<<"  \n";
        }
    }
    r--;
    buf<<"TER   "
       <<resetiosflags(ios::left)<<setw(5)<<i++<<"      "
       <<chain.residues[r].resn<<' '<<chainID<<setw(4)
       <<chain.residues[r].resi<<chain.residues[r].icode
       <<"                                                     "
       <<endl;
    return buf.str();
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_pdb_structure(const char *filename,ChainUnit &chain)
{
    int i=1;
    if (strcmp(filename,"-")==0)
        cout<<write_pdb_structure(chain,i);
    else
    {
        ofstream fp(filename);
        fp<<write_pdb_structure(chain,i);
        fp.close();
    }
}

string write_pdb_structure(ModelUnit &pep,string &header)
{
    string txt;
    string line="";
    int c,r,s;
    vector<string> line_vec;
    Split(header,line_vec,'\n');
    for (s=0;s<line_vec.size();s++)
    {
        line=line_vec[s];
        if (line.size()<80) 
        {
            for (r=line.size();r<80;r++) line+=' ';
        }
        txt+=line+'\n';
    }
    vector<string> ().swap(line_vec);
    stringstream buf;
    line="";
    char chainID;
    for (c=0;c<pep.chains.size();c++)
    {
        s=0;
        chainID=pep.chains[c].chainID;
        if (chainID=='_' || chainID=='.') chainID=' ';
        for (r=0;r<pep.chains[c].residues.size();r++)
        {
            if (line.size()==0)
            {
                s++;
                buf<<"SEQRES "<<setw(3)<<s<<' '<<chainID<<' '
                    <<setw(4)<<pep.chains[c].residues.size()<<" "<<flush;
                line=buf.str();
                buf.str(string());
            }
            line+=" "+pep.chains[c].residues[r].resn;
            if (line.size()>=70)
            {
                txt+=line+"          \n";
                line.clear();
            }
            else if (r+1==pep.chains[c].residues.size())
            {
                while (line.size()<80) line+=' ';
                txt+=line+'\n';
                line.clear();
            }
        }
    }
    txt+=""
"CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          \n"
"ORIGX1      1.000000  0.000000  0.000000        0.00000                         \n"
"ORIGX2      0.000000  1.000000  0.000000        0.00000                         \n"
"ORIGX3      0.000000  0.000000  1.000000        0.00000                         \n"
"SCALE1      1.000000  0.000000  0.000000        0.00000                         \n"
"SCALE2      0.000000  1.000000  0.000000        0.00000                         \n"
"SCALE3      0.000000  0.000000  1.000000        0.00000                         \n"
"MODEL        1                                                                  \n";
    int i=1; // atom index
    for (c=0;c<pep.chains.size();c++)
        txt+=write_pdb_structure(pep.chains[c],i);
    txt+=""
"ENDMDL                                                                          \n"
"END                                                                             \n";
    return txt;
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_pdb_structure(const char *filename,ModelUnit &pep,string &header)
{
    if (strcmp(filename,"-")==0)
        cout<<write_pdb_structure(pep,header)<<flush;
    else
    {
        ofstream fp(filename);
        fp<<write_pdb_structure(pep,header)<<flush;
        fp.close();
    }
}

/* convert pdb structure to fasta sequence 
 * convertX - how to deal with non-standard amino acids
 *            0 - only 20 standard amino acids
 *            1 - 20 standard amino acids + MSE
 *            2 - non-standard amino acid with known parent, 
 *                all to legal amino acid in BLOSUM
 *            3 - non-standard amino acid with known parent
 */
inline char aa3to1(const string resn,const int convertX=2)
{
    // 20 standard amino acid + MSE
    if (resn[0]==' ' && (resn[1]=='D'||resn[1]==' ')) return tolower(resn[2]);
    if (resn=="ALA") return 'A';
    if (resn=="CYS") return 'C';
    if (resn=="ASP") return 'D';
    if (resn=="GLU") return 'E';
    if (resn=="PHE") return 'F';
    if (resn=="GLY") return 'G';
    if (resn=="HIS") return 'H';
    if (resn=="ILE") return 'I';
    if (resn=="LYS") return 'K';
    if (resn=="LEU") return 'L';
    if (resn=="MET") return 'M';
    if (resn=="ASN") return 'N';
    if (resn=="PRO") return 'P';
    if (resn=="GLN") return 'Q';
    if (resn=="ARG") return 'R';
    if (resn=="SER") return 'S';
    if (resn=="THR") return 'T';
    if (resn=="VAL") return 'V'; 
    if (resn=="TRP") return 'W';
    if (resn=="TYR") return 'Y';

    if (resn=="MSE" && convertX>=1) return 'M';

    if (convertX>=2)
    {
        // non-standard amino acid with known parent
        if (resn=="CHG"||resn=="HAC"||resn=="AYA"||resn=="TIH"||resn=="BNN"||
            resn=="ALM"||resn=="TPQ"||resn=="MAA"||resn=="PRR"||resn=="FLA"||
            resn=="AIB"||resn=="DAL"||resn=="CSD"||resn=="DHA"||resn=="DNP") 
            return 'A';
        if (resn=="PR3"||resn=="CCS"||resn=="C6C"||resn=="SMC"||resn=="BCS"||
            resn=="SCY"||resn=="DCY"||resn=="SCS"||resn=="CME"||resn=="CY1"||
            resn=="CYQ"||resn=="CEA"||resn=="CYG"||resn=="BUC"||resn=="PEC"||
            resn=="CYM"||resn=="CY3"||resn=="CSO"||resn=="SOC"||resn=="CSX"||
            resn=="CSW"||resn=="EFC"||resn=="CSP"||resn=="CSS"||resn=="SCH"||
            resn=="OCS"||resn=="SHC"||resn=="C5C") return 'C';
        if (resn=="DGL"||resn=="GGL"||resn=="CGU"||resn=="GMA"||resn=="5HP"||
            resn=="PCA") return 'E';
        if (resn=="ASQ"||resn=="ASB"||resn=="ASA"||resn=="ASK"||resn=="ASL"||
            resn=="2AS"||resn=="DAS"||resn=="DSP"||resn=="BHD") return 'D';
        if (resn=="PHI"||resn=="PHL"||resn=="DPN"||resn=="DAH"||resn=="HPQ")
            return 'F';
        if (resn=="GLZ"||resn=="SAR"||resn=="GSC"||resn=="GL3"||resn=="MSA"||
            resn=="MPQ"||resn=="NMC") return 'G';
        if (resn=="NEM"||resn=="NEP"||resn=="HSD"||resn=="HSP"||resn=="MHS"||
            resn=="3AH"||resn=="HIC"||resn=="HIP"||resn=="DHI"||resn=="HSE") 
            return 'H';
        if (resn=="IIL"||resn=="DIL") return 'I';
        if (resn=="DLY"||resn=="LYZ"||resn=="SHR"||resn=="ALY"||resn=="TRG"||
            resn=="LYM"||resn=="LLY"||resn=="KCX") return 'K';
        if (resn=="NLE"||resn=="CLE"||resn=="NLP"||resn=="DLE"||resn=="BUG"||
            resn=="NLN"||resn=="MLE") return 'L';
        if (resn=="FME"||resn=="CXM"||resn=="OMT") return 'M';
        if (resn=="MEN") return 'N';
        if (resn=="DPR"||resn=="HYP") return 'P';
        if (resn=="DGN") return 'Q';
        if (resn=="AGM"||resn=="ACL"||resn=="DAR"||resn=="HAR"||resn=="HMR"||
            resn=="ARM") return 'R';
        if (resn=="OAS"||resn=="MIS"||resn=="SAC"||resn=="SEL"||resn=="SVA"||
            resn=="SET"||resn=="DSN"||resn=="SEP") return 'S';
        if (resn=="DTH"||resn=="TPO"||resn=="ALO"||resn=="BMT") return 'T';
        if (resn=="DVA"||resn=="MVA"||resn=="DIV") return 'V';
        if (resn=="LTR"||resn=="DTR"||resn=="TRO"||resn=="TPL"||resn=="HTR") 
            return 'W';
        if (resn=="PAQ"||resn=="STY"||resn=="TYQ"||resn=="IYR"||resn=="TYY"||
            resn=="DTY"||resn=="TYB"||resn=="PTR"||resn=="TYS") return 'Y';
        
        // undeterminted amino acid
        if (resn=="ASX") return 'B'; // or D or N
        if (resn=="GLX") return 'Z'; // or Q or E

        // amino acid with no code in BLOSUM62
        if (convertX>=3)
        {
            if (resn=="SEC") return 'U';
            if (resn=="PYL") return 'O';
        }
        if (resn=="SEC") return 'C';
        if (resn=="PYL") return 'K';
    }
    return 'X';
}

/* only residues in 'ATOM' record with CA or C3' atoms are converted */
string pdb2fasta(ChainUnit& chain)
{
    chain.sequence="";
    int r,a;
    for (r=0;r<chain.residues.size();r++)
    {
         for (a=0;a<chain.residues[r].atoms.size();a++)
         {
             if (chain.residues[r].atoms[a].name==" CA "||
                 chain.residues[r].atoms[a].name==" C3'")
                 chain.sequence+=aa3to1(chain.residues[r].resn);
         }
    }
    return chain.sequence;
}

inline string aa1to3(const char aa)
{
    if (aa=='A') return "ALA";
    if (aa=='B') return "ASX";
    if (aa=='C') return "CYS";
    if (aa=='D') return "ASP";
    if (aa=='E') return "GLU";
    if (aa=='F') return "PHE";
    if (aa=='G') return "GLY";
    if (aa=='H') return "HIS";
    if (aa=='I') return "ILE";
    if (aa=='K') return "LYS";
    if (aa=='L') return "LEU";
    if (aa=='M') return "MET";
    if (aa=='N') return "ASN";
    if (aa=='O') return "PYL";
    if (aa=='P') return "PRO";
    if (aa=='Q') return "GLN";
    if (aa=='R') return "ARG";
    if (aa=='S') return "SER";
    if (aa=='T') return "THR";
    if (aa=='U') return "SEC"; 
    if (aa=='V') return "VAL"; 
    if (aa=='W') return "TRP";
    if (aa=='Y') return "TYR";
    if (aa=='Z') return "GLX";
    if ('a'<=aa && aa<='z')
    {
        char aaa[3]={' ',' ',(char)(toupper(aa))};
        return (string)(aaa);
    }
    return "UNK";
}

/* ShowSeqLen - whether to show residue number for each chain */
string pdb2fasta(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=0)
{
    stringstream buf;
    string sequence="";
    for (int c=0;c<pep.chains.size();c++)
    {
        sequence=pdb2fasta(pep.chains[c]);
        buf<<'>'<<PDBid<<':'<<pep.chains[c].chainID;
        if (ShowSeqLen) buf<<'\t'<<sequence.length();
        buf<<'\n'<<sequence<<'\n';
    }
    sequence.clear();
    return buf.str();
}

/* count the number of atoms with specific name in a residue */
int has_atom_name(ResidueUnit residue,string name=" CA ")
{
    int atom_name_count=0;
    for (int a=0;a<residue.atoms.size();a++)
        if (residue.atoms[a].name==name) atom_name_count++;
    return atom_name_count;
}

/* remove sidechain or backbone atoms 
 * atomic_detail - 1: only remove sidechain atoms
 *                 0: remove all non-CA atom*/
void remove_sidechain(ResidueUnit& residue,int atomic_detail=1)
{
    vector<AtomUnit> atoms; // list of atoms
    for (int a=0;a<residue.atoms.size();a++)
    {
        if ((atomic_detail==0 && residue.atoms[a].name==" CA ")||
            (atomic_detail==1 &&(residue.atoms[a].name==" CA " ||
             residue.atoms[a].name==" N  " || residue.atoms[a].name==" C  "
          || residue.atoms[a].name==" O  ")))
            atoms.push_back(residue.atoms[a]);
    }
    residue.atoms=atoms;
    atoms.clear();
}

void remove_sidechain(ChainUnit& chain,int atomic_detail=1)
{
    for (int r=0;r<chain.residues.size();r++)
        remove_sidechain(chain.residues[r],atomic_detail);
}

void remove_sidechain(ModelUnit& pep,int atomic_detail=1)
{
    for (int c=0;c<pep.chains.size();c++)
        remove_sidechain(pep.chains[c],atomic_detail);
}

void initialize_atom_order_map(map<string, map<string,int> > & ordMap)
{
    ordMap.clear();
    map<string,int> ord;
    ord[" P  "]= 0; ord[" OP1"]= 1; ord[" OP2"]= 2; ord[" O5'"]= 3;
    ord[" C5'"]= 4; ord[" C4'"]= 5; ord[" O4'"]= 6; ord[" C3'"]= 7;
    ord[" O3'"]= 8; ord[" C2'"]= 9; ord[" O2'"]=10; ord[" C1'"]=11;
    ordMap["  A"]=ord;
    ordMap["  C"]=ord;
    ordMap["  G"]=ord;
    ordMap["  U"]=ord;
    ord.clear();


    ordMap["  A"][" N9 "]=12; ordMap["  A"][" C8 "]=13; ordMap["  A"][" N7 "]=14;
    ordMap["  A"][" C5 "]=15; ordMap["  A"][" C6 "]=16; ordMap["  A"][" N6 "]=17;
    ordMap["  A"][" N1 "]=18; ordMap["  A"][" C2 "]=19; ordMap["  A"][" N3 "]=20;
    ordMap["  A"][" C4 "]=21;

    ordMap["  C"][" N1 "]=12; ordMap["  C"][" C2 "]=13; ordMap["  C"][" O2 "]=14;
    ordMap["  C"][" N3 "]=15; ordMap["  C"][" C4 "]=16; ordMap["  C"][" N4 "]=17;
    ordMap["  C"][" C5 "]=18; ordMap["  C"][" C6 "]=19;

    ordMap["  G"][" N9 "]=12; ordMap["  G"][" C8 "]=13; ordMap["  G"][" N7 "]=14;
    ordMap["  G"][" C5 "]=15; ordMap["  G"][" C6 "]=16; ordMap["  G"][" O6 "]=17;
    ordMap["  G"][" N1 "]=18; ordMap["  G"][" C2 "]=19; ordMap["  G"][" N2 "]=20;
    ordMap["  G"][" N3 "]=21; ordMap["  G"][" C4 "]=22;

    ordMap["  U"][" N1 "]=12; ordMap["  U"][" C2 "]=13; ordMap["  U"][" O2 "]=14;
    ordMap["  U"][" N3 "]=15; ordMap["  U"][" C4 "]=16; ordMap["  U"][" O4 "]=17;
    ordMap["  U"][" C5 "]=18; ordMap["  U"][" C6 "]=19;

    ord[" N  "]=0;
    ord[" CA "]=1;
    ord[" C  "]=2;
    ord[" O  "]=3;
    ordMap["GLY"]=ord;
    ord[" CB "]=3;
    ord[" O  "]=4;
    ordMap["ALA"]=ord;
    ordMap["CYS"]=ord; ordMap["CYS"][" SG "]=5;
    ordMap["ASP"]=ord; ordMap["ASP"][" CG "]=5; ordMap["ASP"][" OD1"]=6;
                       ordMap["ASP"][" OD2"]=7;
    ordMap["GLU"]=ord; ordMap["GLU"][" CG "]=5; ordMap["GLU"][" CD "]=6;
                       ordMap["GLU"][" OE1"]=7; ordMap["GLU"][" OE2"]=8;
    ordMap["PHE"]=ord; ordMap["PHE"][" CG "]=5; ordMap["PHE"][" CD1"]=6;
                       ordMap["PHE"][" CD2"]=7; ordMap["PHE"][" CE1"]=8;
                       ordMap["PHE"][" CE2"]=9; ordMap["PHE"][" CZ "]=10;
    ordMap["HIS"]=ord; ordMap["HIS"][" CG "]=5; ordMap["HIS"][" CD2"]=6;
                       ordMap["HIS"][" ND1"]=7; ordMap["HIS"][" CE1"]=8;
                       ordMap["HIS"][" NE2"]=9;
    ordMap["ILE"]=ord; ordMap["ILE"][" CG1"]=5; ordMap["ILE"][" CG2"]=6;
                       ordMap["ILE"][" CD1"]=7;
    ordMap["LYS"]=ord; ordMap["LYS"][" CG "]=5; ordMap["LYS"][" CD "]=6;
                       ordMap["LYS"][" CE "]=7; ordMap["LYS"][" NZ "]=8;
    ordMap["LEU"]=ord; ordMap["LEU"][" CG "]=5; ordMap["LEU"][" CD1"]=6;
                       ordMap["LEU"][" CD2"]=7;
    ordMap["MET"]=ord; ordMap["MET"][" CG "]=5; ordMap["MET"][" SD "]=6;
                       ordMap["MET"][" CE "]=7;
    ordMap["ASN"]=ord; ordMap["ASN"][" CG "]=5; ordMap["ASN"][" ND2"]=6;
                       ordMap["ASN"][" OD1"]=7;
    ordMap["PRO"]=ord; ordMap["PRO"][" CG "]=5; ordMap["PRO"][" CD "]=6;
    ordMap["GLN"]=ord; ordMap["GLN"][" CG "]=5; ordMap["GLN"][" CD "]=6;
                       ordMap["GLN"][" NE2"]=7; ordMap["GLN"][" OE1"]=8;
    ordMap["ARG"]=ord; ordMap["ARG"][" CG "]=5; ordMap["ARG"][" CD "]=6;
                       ordMap["ARG"][" NE "]=7; ordMap["ARG"][" NH1"]=8;
                       ordMap["ARG"][" NH2"]=9; ordMap["ARG"][" CZ "]=10;
    ordMap["SER"]=ord; ordMap["SER"][" OG "]=5;
    ordMap["THR"]=ord; ordMap["THR"][" CG2"]=5; ordMap["THR"][" OG1"]=6;
    ordMap["VAL"]=ord; ordMap["VAL"][" CG1"]=5; ordMap["VAL"][" CG2"]=6;
    ordMap["TRP"]=ord; ordMap["TRP"][" CG "]=5; ordMap["TRP"][" CD1"]=6;
                       ordMap["TRP"][" CD2"]=7; ordMap["TRP"][" CE2"]=8;
                       ordMap["TRP"][" CE3"]=9; ordMap["TRP"][" NE1"]=10;
                       ordMap["TRP"][" CH2"]=11;ordMap["TRP"][" CZ2"]=12;
                       ordMap["TRP"][" CZ3"]=13;
    ordMap["TYR"]=ord; ordMap["TYR"][" CG "]=5; ordMap["TYR"][" CD1"]=6;
                       ordMap["TYR"][" CD2"]=7; ordMap["TYR"][" CE1"]=8;
                       ordMap["TYR"][" CE2"]=9; ordMap["TYR"][" OH "]=10;
                       ordMap["TYR"][" CZ "]=11;
}

void standardize_pdb_order(ChainUnit &chain, map<string, map<string,int> >&ordMap)
{
    int r,a,order;
    AtomUnit temp_atom; 
    temp_atom.xyz.assign(3,0);
    ResidueUnit temp_residue;
    temp_residue.atoms.assign(23,temp_atom); // Longest residue length is 23 atoms, and each atom has x,y,z coordinates
    string resn,name;
    for (r=0;r<chain.residues.size();r++)
    {
        resn=chain.residues[r].resn;
        if (ordMap.count(resn)==0)
        {
            cerr<<"ERROR! unknown residue "<<resn<<endl;
            exit(1);
        }
        if (chain.residues[r].atoms.size()!=ordMap[resn].size() &&
            !(r==chain.residues.size()-1 &&
              chain.residues[r].atoms.size()==1+ordMap[resn].size()))
        {
            cerr<<"ERROR! "<<resn
                <<' '<<chain.chainID
                <<' '<<chain.residues[r].resi
                <<" has "<<chain.residues[r].atoms.size()
                <<" atoms != "<<ordMap[resn].size()<<endl;
            exit(1);
        }

        // Loop through atoms in each residue
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            name=chain.residues[r].atoms[a].name;
            if (name==" OXT") order=ordMap[resn].size();
            else if (ordMap[resn].count(name)==0)
            {
                cerr<<"Skip atom "<<resn<<' '<<name<<endl;
                continue;
            }
            else order=ordMap[resn][name];

            temp_residue.atoms[order].xyz[0] = chain.residues[r].atoms[a].xyz[0];
            temp_residue.atoms[order].xyz[1] = chain.residues[r].atoms[a].xyz[1];
            temp_residue.atoms[order].xyz[2] = chain.residues[r].atoms[a].xyz[2];
            temp_residue.atoms[order].name   = name;
            temp_residue.atoms[order].bfactor= chain.residues[r].atoms[a].bfactor;
        }
        
        // Use temp_residue to update the actual peptide with the correct order
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            chain.residues[r].atoms[a].xyz[0] = temp_residue.atoms[a].xyz[0];
            chain.residues[r].atoms[a].xyz[1] = temp_residue.atoms[a].xyz[1];
            chain.residues[r].atoms[a].xyz[2] = temp_residue.atoms[a].xyz[2];
            chain.residues[r].atoms[a].name   = temp_residue.atoms[a].name; 
            chain.residues[r].atoms[a].bfactor= temp_residue.atoms[a].bfactor; 
        }
    
    }
}

void standardize_pdb_order(ModelUnit &pep, map<string, map<string,int> >&ordMap)
{
    int c;
    for (c=0;c<pep.chains.size();c++)
        standardize_pdb_order(pep.chains[c],ordMap);
}

char check_moltype(ChainUnit &chain)
{
    int pro=0;
    int dna=0;
    int rna=0;
    int r;
    for (r=0;r<chain.residues.size();r++)
    {
        if (chain.residues[r].resn[0]!=' ')      pro++;
        else if (chain.residues[r].resn[1]=='D') dna++;
        else if (chain.residues[r].resn[1]==' ') rna++;
    }
    if      (dna>pro && dna>rna) return 'D';
    else if (rna>pro && rna>dna) return 'R';
    return 'P';
}

int check_bfactor_mode(ChainUnit &chain)
{
    int bfactor_mode=-1;
    int32_t bfactor=0;
    int r,a;
    for (r=0;r<chain.residues.size();r++)
    {
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (bfactor_mode==-1)
            {
                bfactor=chain.residues[r].atoms[a].bfactor;
                bfactor_mode=0;
            }
            if (bfactor_mode==0 && bfactor!=
                chain.residues[r].atoms[a].bfactor)
            {
                if (a) return 2;
                else bfactor_mode=1;
            }
            if (a==0) bfactor=chain.residues[r].atoms[a].bfactor;
            if (bfactor!=chain.residues[r].atoms[a].bfactor)
                return 2;
        }
    }
    return bfactor_mode;
}

string write_pdc_structure(ModelUnit &pep,string &header)
{
    string txt=header;
    stringstream buf;
    int a,c,r;
    char moltype;
    int bfactor_mode;
    pdb2fasta(pep);
    int resi_prev=0;
    char icode_prev=' ';
    int L;
    int32_t x,y,z,bfactor;
    int32_t prev_x=pep.chains[c].residues[0].atoms[0].xyz[0];
    int32_t prev_y=pep.chains[c].residues[0].atoms[0].xyz[1];
    int32_t prev_z=pep.chains[c].residues[0].atoms[0].xyz[2];
    int32_t dxf,dyf,dzf,dbf;
    int16_t dx16,dy16,dz16,db16;
    for (c=0;c<pep.chains.size();c++)
    {
        moltype=check_moltype(pep.chains[c]);
        bfactor_mode=check_bfactor_mode(pep.chains[c]);
        L=pep.chains[c].residues.size();
        buf <<"@a\n>"<<pep.chains[c].chainID<<'\t'
            <<moltype<<'\t'<<bfactor_mode<<'\t'<<L<<'\n'
            <<pep.chains[c].sequence<<'\n';
        for (r=0;r<L;r++)
        {
            if (r==0)
            {
                resi_prev =pep.chains[c].residues[0].resi;
                icode_prev=pep.chains[c].residues[0].icode;
                buf<<resi_prev;
                if (icode_prev!=' ') buf<<icode_prev;
            }
            else
            {
                if (icode_prev!=' ')
                {
                    buf<<','<<pep.chains[c].residues[r].resi;
                    if (pep.chains[c].residues[r].icode!=' ')
                        buf<<pep.chains[c].residues[r].icode;
                }
                else
                {
                    if (pep.chains[c].residues[r].icode!=' ')
                        buf<<'~'<<resi_prev
                           <<','<<pep.chains[c].residues[r].resi
                           <<pep.chains[c].residues[r].icode;
                    else if (pep.chains[c].residues[r].resi!=resi_prev+1)
                        buf<<'~'<<resi_prev;
                    else if (r==L-1)
                        buf<<'~'<<pep.chains[c].residues[r].resi;
                }
                resi_prev =pep.chains[c].residues[r].resi;
                icode_prev=pep.chains[c].residues[r].icode;
            }
        }
        buf<<'\n';
        buf.write((char *)&prev_x,sizeof(int32_t));
        buf.write((char *)&prev_y,sizeof(int32_t));
        buf.write((char *)&prev_z,sizeof(int32_t));
        for (r=0;r<pep.chains[c].residues.size();r++)
        {
            for (a=0;a<pep.chains[c].residues[r].atoms.size();a++)
            {
                x=pep.chains[c].residues[r].atoms[a].xyz[0];
                y=pep.chains[c].residues[r].atoms[a].xyz[1];
                z=pep.chains[c].residues[r].atoms[a].xyz[2];
                dx16=x-prev_x;
                dy16=y-prev_y;
                dz16=z-prev_z;
                buf.write((char *)&dx16,sizeof(int16_t));
                buf.write((char *)&dy16,sizeof(int16_t));
                buf.write((char *)&dz16,sizeof(int16_t));
                prev_x=x;
                prev_y=y;
                prev_z=z;
            }
            prev_x=pep.chains[c].residues[r].atoms[2].xyz[0];
            prev_y=pep.chains[c].residues[r].atoms[2].xyz[1];
            prev_z=pep.chains[c].residues[r].atoms[2].xyz[2];
        }
        if (bfactor_mode==0)
        {
            bfactor=pep.chains[c].residues[0].atoms[0].bfactor;
            buf.write((char *)&bfactor,sizeof(int16_t));
        }
        else if (bfactor_mode==1)
        {
            for (r=0;r<pep.chains[c].residues.size();r++)
            {
                bfactor=pep.chains[c].residues[r].atoms[0].bfactor;
                buf.write((char *)&bfactor,sizeof(int16_t));
            }
        }
        else if (bfactor_mode==2)
        {
            for (r=0;r<pep.chains[c].residues.size();r++)
            {
                for (a=0;a<pep.chains[c].residues[r].atoms.size();a++)
                {
                    bfactor=pep.chains[c].residues[r].atoms[a].bfactor;
                    buf.write((char *)&bfactor,sizeof(int16_t));
                }
            }
        }
    }
    buf<<flush;
    txt+=buf.str()+"\nEND\n";
    buf.str(string());
    return txt;
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_pdc_structure(const char *filename,ModelUnit &pep,string &header)
{
    if (strcmp(filename,"-")==0)
        cout<<write_pdc_structure(pep,header)<<flush;
    else
    {
        ofstream fp(filename);
        fp<<write_pdc_structure(pep,header)<<flush;
        fp.close();
    }
}

void initialize_reverse_atom_order_map(map<string, vector<string> > & ordMapR)
{
    map<string, map<string,int> >ordMap;
    initialize_atom_order_map(ordMap);
    vector<string> atomName_vec;
    for ( const auto &myPair : ordMap )
    {
        atomName_vec.assign(myPair.second.size(),"");
        for ( const auto &myPair2 : myPair.second )
            atomName_vec[myPair2.second]=myPair2.first;
        ordMapR[myPair.first]=atomName_vec;
        atomName_vec.clear();
    }
}

ModelUnit read_pdc_structure(const char *filename,string &header,
    map<string, vector<string> > & ordMapR)
{
    ModelUnit pep;
    ChainUnit chain;
    ResidueUnit residue;
    AtomUnit atom;
    atom.xyz.assign(3,0);
    atom.bfactor=0;
    string filename_str=(string) filename;
    string line;
    
    int use_stdin=(filename_str=="-");
    int use_pstream=0; // input is compressed

    ifstream fp;
    redi::ipstream fp_gz; // if file is compressed
    if (filename_str.length()>=3 && 
        filename_str.substr(filename_str.length()-3,3)==".gz")
    {
        // gzip pdb
        fp_gz.open("zcat "+filename_str);
        use_pstream=1;
    }
    else
    {
        fp.open(filename,ios::in); //ifstream fp(filename,ios::in);
    }

    string sequence;
    int L;
    vector<string> line_vec;
    int bfactor_mode;
    char moltype;
    int a,r,c;
    int32_t prev_x,prev_y,prev_z;
    int32_t x,y,z,bfactor;
    int16_t dx16,dy16,dz16,db16;
    int atomNum;
    while(use_stdin?cin.good():(use_pstream?fp_gz.good():fp.good()))
    {
        if (use_stdin)        getline(cin,line);
        else if (use_pstream) getline(fp_gz,line);
        else                  getline(fp,line);

        if (line.substr(0,3)=="END") break;
        if (line.substr(0,2)!="@a")
        {
            header+=line+'\n';
            continue;
        }
        
        if (use_stdin)        getline(cin,line);
        else if (use_pstream) getline(fp_gz,line);
        else                  getline(fp,line);
        Split(line,line_vec,'\t');
        chain.chainID=line_vec[0][1];
        moltype=line_vec[1][0];
        bfactor_mode=atoi(line_vec[2].c_str());
        L=atoi(line_vec[3].c_str());
        line_vec.clear();
        c=pep.chains.size();
        if (use_stdin)        getline(cin,sequence);
        else if (use_pstream) getline(fp_gz,sequence);
        else                  getline(fp,sequence);
        if (use_stdin)        getline(cin,line);
        else if (use_pstream) getline(fp_gz,line);
        else                  getline(fp,line);
        Split(line,line_vec,',');
        vector<int>  resi_vec;
        vector<char> icode_vec;
        int rmin,rmax;
        for (a=0;a<line_vec.size();a++)
        {
            line=line_vec[a];
            vector<string> resi_str_vec;
            Split(line,resi_str_vec,'~');
            if (resi_str_vec.size()!=2)
            {
                if (line.back()<'0' || line.back()>'9')
                {
                    icode_vec.push_back(line.back());
                    line=line.substr(0,line.size()-1);
                }
                else icode_vec.push_back(' ');
                resi_vec.push_back(atoi(line.c_str()));
            }
            else
            {
                rmin=atoi(resi_str_vec[0].c_str());
                rmax=atoi(resi_str_vec[1].c_str());
                for (r=rmin;r<=rmax;r++)
                {
                    resi_vec.push_back(r);
                    icode_vec.push_back(' ');
                }
            }
        }
        line_vec.clear();
        if (use_stdin)
        {
            cin.read((char *)&prev_x,sizeof(int32_t));
            cin.read((char *)&prev_y,sizeof(int32_t));
            cin.read((char *)&prev_z,sizeof(int32_t));
        }
        else if (use_pstream)
        {
            fp_gz.read((char *)&prev_x,sizeof(int32_t));
            fp_gz.read((char *)&prev_y,sizeof(int32_t));
            fp_gz.read((char *)&prev_z,sizeof(int32_t));
        }
        else
        {
            fp.read((char *)&prev_x,sizeof(int32_t));
            fp.read((char *)&prev_y,sizeof(int32_t));
            fp.read((char *)&prev_z,sizeof(int32_t));
        }
        for (r=0;r<L;r++)
        {
            residue.resi=resi_vec[r];
            residue.icode=icode_vec[r];
            residue.resn=aa1to3(sequence[r]);
            for (a=0;a<ordMapR[residue.resn].size()+(r+1==L);a++)
            {
                if (use_stdin)
                {
                    cin.read((char *)&dx16,sizeof(int16_t));
                    cin.read((char *)&dy16,sizeof(int16_t));
                    cin.read((char *)&dz16,sizeof(int16_t));
                }
                else if (use_pstream)
                {
                    fp_gz.read((char *)&dx16,sizeof(int16_t));
                    fp_gz.read((char *)&dy16,sizeof(int16_t));
                    fp_gz.read((char *)&dz16,sizeof(int16_t));
                }
                else
                {
                    fp.read((char *)&dx16,sizeof(int16_t));
                    fp.read((char *)&dy16,sizeof(int16_t));
                    fp.read((char *)&dz16,sizeof(int16_t));
                }
                prev_x+=dx16;
                prev_y+=dy16;
                prev_z+=dz16;
                atom.xyz[0]=prev_x;
                atom.xyz[1]=prev_y;
                atom.xyz[2]=prev_z;
                atom.name=" OXT";
                if (a<ordMapR[residue.resn].size()) atom.name=ordMapR[residue.resn][a];
                residue.atoms.push_back(atom);
            }
            prev_x=residue.atoms[2].xyz[0];
            prev_y=residue.atoms[2].xyz[1];
            prev_z=residue.atoms[2].xyz[2];
            chain.residues.push_back(residue);
            residue.atoms.clear();
        }
        if (bfactor_mode==0)
        {
            if (use_stdin)        cin.read((char *)&bfactor,sizeof(int16_t));
            else if (use_pstream) fp_gz.read((char *)&bfactor,sizeof(int16_t));
            else                  fp.read((char *)&bfactor,sizeof(int16_t));
            for (r=0;r<chain.residues.size();r++)
                for (a=0;a<chain.residues[r].atoms.size();a++)
                    chain.residues[r].atoms[a].bfactor=bfactor;
        }
        else if (bfactor_mode==1)
        {
            for (r=0;r<chain.residues.size();r++)
            {
                if (use_stdin)        cin.read((char *)&bfactor,sizeof(int16_t));
                else if (use_pstream) fp_gz.read((char *)&bfactor,sizeof(int16_t));
                else                  fp.read((char *)&bfactor,sizeof(int16_t));
                for (a=0;a<chain.residues[r].atoms.size();a++)
                    chain.residues[r].atoms[a].bfactor=bfactor;
            }
        }
        else if (bfactor_mode==2)
        {
            for (r=0;r<chain.residues.size();r++)
            {
                for (a=0;a<chain.residues[r].atoms.size();a++)
                {
                    if (use_stdin)        cin.read((char *)&bfactor,sizeof(int16_t));
                    else if (use_pstream) fp_gz.read((char *)&bfactor,sizeof(int16_t));
                    else                  fp.read((char *)&bfactor,sizeof(int16_t));
                    chain.residues[r].atoms[a].bfactor=bfactor;
                }
            }
        }
        pep.chains.push_back(chain);
        chain.residues.clear();
    }
    if (!use_stdin)
    {
        if (use_pstream==0) fp.close();
        else             fp_gz.close();
    }
    atom.xyz.clear();
    residue.atoms.clear();
    chain.residues.clear();
    vector<string>().swap(line_vec);
    return pep;
}

string write_cif_structure(ModelUnit &pep,string &header)
{
    string txt="data_AF\n#\n_entry.id AF\n#\n";
    pdb2fasta(pep);
    stringstream buf;
    int a,r,c,s;
    vector<string> line_vec;
    Split(header,line_vec,'\n');
    for (s=0;s<line_vec.size();s++)
    {
        if (StartsWith(line_vec[s],"DBREF "))
        {
            txt=Trim(line_vec[s].substr(33,8));
            r=atoi(line_vec[s].substr(55,5).c_str());
            r/=200;
            r++;
            buf<<"data_AF-"<<txt<<"-F"<<r
                <<"\n#\n_entry.id AF-"<<txt<<"-F"<<r<<"\n#\n";
            txt=buf.str();
            buf.str(string());
            break;
        }
    }
    vector<string> ().swap(line_vec);


    txt+=""
"loop_\n"
"_atom_type.symbol\n"
"C \n"
"N \n"
"O \n";
    bool hasS=false;
    bool hasP=false;
    string sequence="";
    for (c=0;c<pep.chains.size();c++)
    {
        sequence+=pep.chains[c].sequence;
        char moltype=check_moltype(pep.chains[c]);
        if (moltype=='D'||moltype=='R') hasP=true;
        else if (hasS==false &&
            pep.chains[c].sequence.find('M')!=string::npos ||
            pep.chains[c].sequence.find('C')!=string::npos)
            hasS=true;
    }
    if (hasP) txt+="P \n";
    if (hasS) txt+="S \n";
    txt+="#\n";
    txt+="loop_\n"
"_entity_poly_seq.entity_id\n"
"_entity_poly_seq.hetero\n"
"_entity_poly_seq.mon_id\n"
"_entity_poly_seq.num\n";

    double global_metric=0;
    int L=sequence.size();
    for (c=0;c<pep.chains.size();c++)
    {
        for (r=0;r<pep.chains[c].residues.size();r++)
        {
            buf<<"1 n "<<pep.chains[c].residues[r].resn<<' '
                <<left<<setw(5)<<pep.chains[c].residues[r].resi<<endl;
            txt+=buf.str();
            buf.str(string());
            global_metric+=pep.chains[c].residues[r].atoms[1].bfactor;
        }
    }
    global_metric/=(100*L);
    txt+="#\n"
"loop_\n"
"_ma_data.content_type\n"
"_ma_data.id\n"
"_ma_data.name\n"
"\"model coordinates\" 1 Model             \n"
"\"input structure\"   2 \"Input structure\" \n"
"#\n"
"loop_\n"
"_ma_protocol_step.method_type\n"
"_ma_protocol_step.ordinal_id\n"
"_ma_protocol_step.protocol_id\n"
"_ma_protocol_step.step_id\n"
"\"coevolution MSA\" 1 1 1 \n"
"\"template search\" 2 1 2 \n"
"modeling          3 1 3 \n"
"#\n"
"_ma_software_group.group_id    1\n"
"_ma_software_group.ordinal_id  1\n"
"_ma_software_group.software_id 1\n"
"#\n"
"_ma_template_trans_matrix.id               1\n"
"_ma_template_trans_matrix.rot_matrix[1][1] 1.0\n"
"_ma_template_trans_matrix.rot_matrix[1][2] 0.0\n"
"_ma_template_trans_matrix.rot_matrix[1][3] 0.0\n"
"_ma_template_trans_matrix.rot_matrix[2][1] 0.0\n"
"_ma_template_trans_matrix.rot_matrix[2][2] 1.0\n"
"_ma_template_trans_matrix.rot_matrix[2][3] 0.0\n"
"_ma_template_trans_matrix.rot_matrix[3][1] 0.0\n"
"_ma_template_trans_matrix.rot_matrix[3][2] 0.0\n"
"_ma_template_trans_matrix.rot_matrix[3][3] 1.0\n"
"_ma_template_trans_matrix.tr_vector[1]     0.0\n"
"_ma_template_trans_matrix.tr_vector[2]     0.0\n"
"_ma_template_trans_matrix.tr_vector[3]     0.0\n"
"#\n"
"loop_\n"
"_pdbx_poly_seq_scheme.asym_id\n"
"_pdbx_poly_seq_scheme.auth_seq_num\n"
"_pdbx_poly_seq_scheme.entity_id\n"
"_pdbx_poly_seq_scheme.hetero\n"
"_pdbx_poly_seq_scheme.mon_id\n"
"_pdbx_poly_seq_scheme.pdb_ins_code\n"
"_pdbx_poly_seq_scheme.pdb_mon_id\n"
"_pdbx_poly_seq_scheme.pdb_seq_num\n"
"_pdbx_poly_seq_scheme.pdb_strand_id\n"
"_pdbx_poly_seq_scheme.seq_id\n";
    for (c=0;c<pep.chains.size();c++)
    {
        for (r=0;r<pep.chains[c].residues.size();r++)
        {
            buf<<pep.chains[c].chainID<<' '
                <<left<<setw(4)<<pep.chains[c].residues[r].resi<<" 1 n "
                <<pep.chains[c].residues[r].resn<<" . "
                <<pep.chains[c].residues[r].resn<<' '
                <<left<<setw(4)<<pep.chains[c].residues[r].resi
                <<" "<<pep.chains[c].chainID<<" "
                <<left<<setw(5)<<pep.chains[c].residues[r].resi<<endl;
            txt+=buf.str();
            buf.str(string());
        }
    }
    txt+="#\n"
"loop_\n"
"_atom_site.group_PDB\n"
"_atom_site.id\n"
"_atom_site.type_symbol\n"
"_atom_site.label_atom_id\n"
"_atom_site.label_alt_id\n"
"_atom_site.label_comp_id\n"
"_atom_site.label_asym_id\n"
"_atom_site.label_entity_id\n"
"_atom_site.label_seq_id\n"
"_atom_site.pdbx_PDB_ins_code\n"
"_atom_site.Cartn_x\n"
"_atom_site.Cartn_y\n"
"_atom_site.Cartn_z\n"
"_atom_site.occupancy\n"
"_atom_site.B_iso_or_equiv\n"
"_atom_site.pdbx_formal_charge\n"
"_atom_site.auth_seq_id\n"
"_atom_site.auth_comp_id\n"
"_atom_site.auth_asym_id\n"
"_atom_site.auth_atom_id\n"
"_atom_site.pdbx_PDB_model_num\n";
    size_t i=0;
    char icode='?';
    for (c=0;c<pep.chains.size();c++)
    {
        for (r=0;r<pep.chains[c].residues.size();r++)
        {
            icode=pep.chains[c].residues[r].icode;
            if (icode==' ') icode='?';
            for (a=0;a<pep.chains[c].residues[r].atoms.size();a++)
            {
                i++;

                buf<<"ATOM "
                    <<left<<setw(5)<<i<<" "
                    <<Trim(pep.chains[c].residues[r].atoms[a].name)[0]
                    <<' '<<left<<setw(3)
                    <<Trim(pep.chains[c].residues[r].atoms[a].name)<<" . "
                    <<pep.chains[c].residues[r].resn
                    <<' '<<pep.chains[c].chainID<<" 1 "
                    <<left<<setw(4)<<pep.chains[c].residues[r].resi
                    <<' '<<icode<<' '
                    <<left<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)
                    <<0.001*pep.chains[c].residues[r].atoms[a].xyz[0]<<' '
                    <<left<<setw(7)<<setiosflags(ios::fixed)<<setprecision(3)
                    <<0.001*pep.chains[c].residues[r].atoms[a].xyz[1]<<' '
                    <<left<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)
                    <<0.001*pep.chains[c].residues[r].atoms[a].xyz[2]<<" 1.0 "
                    <<left<<setw(5)<<setiosflags(ios::fixed)<<setprecision(2)
                    <<0.01*pep.chains[c].residues[r].atoms[a].bfactor
                    <<' '<<icode<<' '
                    <<left<<setw(4)<<pep.chains[c].residues[r].resi
                    <<' '<<pep.chains[c].residues[r].resn
                    <<' '<<pep.chains[c].chainID
                    <<' '<<left<<setw(3)
                    <<Trim(pep.chains[c].residues[r].atoms[a].name)<<" 1 "
                    <<endl;
                txt+=buf.str();
                buf.str(string());
            }
        }
    }
    txt+="#\n";
    return txt;
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_cif_structure(const char *filename,ModelUnit &pep,string &header)
{
    if (strcmp(filename,"-")==0)
        cout<<write_cif_structure(pep,header)<<flush;
    else
    {
        ofstream fp(filename);
        fp<<write_cif_structure(pep,header)<<flush;
        fp.close();
    }
}

void deepClean(AtomUnit &atom)
{
    string ().swap(atom.name);
    vector<int32_t>().swap(atom.xyz);
}

void deepClean(ResidueUnit &residue)
{
    int a;
    for (a=0;a<residue.atoms.size();a++) deepClean(residue.atoms[a]);
    residue.atoms.clear();
    string ().swap(residue.resn);
}

void deepClean(ChainUnit &chain)
{
    int r;
    for (r=0;r<chain.residues.size();r++) deepClean(chain.residues[r]);
    chain.residues.clear();
}

void deepClean(ModelUnit &pep)
{
    int c;
    for (c=0;c<pep.chains.size();c++) deepClean(pep.chains[c]);
}

#endif
