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
#include "GeometryTools.hpp"
#include "Superpose.hpp"

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
    while (line.size()<8) line=' '+line;
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
    while (line.size()<6) line=' '+line;
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
    if (StartsWith(line,"TITLE ") || StartsWith(line,"COMPND") || 
        StartsWith(line,"SOURCE") || StartsWith(line,"DBREF")) return 2;
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
        fp_gz.open("gunzip -c "+filename_str);
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

ModelUnit read_cif_structure(const char *filename,string &header,
    const int atomic_detail=2,const int allowX=1)
{
    ModelUnit pep;
    ChainUnit chain;
    ResidueUnit residue;
    AtomUnit atom;
    atom.xyz.assign(3,0);
    atom.bfactor=0;
    int s,a,r,c;

    map<string,int> _atom_site;
    int _atom_site_count=0;
    stringstream buf;

    string line="";
    string filename_str=(string) filename;
    
    int use_stdin=(filename_str=="-");
    int use_pstream=0; // input is compressed
    ifstream fp;
    redi::ipstream fp_gz; // if file is compressed
    if (filename_str.length()>=3 && 
        filename_str.substr(filename_str.length()-3,3)==".gz")
    {
        // gzip pdb
        fp_gz.open("gunzip -c "+filename_str);
        use_pstream=1;
    }
    else
    {
        fp.open(filename,ios::in); //ifstream fp(filename,ios::in);
    }

    string model_group_name="ALPHAFOLD MONOMER V2.0 PREDICTION";
    string db_accession="";
    string db_code="UNP";
    string db_name="";
    string ncbi_taxonomy_id="";
    string organism_scientific="";
    string seq_db_align_begin="";
    string seq_db_align_end="";
    string pdbx_description="";
    string asym_id="";

    string model_num="";
    string resi="";
    string icode="";
    string chainID="";
    string seq_id="";
    string resn="";
    string group_PDB="";
    string atom_id="";
    vector<string> pdb_line_vec;
    
    while(use_stdin?cin.good():(use_pstream?fp_gz.good():fp.good()))
    {
        if (use_stdin)        getline(cin,line);
        else if (use_pstream) getline(fp_gz,line);
        else                  getline(fp,line);

        if (StartsWith(line,"_ma_model_list.model_group_name"))
            model_group_name=Trim(line.substr(32)," \"");
        else if (StartsWith(line,"_ma_target_ref_db_details.db_accession"))
            db_accession=Trim(line.substr(39));
        else if (StartsWith(line,"_ma_target_ref_db_details.db_code"))
            db_code=Trim(line.substr(34));
        else if (StartsWith(line,"_ma_target_ref_db_details.db_name"))
            db_name=Trim(line.substr(34));
        else if (StartsWith(line,"_ma_target_ref_db_details.ncbi_taxonomy_id"))
            ncbi_taxonomy_id=Trim(line.substr(43));
        else if (StartsWith(line,"_ma_target_ref_db_details.organism_scientific"))
            organism_scientific=Trim(line.substr(46)," \"");
        else if (StartsWith(line,"_ma_target_ref_db_details.seq_db_align_begin"))
            seq_db_align_begin=Trim(line.substr(45));
        else if (StartsWith(line,"_ma_target_ref_db_details.seq_db_align_end"))
            seq_db_align_end=Trim(line.substr(43));
        else if (StartsWith(line,"_entity.pdbx_description"))
            pdbx_description=Trim(line.substr(25)," \"");
        else if (StartsWith(line,"_ma_target_entity_instance.asym_id"))
            asym_id=Trim(line.substr(35));
        else if (StartsWith(line,"loop_") || StartsWith(line,"#"))
        {
            _atom_site_count=0;
            _atom_site["group_PDB"]=-1; // *
            _atom_site["atom_id"]=-1;   // * atom name
            _atom_site["alt_id"]=-1;
            _atom_site["comp_id"]=-1;   // * resn
            _atom_site["asym_id"]=-1;   // * chainID
            _atom_site["seq_id"]=-1;    // * resi
            _atom_site["icode"]=-1;     // insertion code
            _atom_site["model_num"]=-1; // model index
            _atom_site["Cartn_x"]=-1;   // *
            _atom_site["Cartn_y"]=-1;   // *
            _atom_site["Cartn_z"]=-1;   // *
            _atom_site["bfactor"]=-1;
            model_num=-1;
        }
        else if (StartsWith(line,"_atom_site."))
        {
            if (StartsWith(line,"_atom_site.group_PDB"))
                _atom_site["group_PDB"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.label_atom_id"))
                _atom_site["atom_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.label_alt_id"))
                _atom_site["alt_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.label_comp_id") && 
                _atom_site["comp_id"]<0)
                _atom_site["comp_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.auth_comp_id"))
                _atom_site["comp_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.label_asym_id") &&
                _atom_site["asym_id"]<0)
                _atom_site["asym_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.auth_asym_id"))
                _atom_site["asym_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.label_seq_id") &&
                 _atom_site["seq_id"]<0)
                _atom_site["seq_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.auth_seq_id"))
                _atom_site["seq_id"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.pdbx_PDB_ins_code"))
                _atom_site["icode"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.pdbx_PDB_model_num"))
                _atom_site["model_num"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.Cartn_x"))
                _atom_site["Cartn_x"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.Cartn_y"))
                _atom_site["Cartn_y"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.Cartn_z"))
                _atom_site["Cartn_z"]=_atom_site_count;
            else if (StartsWith(line,"_atom_site.B_iso_or_equiv"))
                _atom_site["bfactor"]=_atom_site_count;
            _atom_site_count++;
        }
        else if (_atom_site_count>=8 && !StartsWith(line,"_atom_site."))
        {
            if (_atom_site["group_PDB"]<0 || _atom_site["atom_id"]<0 ||
                _atom_site["comp_id"]<0   || _atom_site["asym_id"]<0 ||
                _atom_site["seq_id"]<0    || _atom_site["Cartn_x"]<0 ||
                _atom_site["Cartn_y"]<0   || _atom_site["Cartn_z"]<0)
                continue;
            Split(line,pdb_line_vec);
            group_PDB=pdb_line_vec[_atom_site["group_PDB"]];
            resn=pdb_line_vec[_atom_site["comp_id"]];
            if ((allowX==0 && group_PDB!="ATOM")||
                (allowX==1 && group_PDB!="ATOM" &&  
               !(group_PDB=="HETATM" && resn=="MSE"))||
                (group_PDB!="ATOM" && group_PDB!="HETATM"))
                continue;
            if (_atom_site["model_num"]>=0)
                if (model_num.size()>=0) model_num=pdb_line_vec[_atom_site["model_num"]];
                else if (model_num!=pdb_line_vec[_atom_site["model_num"]]) continue;
            if (_atom_site["alt_id"]>=0 && pdb_line_vec[_atom_site["alt_id"]]!="."
                && pdb_line_vec[_atom_site["alt_id"]]!="A") continue;
            atom.name=pdb_line_vec[_atom_site["atom_id"]];
            if (resn=="MSE" && atom.name=="SE") atom.name=" SD ";
            if      (resn.size()==1) resn="  "+resn;
            else if (resn.size()==2) resn=" "+resn;
            else if (resn=="MSE") resn="MET";
            if      (atom.name.size()==1) atom.name+=' ';
            if      (atom.name.size()==2) atom.name+=' ';
            if      (atom.name.size()==3) atom.name=' '+atom.name;
            if      (atom.name==" O1P")   atom.name=" OP1";
            else if (atom.name==" O2P")   atom.name=" OP2";
            else if (atom.name[3]=='*')   atom.name=atom.name.substr(0,3)+"'";
            atom.xyz[0]=XYZtoint32(pdb_line_vec[_atom_site["Cartn_x"]]);
            atom.xyz[1]=XYZtoint32(pdb_line_vec[_atom_site["Cartn_y"]]);
            atom.xyz[2]=XYZtoint32(pdb_line_vec[_atom_site["Cartn_z"]]);
            if (_atom_site["bfactor"]>=0) 
                atom.bfactor=Btoint16(pdb_line_vec[_atom_site["bfactor"]]);
            if (chainID!=pdb_line_vec[_atom_site["asym_id"]])
            {
                chain.residues.clear();
                chainID=pdb_line_vec[_atom_site["asym_id"]];
                if (chainID==".") chainID="_";
                chain.chainID=chainID[0];
                pep.chains.push_back(chain);
                residue.atoms.clear();
                residue.resn=resn;
                resi=pdb_line_vec[_atom_site["seq_id"]];
                residue.resi=atoi(resi.c_str());
                residue.icode=' ';
                if (_atom_site["icode"]>=0)
                {
                    icode=pdb_line_vec[_atom_site["icode"]];
                    residue.icode=icode[0];
                    if (icode=="." || icode=="?") residue.icode=' ';
                }
                pep.chains[c].residues.push_back(residue);
            }
            c=pep.chains.size()-1;
            r=pep.chains[c].residues.size()-1;
            if (resi!=pdb_line_vec[_atom_site["seq_id"]] || (
                _atom_site["icode"]>=0 && icode!=pdb_line_vec[_atom_site["icode"]]))
            {
                residue.atoms.clear();
                residue.resn=resn;
                resi=pdb_line_vec[_atom_site["seq_id"]];
                residue.resi=atoi(resi.c_str());
                residue.icode=' ';
                if (_atom_site["icode"]>=0)
                {
                    icode=pdb_line_vec[_atom_site["icode"]];
                    residue.icode=icode[0];
                    if (icode=="." || icode=="?") residue.icode=' ';
                }
                pep.chains[c].residues.push_back(residue);
                r=pep.chains[c].residues.size()-1;
            }
            pep.chains[c].residues[r].atoms.push_back(atom);
            for (s=0;s<pdb_line_vec.size();s++) pdb_line_vec[s].clear();
            pdb_line_vec.clear();
        }
    }
    vector<string> ().swap(pdb_line_vec);
    if (model_group_name.size())
        header+="TITLE     "+model_group_name+" FOR "+pdbx_description+
            " ("+db_accession+")\n";
    header+="COMPND    MOL_ID: 1;\n";
    if (pdbx_description.size()) header+="COMPND   2 MOLECULE: "+
        pdbx_description+";\n";
    if (asym_id.size())
        header+="COMPND   3 CHAIN: "+asym_id+"\n";
    header+="SOURCE    MOL_ID: 1;\n";
    if (organism_scientific.size())
        header+="SOURCE   2 ORGANISM_SCIENTIFIC: "+organism_scientific+";\n";
    if (ncbi_taxonomy_id.size())
        header+="SOURCE   3 ORGANISM_TAXID: "+ncbi_taxonomy_id+"\n";
    if (seq_db_align_begin.size() && seq_db_align_end.size() && 
        db_accession.size() && db_code.size())
    {
        buf<<"DBREF  XXXX "<<asym_id[0]<<"    1  "<<setw(4)
            <<1+atoi(seq_db_align_end.c_str())-atoi(seq_db_align_begin.c_str())
            <<"  "<<setw(6)<<left<<db_name
            <<" "<<setw(8)<<left<<db_accession
            <<" "<<setw(12)<<left<<db_code
            <<" "<<setw(5)<<right<<seq_db_align_begin
            <<"  "<<setw(5)<<seq_db_align_end<<endl;
        header+=buf.str();
        buf.str(string());
    }

    if (!use_stdin)
    {
        if (use_pstream==0) fp.close();
        else fp_gz.close();
    }
    
    /* clean up */
    map<string,int> ().swap(_atom_site);
    atom.xyz.clear();
    residue.atoms.clear();
    chain.residues.clear();
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
void write_pdb_structure(const string &filename,ModelUnit &pep,string &header)
{
    if (filename=="-")
        cout<<write_pdb_structure(pep,header)<<flush;
    else
    {
        string filename_str=filename;
        if (EndsWith(filename,".gz")) filename_str=filename.substr(0,filename.size()-3);
        ofstream fp(filename_str.c_str());
        fp<<write_pdb_structure(pep,header)<<flush;
        fp.close();
        if (EndsWith(filename,".gz"))
            int r=system(((string)("gzip -f "+filename_str)).c_str());
        filename_str.clear();
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

void get_bb_angle(const string &resn, double &N_CA_C_angle,
    double &CA_C_O_angle, double &N_C_CA_CB_diangle, double &N_CA_C_O_diangle)
{
    char aa=aa3to1(resn);

    N_CA_C_angle= 111.;
    CA_C_O_angle=120.5;
    N_C_CA_CB_diangle=122.6;
    N_CA_C_O_diangle =120.0; // LTRKDEQMHFYW

    if      (aa=='G')
    {
        N_CA_C_angle=110.8914;
        CA_C_O_angle=120.5117;
    }
    else if (aa=='A')
    {
        N_CA_C_angle=111.068;
    }
    else if (aa=='P')
    {
        N_CA_C_O_diangle=-45.0;
        CA_C_O_angle=120.2945;
        N_CA_C_angle=112.7499;
        N_C_CA_CB_diangle=115.2975;
    }
    else if (string("SCUVIN").find(aa)!=string::npos)
    {
        N_CA_C_O_diangle= -60.0;
        if (aa=='S')
        {
            N_CA_C_angle=111.2812;
            N_C_CA_CB_diangle=122.6618;
        }
        else if (aa=='C' || aa=='U')
        {
            N_CA_C_angle= 110.8856;
            N_C_CA_CB_diangle=122.5037;
        }
        else if (aa=='V')
        {
            CA_C_O_angle=120.5686;
            N_CA_C_angle=109.7698;
            N_C_CA_CB_diangle=123.2347;
        }
        else if (aa=='I')
        {
            CA_C_O_angle=120.5403;
            N_CA_C_angle=109.7202;
            N_C_CA_CB_diangle=123.2347;
        }
        else if (aa=='N') 
        {
            CA_C_O_angle=120.4826;
            N_CA_C_angle=111.5;
            N_C_CA_CB_diangle=123.2254;
        }
    }
    else if (aa=='L')
    {
        CA_C_O_angle=120.4647;
        N_CA_C_angle=110.8652;
        N_C_CA_CB_diangle=122.4948;
    }
    else if (aa=='T')
    {
        CA_C_O_angle=120.5359;
        N_CA_C_angle=110.7014;
        N_C_CA_CB_diangle=123.0953;
    }
    else if (aa=='R'||aa=='K'||aa=='O')
    {
        CA_C_O_angle=120.54;
        N_CA_C_angle=111.08;
        N_C_CA_CB_diangle=122.76;
        if (aa=='R') N_CA_C_angle=110.98;
    }
    else if (aa=='D')
    {
        CA_C_O_angle=120.51;
        N_CA_C_angle=111.03;
        N_C_CA_CB_diangle=122.82;
    }
    else if (aa=='Q')
    {
        CA_C_O_angle=120.5029;
        N_CA_C_angle=111.0849;
        N_C_CA_CB_diangle=122.8134;
    }
    else if (aa=='E')
    {
        CA_C_O_angle=120.511;
        N_CA_C_angle=111.1703;
        N_C_CA_CB_diangle=122.8702;
    }
    else if (aa=='M')
    {
        CA_C_O_angle=120.4816;
        N_CA_C_angle=110.9416;
        N_C_CA_CB_diangle=122.6733;
    }
    else if (aa=='H')
    {
        CA_C_O_angle=120.4732;
        N_CA_C_angle=111.0859;
        N_C_CA_CB_diangle=122.6711;
    }
    else if (aa=='F')
    {
        CA_C_O_angle=120.5316;
        N_CA_C_angle=110.7528;
        N_C_CA_CB_diangle=122.6054;
    }
    else if (aa=='Y')
    {
        CA_C_O_angle=120.5434;
        N_CA_C_angle=110.9288;
        N_C_CA_CB_diangle=122.6023;
    }
    else if (aa=='W')
    {
        CA_C_O_angle=120.5117;
        N_CA_C_angle=110.8914;
        N_C_CA_CB_diangle=122.6112;
    }
}


string write_pdc_structure(ModelUnit &pep,string &header)
{
    string txt=header;
    if (pep.chains.size()==0) return txt;
    stringstream buf;
    int a,c,r;
    char moltype;
    int bfactor_mode;
    pdb2fasta(pep);
    int resi_prev=0;
    char icode_prev=' ';
    int L;
    int32_t x,y,z;
    int16_t bfactor;
    int32_t prev_x,prev_y,prev_z;
    int16_t dx16,dy16,dz16;
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
        prev_x=pep.chains[c].residues[0].atoms[0].xyz[0];
        prev_y=pep.chains[c].residues[0].atoms[0].xyz[1];
        prev_z=pep.chains[c].residues[0].atoms[0].xyz[2];
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
void write_pdc_structure(const string &filename,ModelUnit &pep,string &header)
{
    if (filename=="-")
        cout<<write_pdc_structure(pep,header)<<flush;
    else
    {
        string filename_str=filename;
        if (EndsWith(filename,".gz")) filename_str=filename.substr(0,filename.size()-3);
        ofstream fp(filename_str.c_str());
        fp<<write_pdc_structure(pep,header)<<flush;
        fp.close();
        if (EndsWith(filename,".gz"))
            int r=system(((string)("gzip -f "+filename_str)).c_str());
        filename_str.clear();
    }
}

void xyz2double(AtomUnit &atom, double *c1)
{
    c1[0]=atom.xyz[0]*0.001;
    c1[1]=atom.xyz[1]*0.001;
    c1[2]=atom.xyz[2]*0.001;
}

/* lossy compression */
string write_pdc_lossy_structure(ModelUnit &pep,string &header,const int lossy)
{
    string txt=header;
    if (pep.chains.size()==0) return txt;
    stringstream buf;
    int a,c,r;
    char moltype;
    int bfactor_mode;
    pdb2fasta(pep);
    int resi_prev=0;
    char icode_prev=' ';
    int L;
    int32_t x,y,z;
    int16_t bfactor,prev_b;
    int32_t prev_x,prev_y,prev_z;
    int16_t dx16,dy16,dz16,db16;
    int8_t dx8,dy8,dz8,db8;
    double tor;
    int8_t tor8;
    double c1[3],c2[3],c3[3],c4[3];
    string resn;


    AtomUnit atom;
    atom.xyz.assign(3,0);
    atom.name=" CA ";
    atom.bfactor=0;
    ResidueUnit residue;
    ChainUnit chain;

    double **allCA_arr;
    double **allNCAC_arr;
    double RotMatix[3][3];
    double TranVect[3];
    double **target_stru;
    double **refrence_stru;
    double v1[3],v2[3],v3[3];
    for (c=0;c<pep.chains.size();c++)
    {
        moltype=check_moltype(pep.chains[c]);
        bfactor_mode=check_bfactor_mode(pep.chains[c]);
        if (lossy>=3 && bfactor_mode==2) bfactor_mode=1;
        L=pep.chains[c].residues.size();
        if (L<3 && 1<=lossy && lossy<=2)
        {
            buf <<"@a\n>"<<pep.chains[c].chainID<<'\t'
                <<moltype<<'\t'<<bfactor_mode<<'\t'<<L<<'\n'
                <<pep.chains[c].sequence<<'\n';
        }
        buf <<"@a\n>"<<pep.chains[c].chainID<<'\t'
            <<moltype<<'\t'<<bfactor_mode<<'\t'<<L<<'\t'<<lossy<<'\n'
            <<pep.chains[c].sequence<<'\n';
        if (lossy==0)
        {
            NewArray(&allCA_arr,L,3);
            NewArray(&allNCAC_arr,L*5+1,3);
            NewArray(&target_stru,3,3);
            NewArray(&refrence_stru,3,3);
            for (r=0;r<L;r++)
                chain.residues.push_back(pep.chains[c].residues[r]);
        }
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
        if (L<3 && 1<=lossy && lossy<=2)
        {
            prev_x=pep.chains[c].residues[0].atoms[0].xyz[0];
            prev_y=pep.chains[c].residues[0].atoms[0].xyz[1];
            prev_z=pep.chains[c].residues[0].atoms[0].xyz[2];
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
            continue;
        }
        for (r=0;r<L;r++)
        {
            x=pep.chains[c].residues[r].atoms[1].xyz[0]; // CA
            y=pep.chains[c].residues[r].atoms[1].xyz[1];
            z=pep.chains[c].residues[r].atoms[1].xyz[2];
            if (r==0)
            {
                buf.write((char *)&x,sizeof(int32_t));
                buf.write((char *)&y,sizeof(int32_t));
                buf.write((char *)&z,sizeof(int32_t));
                prev_x=x;
                prev_y=y;
                prev_z=z;
            }
            else
            {
                if (lossy<=1 || lossy==3)
                {
                    dx16=(x-prev_x);
                    dy16=(y-prev_y);
                    dz16=(z-prev_z);
                    buf.write((char *)&dx16,sizeof(int16_t));
                    buf.write((char *)&dy16,sizeof(int16_t));
                    buf.write((char *)&dz16,sizeof(int16_t));
                    prev_x+=dx16;
                    prev_y+=dy16;
                    prev_z+=dz16;
                }
                else
                {
                    dx8=(x-prev_x)/100;
                    dy8=(y-prev_y)/100;
                    dz8=(z-prev_z)/100;
                    buf.write((char *)&dx8,sizeof(int8_t));
                    buf.write((char *)&dy8,sizeof(int8_t));
                    buf.write((char *)&dz8,sizeof(int8_t));
                    prev_x+=dx8*100;
                    prev_y+=dy8*100;
                    prev_z+=dz8*100;
                }
            }
            if (lossy==0)
            {
                allCA_arr[r][0]=0.001*prev_x;
                allCA_arr[r][1]=0.001*prev_y;
                allCA_arr[r][2]=0.001*prev_z;
                chain.residues[r].atoms[1].xyz[0]=prev_x;
                chain.residues[r].atoms[1].xyz[1]=prev_y;
                chain.residues[r].atoms[1].xyz[2]=prev_z;
            }
        }
        if (lossy==3)
        {
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
            continue;
        }
        if (lossy>=3) continue;
        else if (lossy==0)
        {
            allNCAC_arr[0][0]=1.46 * cos(deg2rad(110.8914));
            allNCAC_arr[0][1]=allNCAC_arr[0][2]=0;
            allNCAC_arr[1][0]=allNCAC_arr[1][1]=allNCAC_arr[1][2]=0;
            allNCAC_arr[0][1]=1.46 * sin(deg2rad(110.8914));
            allNCAC_arr[2][0]=1.52;
            allNCAC_arr[2][1]=allNCAC_arr[2][2]=0;
            allNCAC_arr[3][0]=allNCAC_arr[3][1]=allNCAC_arr[3][2]=0;
            allNCAC_arr[4][0]=allNCAC_arr[4][1]=allNCAC_arr[4][2]=0;
        }
        /* phi, psi, omega */
        double N_CA_C_angle=111;
        double CA_C_O_angle=120.5;
        double N_C_CA_CB_diangle=122.6;
        double N_CA_C_O_diangle =120.0; // LTRKDEQMHFYW
        string resn;
        for (r=0;r<L;r++)
        {
            get_bb_angle(resn, N_CA_C_angle, CA_C_O_angle,
                N_C_CA_CB_diangle, N_CA_C_O_diangle);
            
            /* N-CA-C-N[+1] or N-CA-C-OXT */
            xyz2double(pep.chains[c].residues[r].atoms[0],c1); //N
            xyz2double(pep.chains[c].residues[r].atoms[1],c2); //CA
            xyz2double(pep.chains[c].residues[r].atoms[2],c3); //C
            if (r+1<L) xyz2double(pep.chains[c].residues[r+1].atoms[0],c4);//N[+1]
            else     xyz2double(pep.chains[c].residues[r].atoms.back(),c4);//OXT
            tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
            tor8=tor*INT8_MAX/180.+0.5;
            buf.write((char *)&tor8,sizeof(int8_t));
            if (lossy==0)
            {
                tor=tor8*180./INT8_MAX;
                calculateCoordinates(allNCAC_arr[(r+1)*5],
                    allNCAC_arr[r*5],allNCAC_arr[r*5+1],allNCAC_arr[r*5+2],
                    1.33,116.642992978143,tor);
            }

            if (r+1<L)
            {
                /* omega: CA-C-N[+1]-CA[+1] */
                tor=180;
                if (lossy==1)
                {
                    xyz2double(pep.chains[c].residues[r].atoms[1],c1); //CA
                    xyz2double(pep.chains[c].residues[r].atoms[2],c2); //C
                    xyz2double(pep.chains[c].residues[r+1].atoms[0],c3); //N[+1]
                    xyz2double(pep.chains[c].residues[r+1].atoms[1],c4); //CA[+1]
                    tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
                    tor8=tor*INT8_MAX/180.+0.5;
                    buf.write((char *)&tor8,sizeof(int8_t));
                    tor=tor8*180./INT8_MAX;
                }
                if (lossy==0) calculateCoordinates(allNCAC_arr[(r+1)*5+1],
                    allNCAC_arr[r*5+1],allNCAC_arr[r*5+2],allNCAC_arr[(r+1)*5],
                    1.46,121.382215820277,tor);

                /* C-N[+1]-CA[+1]-C[+1] */
                xyz2double(pep.chains[c].residues[r].atoms[2],c1); //C
                xyz2double(pep.chains[c].residues[r+1].atoms[0],c2); //N[+1}
                xyz2double(pep.chains[c].residues[r+1].atoms[1],c3); //CA[+1]
                xyz2double(pep.chains[c].residues[r+1].atoms[2],c4); //C[+1]
                tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
                tor8=tor*INT8_MAX/180.+0.5;
                buf.write((char *)&tor8,sizeof(int8_t));
                if (lossy==0)
                {
                    tor=tor8*180./INT8_MAX;
                    calculateCoordinates(allNCAC_arr[(r+1)*5+2],
                        allNCAC_arr[r*5+2],allNCAC_arr[(r+1)*5],allNCAC_arr[(r+1)*5+1],
                        1.52,N_CA_C_angle,tor);
                }
            }

            if (lossy) continue;
            /* N-C-CA-CB */
            calculateCoordinates(allNCAC_arr[r*5+3],
                allNCAC_arr[r*5],allNCAC_arr[r*5+2],allNCAC_arr[r*5+1],
                1.52,109.5,N_C_CA_CB_diangle);
            
            /* N[+1]-CA-C-O */
            calculateCoordinates(allNCAC_arr[r*5+4],
                allNCAC_arr[(r+1)*5],allNCAC_arr[r*5+1],allNCAC_arr[r*5+2],
                1.23,CA_C_O_angle,180);
        }
        /* coordinate of N, C, O, CB, OXT */
        if (lossy==0)
        {
            int i,j;
            v3[0]=v3[1]=v3[2]=0;
            for (r=1;r<L-1;r++)
            {
                for (i=0;i<3;i++)
                {
                    refrence_stru[i][0]=allCA_arr[i+r-1][0];
                    refrence_stru[i][1]=allCA_arr[i+r-1][1];
                    refrence_stru[i][2]=allCA_arr[i+r-1][2];

                    target_stru[i][0]=allNCAC_arr[(i+r-1)*5+1][0];
                    target_stru[i][1]=allNCAC_arr[(i+r-1)*5+1][1];
                    target_stru[i][2]=allNCAC_arr[(i+r-1)*5+1][2];
                }
                RotateCoor(target_stru, refrence_stru, 3, RotMatix, TranVect);
                if (r==1) /* update first residue */
                {
                    v1[0]=target_stru[0][0];
                    v1[1]=target_stru[0][1];
                    v1[2]=target_stru[0][2];
                    ChangeCoor(v1, RotMatix, TranVect, v2);
                    v3[0]=refrence_stru[0][0]-v2[0];
                    v3[1]=refrence_stru[0][1]-v2[1];
                    v3[2]=refrence_stru[0][2]-v2[2];
                    for (a=0;a<5;a++)
                    {
                        //if (a==1 || (a==3 && chain.residues[0].resn=="GLY"))
                        if ((a==3 && chain.residues[0].resn=="GLY"))
                            continue;
                        v1[0]=allNCAC_arr[a][0];
                        v1[1]=allNCAC_arr[a][1];
                        v1[2]=allNCAC_arr[a][2];
                        ChangeCoor(v1, RotMatix, TranVect, v2);
                        if (a==4 && chain.residues[0].resn=="GLY")
                        {
                            chain.residues[0].atoms[3].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[0].atoms[3].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[0].atoms[3].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                        else
                        {
                            chain.residues[0].atoms[a].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[0].atoms[a].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[0].atoms[a].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                    }
                }
                if (r+1==L-1) /* update last residue */
                {
                    v1[0]=target_stru[2][0];
                    v1[1]=target_stru[2][1];
                    v1[2]=target_stru[2][2];
                    ChangeCoor(v1, RotMatix, TranVect, v2);
                    v3[0]=refrence_stru[2][0]-v2[0];
                    v3[1]=refrence_stru[2][1]-v2[1];
                    v3[2]=refrence_stru[2][2]-v2[2];
                    for (a=0;a<6;a++)
                    {
                        if (a==1 || (a==3 && chain.residues[r+1].resn=="GLY"))
                            continue;
                        v1[0]=allNCAC_arr[(r+1)*5+a][0];
                        v1[1]=allNCAC_arr[(r+1)*5+a][1];
                        v1[2]=allNCAC_arr[(r+1)*5+a][2];
                        ChangeCoor(v1, RotMatix, TranVect, v2);
                        if (a==4 && chain.residues[r+1].resn=="GLY")
                        {
                            chain.residues[r+1].atoms[3].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[r+1].atoms[3].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[r+1].atoms[3].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                        else if (a==5)
                        {
                            chain.residues[r+1].atoms.back().xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[r+1].atoms.back().xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[r+1].atoms.back().xyz[2]=1000*(v2[2]+v3[2]);
                        }
                        else
                        {
                            chain.residues[r+1].atoms[a].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[r+1].atoms[a].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[r+1].atoms[a].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                    }
                }

                /* update current residue */
                v1[0]=target_stru[1][0];
                v1[1]=target_stru[1][1];
                v1[2]=target_stru[1][2];
                ChangeCoor(v1, RotMatix, TranVect, v2);
                TranVect[0]+=refrence_stru[1][0]-v2[0];
                TranVect[1]+=refrence_stru[1][1]-v2[1];
                TranVect[2]+=refrence_stru[1][2]-v2[2];
                for (a=0;a<5;a++)
                {
                    if (a==1 || (a==3 && chain.residues[r].resn=="GLY"))
                        continue;
                    v1[0]=allNCAC_arr[r*5+a][0];
                    v1[1]=allNCAC_arr[r*5+a][1];
                    v1[2]=allNCAC_arr[r*5+a][2];
                    ChangeCoor(v1, RotMatix, TranVect, v2);
                    if (a==4 && chain.residues[r].resn=="GLY")
                    {
                        chain.residues[r].atoms[3].xyz[0]=1000*v2[0];
                        chain.residues[r].atoms[3].xyz[1]=1000*v2[1];
                        chain.residues[r].atoms[3].xyz[2]=1000*v2[2];
                    }
                    else
                    {
                        chain.residues[r].atoms[a].xyz[0]=1000*v2[0];
                        chain.residues[r].atoms[a].xyz[1]=1000*v2[1];
                        chain.residues[r].atoms[a].xyz[2]=1000*v2[2];
                    }
                }
            }
            if (r+1==L-1) /* update last residue */
            {
                v1[0]=target_stru[2][0];
                v1[1]=target_stru[2][1];
                v1[2]=target_stru[2][2];
                ChangeCoor(v1, RotMatix, TranVect, v2);
                v3[0]=refrence_stru[2][0]-v2[0];
                v3[1]=refrence_stru[2][1]-v2[1];
                v3[2]=refrence_stru[2][2]-v2[2];
                for (a=0;a<6;a++)
                {
                    if (a==1 || (a==3 && chain.residues[r+1].resn=="GLY"))
                        continue;
                    v1[0]=allNCAC_arr[(r+1)*5+a][0];
                    v1[1]=allNCAC_arr[(r+1)*5+a][1];
                    v1[2]=allNCAC_arr[(r+1)*5+a][2];
                    ChangeCoor(v1, RotMatix, TranVect, v2);
                    if (a==4 && chain.residues[r+1].resn=="GLY")
                    {
                        chain.residues[r+1].atoms[3].xyz[0]=1000*(v2[0]+v3[0]);
                        chain.residues[r+1].atoms[3].xyz[1]=1000*(v2[1]+v3[1]);
                        chain.residues[r+1].atoms[3].xyz[2]=1000*(v2[2]+v3[2]);
                    }
                    else if (a==5)
                    {
                        chain.residues[r+1].atoms.back().xyz[0]=1000*(v2[0]+v3[0]);
                        chain.residues[r+1].atoms.back().xyz[1]=1000*(v2[1]+v3[1]);
                        chain.residues[r+1].atoms.back().xyz[2]=1000*(v2[2]+v3[2]);
                    }
                    else
                    {
                        chain.residues[r+1].atoms[a].xyz[0]=1000*(v2[0]+v3[0]);
                        chain.residues[r+1].atoms[a].xyz[1]=1000*(v2[1]+v3[1]);
                        chain.residues[r+1].atoms[a].xyz[2]=1000*(v2[2]+v3[2]);
                    }
                }
            }
        }
        /* side chain torsion */
        for (r=0;r<L;r++)
        {
            resn=pep.chains[c].residues[r].resn;
            if (resn=="GLY" || resn=="ALA") continue; //2

            /* chi-1 */
            xyz2double(pep.chains[c].residues[r].atoms[0],c1); //N
            xyz2double(pep.chains[c].residues[r].atoms[1],c2); //CA
            xyz2double(pep.chains[c].residues[r].atoms[3],c3); //CB
            if (resn=="THR") xyz2double(pep.chains[c].residues[r].atoms[6],c4); //OG1
            else             xyz2double(pep.chains[c].residues[r].atoms[5],c4); //CG
            tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
            tor8=tor*INT8_MAX/180.+0.5;
            buf.write((char *)&tor8,sizeof(int8_t));
            if (lossy==0)
            {
                tor=tor8*180./INT8_MAX;
                if (resn=="SER")
                    calculateCoordinates(c4,c1,c2,c3,1.417,110.773,tor);
                else if (resn=="CYS")
                    calculateCoordinates(c4,c1,c2,c3,1.808,113.8169,tor);
                else if (resn=="VAL" || resn=="ILE")
                    calculateCoordinates(c4,c1,c2,c3,1.527,110.7,tor);
                else if (resn=="LEU")
                    calculateCoordinates(c4,c1,c2,c3,1.53,116.10,tor);
                else if (resn=="THR")
                    calculateCoordinates(c4,c1,c2,c3,1.43,109.18,tor); //OG1
                else if (resn=="ARG" || resn=="LYS")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.83,tor);
                else if (resn=="ASP")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.06,tor);
                else if (resn=="ASN")
                    calculateCoordinates(c4,c1,c2,c3,1.52,112.62,tor);
                else if (resn=="GLU")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.82,tor);
                else if (resn=="GLN")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.75,tor);
                else if (resn=="MET")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.68,tor);
                else if (resn=="HIS")
                    calculateCoordinates(c4,c1,c2,c3,1.49,113.74,tor);
                else if (resn=="PRO")
                    calculateCoordinates(c4,c1,c2,c3,1.49,104.21,tor);
                else if (resn=="PHE")
                    calculateCoordinates(c4,c1,c2,c3,1.50,113.85,tor);
                else if (resn=="TYR")
                    calculateCoordinates(c4,c1,c2,c3,1.51,113.8,tor);
                else if (resn=="TRP")
                    calculateCoordinates(c4,c1,c2,c3,1.50,114.10,tor);
                if (resn=="THR")
                {
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //OG1
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //OG1
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //OG1
                    calculateCoordinates(c4,c1,c2,c3,1.53,111.13,tor-120);
                }
                chain.residues[r].atoms[5].xyz[0]=1000.*c4[0]; //CG
                chain.residues[r].atoms[5].xyz[1]=1000.*c4[1]; //CG
                chain.residues[r].atoms[5].xyz[2]=1000.*c4[2]; //CG
                if (resn=="ILE")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.527,110.7,tor-120);
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //CG2
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //CG2
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //CG2
                }
                else if (resn=="VAL")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.527,110.7,tor+120);
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //CG2
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //CG2
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //CG2
                }
            }
            if (resn=="CYS" || resn=="SER" || resn=="THR" || resn=="VAL" )
                continue; // 2+4=6

            /* chi-2: CA-CB-CG-XD: ARG LYS MET GLU GLN ILE LEU
             * HIS TYR TRP PHE PRO ASP ASN */
            xyz2double(pep.chains[c].residues[r].atoms[1],c1); //CA
            xyz2double(pep.chains[c].residues[r].atoms[3],c2); //CB
            xyz2double(pep.chains[c].residues[r].atoms[5],c3); //CG
            if (resn=="ILE") xyz2double(pep.chains[c].residues[r].atoms[7],c4); //CD1
            else             xyz2double(pep.chains[c].residues[r].atoms[6],c4); //XD
            tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
            tor8=tor*INT8_MAX/180.+0.5;
            buf.write((char *)&tor8,sizeof(int8_t));
            if (lossy==0)
            {
                tor=tor8*180./INT8_MAX;
                if (resn=="ILE")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.97,tor);
                else if (resn=="LEU")
                    calculateCoordinates(c4,c1,c2,c3,1.524,110.27,tor);
                else if (resn=="ARG" || resn=="LYS")
                    calculateCoordinates(c4,c1,c2,c3,1.52,111.79,tor);
                else if (resn=="ASP")
                    calculateCoordinates(c4,c1,c2,c3,1.25,119.22,tor);//OD1
                else if (resn=="ASN")
                    calculateCoordinates(c4,c1,c2,c3,1.23,120.85,tor);//OD1
                else if (resn=="GLU")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.31,tor);
                else if (resn=="GLN")
                    calculateCoordinates(c4,c1,c2,c3,1.52,112.78,tor);
                else if (resn=="MET")
                    calculateCoordinates(c4,c1,c2,c3,1.81,112.69,tor);
                else if (resn=="HIS")
                    calculateCoordinates(c4,c1,c2,c3,1.35,130.61,tor); //CD2
                else if (resn=="PRO")
                    calculateCoordinates(c4,c1,c2,c3,1.50,105.03,tor);
                else if (resn=="PHE")
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,tor); //CD1
                else if (resn=="TYR")
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.98,tor); //CD1
                else if (resn=="TRP")
                    calculateCoordinates(c4,c1,c2,c3,1.37,127.07,tor); //CD1
                
                if (resn=="ILE")
                {
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD1
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD1
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD1
                }
                else
                {
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //XD
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //XD
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //XD
                }
                if (resn=="LEU")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.525,110.58,tor+120);//CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2
                }
                else if (resn=="ASP")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.25,118.218,tor+180);//OD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //OD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //OD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //OD2
                }
                else if (resn=="ASN")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.33,116.48,tor+180);//ND2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //ND2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //ND2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //ND2
                }
                else if (resn=="HIS")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.38,122.85,tor+180);//ND1
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //ND1
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //ND1
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //ND1
                    
                    // CB_CG_ND1_CE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[7],c3); //ND1
                    calculateCoordinates(c4,c1,c2,c3,1.32,108.5,180);//CE1
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE1
                    
                    // CB_CG_CD2_NE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[6],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.35,108.5,180);//NE2
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //NE2
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //NE2
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //NE2
                }
                else if (resn=="PHE")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,tor+180); //CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2

                    // CB_CG_CD1_CE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[6],c3); //CD1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE1
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE1

                    // CB_CG_CD2_CE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[7],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE2
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //CE2
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //CE2
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //CE2
                    
                    // CG_CD1_CE1_CZ_diangle = 0.0
                    xyz2double(chain.residues[r].atoms[5],c1); //CG
                    xyz2double(chain.residues[r].atoms[6],c2); //CD1
                    xyz2double(chain.residues[r].atoms[8],c3); //CE1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,0); //CZ
                    chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //CZ
                    chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //CZ
                    chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //CZ
                }
                else if (resn=="TYR")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.82,tor+180); //CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2

                    // CB_CG_CD1_CE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[6],c3); //CD1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE1
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE1

                    // CB_CG_CD2_CE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[7],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE2
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //CE2
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //CE2
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //CE2
                    
                    // CG_CD1_CE1_CZ_diangle = 0.0
                    xyz2double(chain.residues[r].atoms[5],c1); //CG
                    xyz2double(chain.residues[r].atoms[6],c2); //CD1
                    xyz2double(chain.residues[r].atoms[8],c3); //CE1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,0); //CZ
                    chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //CZ
                    chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //CZ
                    chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //CZ

                    // CD1_CE1_CZ_OH_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[6],c1); //CD1
                    xyz2double(chain.residues[r].atoms[8],c2); //CE1
                    xyz2double(chain.residues[r].atoms[10],c3); //CZ
                    calculateCoordinates(c4,c1,c2,c3,1.39,119.78,180); //OH
                    chain.residues[r].atoms[11].xyz[0]=1000.*c4[0]; //OH
                    chain.residues[r].atoms[11].xyz[1]=1000.*c4[1]; //OH
                    chain.residues[r].atoms[11].xyz[2]=1000.*c4[2]; //OH
                }
                else if (resn=="TRP")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.43,126.66,tor+180); //CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2

                    // CB_CG_CD1_NE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[6],c3); //CD1
                    calculateCoordinates(c4,c1,c2,c3,1.38,108.5,180); //NE1
                    chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //NE1
                    chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //NE1
                    chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //NE1

                    // CB_CG_CD2_CE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[7],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.40,108.5,180); //CE2
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE2
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE2
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE2

                    // CB_CG_CD2_CE3_diangle = 0.0
                    calculateCoordinates(c4,c1,c2,c3,1.40,133.83,0); //CE3
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //CE3
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //CE3
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //CE3

                    // CG_CD2_CE2_CZ2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[5],c1); //CG
                    xyz2double(chain.residues[r].atoms[7],c2); //CD2
                    xyz2double(chain.residues[r].atoms[8],c3); //CE2
                    calculateCoordinates(c4,c1,c2,c3,1.40,120.0,180); //CZ2
                    chain.residues[r].atoms[12].xyz[0]=1000.*c4[0]; //CZ2
                    chain.residues[r].atoms[12].xyz[1]=1000.*c4[1]; //CZ2
                    chain.residues[r].atoms[12].xyz[2]=1000.*c4[2]; //CZ2

                    // CG_CD2_CE3_CZ3_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[9],c3); //CE3
                    calculateCoordinates(c4,c1,c2,c3,1.40,120.0,180); //CZ3
                    chain.residues[r].atoms[13].xyz[0]=1000.*c4[0]; //CZ3
                    chain.residues[r].atoms[13].xyz[1]=1000.*c4[1]; //CZ3
                    chain.residues[r].atoms[13].xyz[2]=1000.*c4[2]; //CZ3

                    // CD2_CE2_CZ2_CH2_diangle = 0.0
                    xyz2double(chain.residues[r].atoms[7],c1); //CD2
                    xyz2double(chain.residues[r].atoms[8],c2); //CE2
                    xyz2double(chain.residues[r].atoms[12],c3); //CZ2
                    calculateCoordinates(c4,c1,c2,c3,1.40,120.0,0); //CH2
                    chain.residues[r].atoms[11].xyz[0]=1000.*c4[0]; //CH2
                    chain.residues[r].atoms[11].xyz[1]=1000.*c4[1]; //CH2
                    chain.residues[r].atoms[11].xyz[2]=1000.*c4[2]; //CH2
                }
            }
            if (resn=="HIS" || resn=="ILE" || resn=="PHE" || resn=="LEU" ||
                resn=="TRP" || resn=="PRO" || resn=="ASN" || resn=="TYR" ||
                resn=="ASP") continue; // 2+4+9=15

            /* chi-3: CB-CG-CD-XE: ARG LYS GLN GLU MET */
            xyz2double(pep.chains[c].residues[r].atoms[3],c1); //CB
            xyz2double(pep.chains[c].residues[r].atoms[5],c2); //CG
            xyz2double(pep.chains[c].residues[r].atoms[6],c3); //CD
            if (resn=="GLN") xyz2double(pep.chains[c].residues[r].atoms[8],c4); //OE1
            else             xyz2double(pep.chains[c].residues[r].atoms[7],c4); //XE
            tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
            tor8=tor*INT8_MAX/180.+0.5;
            buf.write((char *)&tor8,sizeof(int8_t));
            if (lossy==0)
            {
                tor=tor8*180./INT8_MAX;
                if (resn=="ARG"||resn=="LYS")
                    calculateCoordinates(c4,c1,c2,c3,1.46,111.68,tor);
                else if (resn=="GLU")
                    calculateCoordinates(c4,c1,c2,c3,1.25,119.02,tor);
                else if (resn=="GLN")
                    calculateCoordinates(c4,c1,c2,c3,1.24,120.86,tor);
                else if (resn=="MET")
                    calculateCoordinates(c4,c1,c2,c3,1.79,100.61,tor);
                
                if (resn=="GLN")
                {
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //OE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //OE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //OE1
                    calculateCoordinates(c4,c1,c2,c3,1.33,116.50,tor+180);
                }
                chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //XE
                chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //XE
                chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //XE
                if (resn=="GLU")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.25,119.02,tor+180);
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //OE2
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //OE2
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //OE2
                }
            }
            if (resn=="GLN" || resn=="GLU" || resn=="MET") continue; // 2+4+9+3=18
            
            /* chi-4: CG-CD-XE-XZ: ARG LYS */
            xyz2double(pep.chains[c].residues[r].atoms[5],c1); //CG
            xyz2double(pep.chains[c].residues[r].atoms[6],c2); //CD
            xyz2double(pep.chains[c].residues[r].atoms[7],c3); //XE
            if (resn=="LYS") xyz2double(pep.chains[c].residues[r].atoms[8],c4); //NZ
            else            xyz2double(pep.chains[c].residues[r].atoms[10],c4); //CZ
            tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
            tor8=tor*INT8_MAX/180.+0.5;
            buf.write((char *)&tor8,sizeof(int8_t));
            if (lossy==0)
            {
                tor=tor8*180./INT8_MAX;
                calculateCoordinates(c4,c1,c2,c3,1.33,124.79,tor);
                if (resn=="LYS")
                {
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //NZ
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //NZ
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //NZ
                }
                else
                {
                    chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //CZ
                    chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //CZ
                    chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //CZ
                }
            }
            if (resn=="LYS") continue; // 2+4+9+3+1=19
            
            /* chi-5: CD-NE-CZ-NH1: ARG */
            xyz2double(pep.chains[c].residues[r].atoms[6],c1); //CD
            xyz2double(pep.chains[c].residues[r].atoms[7],c2); //NE
            xyz2double(pep.chains[c].residues[r].atoms[10],c3); //CZ
            xyz2double(pep.chains[c].residues[r].atoms[8],c4); //NH1
            tor=rad2deg(Points2Dihedral(c1, c2, c3, c4));
            tor8=tor*INT8_MAX/180.+0.5;
            buf.write((char *)&tor8,sizeof(int8_t));
            if (lossy==0)
            {
                tor=tor8*180./INT8_MAX;
                calculateCoordinates(c4,c1,c2,c3,1.33,120.64,tor);
                chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //NH1
                chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //NH1
                chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //NH1
                calculateCoordinates(c4,c1,c2,c3,1.33,120.64,tor+180);
                chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //NH2
                chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //NH2
                chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //NH2
            }
        }
        if (lossy==0)
        {
            resn.clear();
            DeleteArray(&allCA_arr,L);
            DeleteArray(&allNCAC_arr,L*5+1);
            DeleteArray(&target_stru,3);
            DeleteArray(&refrence_stru,3);

            vector<int16_t> d16_vec;
            vector<int8_t> d8_vec;
            int16_t d16;
            int8_t d8;
            int i;
            for (r=0;r<L;r++)
            {
                for (a=0;a<chain.residues[r].atoms.size();a++)
                {
                    if (a==1) continue;
                    for (i=0;i<3;i++)
                    {
                        d16=pep.chains[c].residues[r].atoms[a].xyz[i]-
                                    chain.residues[r].atoms[a].xyz[i];
                        d8=0;
                        if (INT8_MAX<=d16 && d16<=INT8_MAX) d8=d16;
                        d8_vec.push_back(d8);
                        if (d8==0) d16_vec.push_back(d16);
                    }
                }
            }
            size_t d16_count=d16_vec.size();
            buf.write((char *)&d16_count,sizeof(size_t));
            for (i=0;i<d16_count;i++)
            {
                d16=d16_vec[i];
                buf.write((char *)&d16,sizeof(int16_t));
            }
            for (i=0;i<d8_vec.size();i++)
            {
                d8=d8_vec[i];
                buf.write((char *)&d8,sizeof(int8_t));
            }
            for (r=0;r<L;r++)
                for (a=0;a<chain.residues[r].atoms.size();a++)
                    chain.residues[r].atoms[a].xyz.clear();
            chain.residues.clear();
            d16_vec.clear();
            d8_vec.clear();
        }
        prev_b=pep.chains[c].residues[0].atoms[0].bfactor;
        buf.write((char *)&prev_b,sizeof(int16_t));
        if (bfactor_mode>=1)
        {
            vector<int16_t> bfactor_vec(L,prev_b);
            for (r=0;r<L;r++)
            {
                if (r>0) prev_b=bfactor_vec[r-1];
                bfactor=pep.chains[c].residues[r].atoms[1].bfactor; // CA
                db16=0.1*(bfactor-prev_b);
                if (db16>=INT8_MAX) db8=INT8_MAX;
                else if (db16<=INT8_MIN) db8=INT8_MIN;
                else db8=db16;
                buf.write((char *)&db8,sizeof(int8_t));
                bfactor_vec[r]=prev_b+10*db8;
            }
            if (bfactor_mode==2)
            {
                for (r=0;r<L;r++)
                {
                    bfactor=bfactor_vec[r];
                    for (a=0;a<pep.chains[c].residues[r].atoms.size();a++)
                    {
                        if (a==1) continue;
                        db16=0.1*(pep.chains[c].residues[r].atoms[a].bfactor-bfactor+0.5);
                        if (db16>=INT8_MAX) db8=INT8_MAX;
                        else if (db16<=INT8_MIN) db8=INT8_MIN;
                        else db8=db16;
                        buf.write((char *)&db8,sizeof(int8_t));
                    }
                }
            }
            bfactor_vec.clear();
        }
    }
    buf<<flush;
    txt+=buf.str()+"\nEND\n";
    buf.str(string());
    return txt;
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_pdc_lossy_structure(const string &filename,ModelUnit &pep,string &header,const int lossy)
{
    if (filename=="-")
        cout<<write_pdc_lossy_structure(pep,header,lossy)<<flush;
    else
    {
        string filename_str=filename;
        if (EndsWith(filename,".gz")) filename_str=filename.substr(0,filename.size()-3);
        ofstream fp(filename_str.c_str());
        fp<<write_pdc_lossy_structure(pep,header,lossy)<<flush;
        fp.close();
        if (EndsWith(filename,".gz"))
            int r=system(((string)("gzip -f "+filename_str)).c_str());
        filename_str.clear();
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

int read_pdc_structure(const char *filename,ModelUnit &pep,string &header,
    map<string, vector<string> > & ordMapR)
{
    ChainUnit chain;
    ResidueUnit residue;
    AtomUnit atom;
    atom.xyz.assign(3,0);
    atom.bfactor=0;
    atom.name=" CA ";
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
        fp_gz.open("gunzip -c "+filename_str);
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
    int16_t dx16,dy16,dz16;
    int8_t dx8,dy8,dz8,db8;
    int atomNum;
    int lossy=0;
    double tor,ang,len;
    int8_t tor8;
    char aa;
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
        lossy=-1;
        if (line_vec.size()>4) lossy=atoi(line_vec[4].c_str());
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
        if (lossy>=0)
        {
            double **allCA_arr;
            double **allNCAC_arr;
            NewArray(&allCA_arr,L,3);
            NewArray(&allNCAC_arr,L*5+1,3);
            double RotMatix[3][3];
            double TranVect[3];
            double **target_stru;
            double **refrence_stru;
            NewArray(&target_stru,3,3);
            NewArray(&refrence_stru,3,3);
            double v1[3],v2[3],v3[3];
            double c1[3],c2[3],c3[3],c4[3];
            for (r=0;r<L;r++)
            {
                residue.resi=resi_vec[r];
                residue.icode=icode_vec[r];
                residue.resn=aa1to3(sequence[r]);
                if (lossy>=3) residue.atoms.push_back(atom);
                else
                {
                    for (a=0;a<ordMapR[residue.resn].size()+(r+1==L);a++)
                    {
                        atom.name=" OXT";
                        if (a<ordMapR[residue.resn].size())
                            atom.name=ordMapR[residue.resn][a];
                        residue.atoms.push_back(atom);
                    }
                }
                if (r==0)
                {
                    if (use_stdin)
                    {
                        cin.read((char *)&x,sizeof(int32_t));
                        cin.read((char *)&y,sizeof(int32_t));
                        cin.read((char *)&z,sizeof(int32_t));
                    }
                    else if (use_pstream)
                    {
                        fp_gz.read((char *)&x,sizeof(int32_t));
                        fp_gz.read((char *)&y,sizeof(int32_t));
                        fp_gz.read((char *)&z,sizeof(int32_t));
                    }
                    else
                    {
                        fp.read((char *)&x,sizeof(int32_t));
                        fp.read((char *)&y,sizeof(int32_t));
                        fp.read((char *)&z,sizeof(int32_t));
                    }
                }
                else
                {
                    if (lossy<=1||lossy==3)
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
                        x+=dx16;
                        y+=dy16;
                        z+=dz16;
                    }
                    else
                    {
                        if (use_stdin)
                        {
                            cin.read((char *)&dx8,sizeof(int8_t));
                            cin.read((char *)&dy8,sizeof(int8_t));
                            cin.read((char *)&dz8,sizeof(int8_t));
                        }
                        else if (use_pstream)
                        {
                            fp_gz.read((char *)&dx8,sizeof(int8_t));
                            fp_gz.read((char *)&dy8,sizeof(int8_t));
                            fp_gz.read((char *)&dz8,sizeof(int8_t));
                        }
                        else
                        {
                            fp.read((char *)&dx8,sizeof(int8_t));
                            fp.read((char *)&dy8,sizeof(int8_t));
                            fp.read((char *)&dz8,sizeof(int8_t));
                        }
                        x+=dx8*100;
                        y+=dy8*100;
                        z+=dz8*100;
                    }
                }
                if (lossy>=3)
                {
                    residue.atoms[0].xyz[0]=x;
                    residue.atoms[0].xyz[1]=y;
                    residue.atoms[0].xyz[2]=z;
                }
                else
                {
                    residue.atoms[1].xyz[0]=x;
                    residue.atoms[1].xyz[1]=y;
                    residue.atoms[1].xyz[2]=z;
                }
                chain.residues.push_back(residue);
                residue.atoms.clear();
                allCA_arr[r][0]=0.001*x;
                allCA_arr[r][1]=0.001*y;
                allCA_arr[r][2]=0.001*z;
            }
            if (lossy>=3)
            {
                if (lossy==3)
                {
                    if (bfactor_mode==0)
                    {
                        if (use_stdin)        cin.read((char *)&bfactor,sizeof(int16_t));
                        else if (use_pstream) fp_gz.read((char *)&bfactor,sizeof(int16_t));
                        else                  fp.read((char *)&bfactor,sizeof(int16_t));
                        for (r=0;r<L;r++)
                            for (a=0;a<chain.residues[r].atoms.size();a++)
                                chain.residues[r].atoms[a].bfactor=bfactor;
                    }
                    else
                    {
                        for (r=0;r<L;r++)
                        {
                            if (use_stdin)        cin.read((char *)&bfactor,sizeof(int16_t));
                            else if (use_pstream) fp_gz.read((char *)&bfactor,sizeof(int16_t));
                            else                  fp.read((char *)&bfactor,sizeof(int16_t));
                            for (a=0;a<chain.residues[r].atoms.size();a++)
                                chain.residues[r].atoms[a].bfactor=bfactor;
                        }
                    }
                }
                DeleteArray(&allCA_arr,L);
                DeleteArray(&allNCAC_arr,L*5+1);
                DeleteArray(&target_stru,3);
                DeleteArray(&refrence_stru,3);
                pep.chains.push_back(chain);
                chain.residues.clear();
                continue;
            }
            allNCAC_arr[0][0]=1.46 * cos(deg2rad(110.8914));
            allNCAC_arr[0][1]=allNCAC_arr[0][2]=0;
            allNCAC_arr[1][0]=allNCAC_arr[1][1]=allNCAC_arr[1][2]=0;
            allNCAC_arr[0][1]=1.46 * sin(deg2rad(110.8914));
            allNCAC_arr[2][0]=1.52;
            allNCAC_arr[2][1]=allNCAC_arr[2][2]=0;
            allNCAC_arr[3][0]=allNCAC_arr[3][1]=allNCAC_arr[3][2]=0;
            allNCAC_arr[4][0]=allNCAC_arr[4][1]=allNCAC_arr[4][2]=0;
            /* phi, psi, omega */
            double N_CA_C_angle=111;
            double CA_C_O_angle=120.5;
            double N_C_CA_CB_diangle=122.6;
            double N_CA_C_O_diangle =120.0; // LTRKDEQMHFYW
            string resn;
            for (r=0;r<L;r++)
            {
                get_bb_angle(resn, N_CA_C_angle, CA_C_O_angle,
                    N_C_CA_CB_diangle, N_CA_C_O_diangle);

                /* N-CA-C-N[+1] or N-CA-C-OXT */
                if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                else                     fp.read((char *)&tor8,sizeof(int8_t));
                tor=tor8*180./INT8_MAX;
                calculateCoordinates(allNCAC_arr[(r+1)*5],
                    allNCAC_arr[r*5],allNCAC_arr[r*5+1],allNCAC_arr[r*5+2],
                    1.33,116.642992978143,tor);

                if (r+1<L)
                {
                    /* CA-C-N[+1]-CA[+1] */
                    tor=180;
                    if (lossy==1)
                    {
                        if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                        else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                        else                     fp.read((char *)&tor8,sizeof(int8_t));
                        tor=tor8*180./INT8_MAX;
                    }
                    calculateCoordinates(allNCAC_arr[(r+1)*5+1],
                        allNCAC_arr[r*5+1],allNCAC_arr[r*5+2],allNCAC_arr[(r+1)*5],
                        1.46,121.382215820277,tor);

                    /* C-N[+1]-CA[+1]-C[+1] */
                    if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                    else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                    else                     fp.read((char *)&tor8,sizeof(int8_t));
                    tor=tor8*180./INT8_MAX;
                    calculateCoordinates(allNCAC_arr[(r+1)*5+2],
                        allNCAC_arr[r*5+2],allNCAC_arr[(r+1)*5],allNCAC_arr[(r+1)*5+1],
                        1.52,N_CA_C_angle,tor);
                }
                
                /* N-C-CA-CB */
                calculateCoordinates(allNCAC_arr[r*5+3],
                    allNCAC_arr[r*5],allNCAC_arr[r*5+2],allNCAC_arr[r*5+1],
                    1.52,109.5,N_C_CA_CB_diangle);
                
                /* N[+1]-CA-C-O */
                calculateCoordinates(allNCAC_arr[r*5+4],
                    allNCAC_arr[(r+1)*5],allNCAC_arr[r*5+1],allNCAC_arr[r*5+2],
                    1.23,CA_C_O_angle,180);
            }
            /* coordinate of N, C, O, CB, OXT */
            int i,j;
            v3[0]=v3[1]=v3[2]=0;
            for (r=1;r<L-1;r++)
            {
                for (i=0;i<3;i++)
                {
                    refrence_stru[i][0]=allCA_arr[i+r-1][0];
                    refrence_stru[i][1]=allCA_arr[i+r-1][1];
                    refrence_stru[i][2]=allCA_arr[i+r-1][2];

                    target_stru[i][0]=allNCAC_arr[(i+r-1)*5+1][0];
                    target_stru[i][1]=allNCAC_arr[(i+r-1)*5+1][1];
                    target_stru[i][2]=allNCAC_arr[(i+r-1)*5+1][2];
                }
                RotateCoor(target_stru, refrence_stru, 3, RotMatix, TranVect);
                if (r==1) /* update first residue */
                {
                    v1[0]=target_stru[0][0];
                    v1[1]=target_stru[0][1];
                    v1[2]=target_stru[0][2];
                    ChangeCoor(v1, RotMatix, TranVect, v2);
                    v3[0]=refrence_stru[0][0]-v2[0];
                    v3[1]=refrence_stru[0][1]-v2[1];
                    v3[2]=refrence_stru[0][2]-v2[2];
                    for (a=0;a<5;a++)
                    {
                        //if (a==1 || (a==3 && chain.residues[0].resn=="GLY"))
                        if ((a==3 && chain.residues[0].resn=="GLY"))
                            continue;
                        v1[0]=allNCAC_arr[a][0];
                        v1[1]=allNCAC_arr[a][1];
                        v1[2]=allNCAC_arr[a][2];
                        ChangeCoor(v1, RotMatix, TranVect, v2);
                        if (a==4 && chain.residues[0].resn=="GLY")
                        {
                            chain.residues[0].atoms[3].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[0].atoms[3].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[0].atoms[3].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                        else
                        {
                            chain.residues[0].atoms[a].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[0].atoms[a].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[0].atoms[a].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                    }
                }
                if (r+1==L-1) /* update last residue */
                {
                    v1[0]=target_stru[2][0];
                    v1[1]=target_stru[2][1];
                    v1[2]=target_stru[2][2];
                    ChangeCoor(v1, RotMatix, TranVect, v2);
                    v3[0]=refrence_stru[2][0]-v2[0];
                    v3[1]=refrence_stru[2][1]-v2[1];
                    v3[2]=refrence_stru[2][2]-v2[2];
                    for (a=0;a<6;a++)
                    {
                        if (a==1 || (a==3 && chain.residues[r+1].resn=="GLY"))
                            continue;
                        v1[0]=allNCAC_arr[(r+1)*5+a][0];
                        v1[1]=allNCAC_arr[(r+1)*5+a][1];
                        v1[2]=allNCAC_arr[(r+1)*5+a][2];
                        ChangeCoor(v1, RotMatix, TranVect, v2);
                        if (a==4 && chain.residues[r+1].resn=="GLY")
                        {
                            chain.residues[r+1].atoms[3].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[r+1].atoms[3].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[r+1].atoms[3].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                        else if (a==5)
                        {
                            chain.residues[r+1].atoms.back().xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[r+1].atoms.back().xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[r+1].atoms.back().xyz[2]=1000*(v2[2]+v3[2]);
                        }
                        else
                        {
                            chain.residues[r+1].atoms[a].xyz[0]=1000*(v2[0]+v3[0]);
                            chain.residues[r+1].atoms[a].xyz[1]=1000*(v2[1]+v3[1]);
                            chain.residues[r+1].atoms[a].xyz[2]=1000*(v2[2]+v3[2]);
                        }
                    }
                }

                /* update current residue */
                v1[0]=target_stru[1][0];
                v1[1]=target_stru[1][1];
                v1[2]=target_stru[1][2];
                ChangeCoor(v1, RotMatix, TranVect, v2);
                TranVect[0]+=refrence_stru[1][0]-v2[0];
                TranVect[1]+=refrence_stru[1][1]-v2[1];
                TranVect[2]+=refrence_stru[1][2]-v2[2];
                for (a=0;a<5;a++)
                {
                    if (a==1 || (a==3 && chain.residues[r].resn=="GLY"))
                        continue;
                    v1[0]=allNCAC_arr[r*5+a][0];
                    v1[1]=allNCAC_arr[r*5+a][1];
                    v1[2]=allNCAC_arr[r*5+a][2];
                    ChangeCoor(v1, RotMatix, TranVect, v2);
                    if (a==4 && chain.residues[r].resn=="GLY")
                    {
                        chain.residues[r].atoms[3].xyz[0]=1000*v2[0];
                        chain.residues[r].atoms[3].xyz[1]=1000*v2[1];
                        chain.residues[r].atoms[3].xyz[2]=1000*v2[2];
                    }
                    else
                    {
                        chain.residues[r].atoms[a].xyz[0]=1000*v2[0];
                        chain.residues[r].atoms[a].xyz[1]=1000*v2[1];
                        chain.residues[r].atoms[a].xyz[2]=1000*v2[2];
                    }
                }
            }
            /* side chain atoms other than CB */
            for (r=0;r<L;r++)
            {
                resn=chain.residues[r].resn;
                if (resn=="GLY" || resn=="ALA") continue; //2

                /* chi-1 */
                if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                else                     fp.read((char *)&tor8,sizeof(int8_t));
                tor=tor8*180./INT8_MAX;
                xyz2double(chain.residues[r].atoms[0],c1); //N
                xyz2double(chain.residues[r].atoms[1],c2); //CA
                xyz2double(chain.residues[r].atoms[3],c3); //CB
                if (resn=="SER")
                    calculateCoordinates(c4,c1,c2,c3,1.417,110.773,tor);
                else if (resn=="CYS")
                    calculateCoordinates(c4,c1,c2,c3,1.808,113.8169,tor);
                else if (resn=="VAL" || resn=="ILE")
                    calculateCoordinates(c4,c1,c2,c3,1.527,110.7,tor);
                else if (resn=="LEU")
                    calculateCoordinates(c4,c1,c2,c3,1.53,116.10,tor);
                else if (resn=="THR")
                    calculateCoordinates(c4,c1,c2,c3,1.43,109.18,tor); //OG1
                else if (resn=="ARG" || resn=="LYS")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.83,tor);
                else if (resn=="ASP")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.06,tor);
                else if (resn=="ASN")
                    calculateCoordinates(c4,c1,c2,c3,1.52,112.62,tor);
                else if (resn=="GLU")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.82,tor);
                else if (resn=="GLN")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.75,tor);
                else if (resn=="MET")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.68,tor);
                else if (resn=="HIS")
                    calculateCoordinates(c4,c1,c2,c3,1.49,113.74,tor);
                else if (resn=="PRO")
                    calculateCoordinates(c4,c1,c2,c3,1.49,104.21,tor);
                else if (resn=="PHE")
                    calculateCoordinates(c4,c1,c2,c3,1.50,113.85,tor);
                else if (resn=="TYR")
                    calculateCoordinates(c4,c1,c2,c3,1.51,113.8,tor);
                else if (resn=="TRP")
                    calculateCoordinates(c4,c1,c2,c3,1.50,114.10,tor);
                if (resn=="THR")
                {
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //OG1
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //OG1
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //OG1
                    calculateCoordinates(c4,c1,c2,c3,1.53,111.13,tor-120);
                }
                chain.residues[r].atoms[5].xyz[0]=1000.*c4[0]; //CG
                chain.residues[r].atoms[5].xyz[1]=1000.*c4[1]; //CG
                chain.residues[r].atoms[5].xyz[2]=1000.*c4[2]; //CG
                if (resn=="ILE")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.527,110.7,tor-120);
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //CG2
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //CG2
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //CG2
                }
                else if (resn=="VAL")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.527,110.7,tor+120);
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //CG2
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //CG2
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //CG2
                    continue;
                }
                else if (resn=="CYS" || resn=="SER" || resn=="THR")
                    continue; // 2+4=6

                /* chi-2: CA-CB-CG-XD: ARG LYS MET GLU GLN ILE LEU
                 * HIS TYR TRP PHE PRO ASP ASN */
                if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                else                     fp.read((char *)&tor8,sizeof(int8_t));
                tor=tor8*180./INT8_MAX;
                xyz2double(chain.residues[r].atoms[1],c1); //CA
                xyz2double(chain.residues[r].atoms[3],c2); //CB
                xyz2double(chain.residues[r].atoms[5],c3); //CG
                if (resn=="ILE")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.97,tor);
                else if (resn=="LEU")
                    calculateCoordinates(c4,c1,c2,c3,1.524,110.27,tor);
                else if (resn=="ARG" || resn=="LYS")
                    calculateCoordinates(c4,c1,c2,c3,1.52,111.79,tor);
                else if (resn=="ASP")
                    calculateCoordinates(c4,c1,c2,c3,1.25,119.22,tor);//OD1
                else if (resn=="ASN")
                    calculateCoordinates(c4,c1,c2,c3,1.23,120.85,tor);//OD1
                else if (resn=="GLU")
                    calculateCoordinates(c4,c1,c2,c3,1.52,113.31,tor);
                else if (resn=="GLN")
                    calculateCoordinates(c4,c1,c2,c3,1.52,112.78,tor);
                else if (resn=="MET")
                    calculateCoordinates(c4,c1,c2,c3,1.81,112.69,tor);
                else if (resn=="HIS")
                    calculateCoordinates(c4,c1,c2,c3,1.35,130.61,tor); //CD2
                else if (resn=="PRO")
                    calculateCoordinates(c4,c1,c2,c3,1.50,105.03,tor);
                else if (resn=="PHE")
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,tor); //CD1
                else if (resn=="TYR")
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.98,tor); //CD1
                else if (resn=="TRP")
                    calculateCoordinates(c4,c1,c2,c3,1.37,127.07,tor); //CD1
                
                if (resn=="ILE")
                {
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD1
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD1
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD1
                }
                else
                {
                    chain.residues[r].atoms[6].xyz[0]=1000.*c4[0]; //XD
                    chain.residues[r].atoms[6].xyz[1]=1000.*c4[1]; //XD
                    chain.residues[r].atoms[6].xyz[2]=1000.*c4[2]; //XD
                }
                if (resn=="LEU")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.525,110.58,tor+120);//CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2
                    continue;
                }
                else if (resn=="ASP")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.25,118.218,tor+180);//OD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //OD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //OD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //OD2
                    continue;
                }
                else if (resn=="ASN")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.33,116.48,tor+180);//ND2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //ND2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //ND2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //ND2
                    continue;
                }
                else if (resn=="HIS")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.38,122.85,tor+180);//ND1
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //ND1
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //ND1
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //ND1
                    
                    // CB_CG_ND1_CE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[7],c3); //ND1
                    calculateCoordinates(c4,c1,c2,c3,1.32,108.5,180);//CE1
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE1
                    
                    // CB_CG_CD2_NE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[6],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.35,108.5,180);//NE2
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //NE2
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //NE2
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //NE2
                    continue;
                }
                else if (resn=="PHE")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,tor+180); //CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2

                    // CB_CG_CD1_CE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[6],c3); //CD1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE1
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE1

                    // CB_CG_CD2_CE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[7],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE2
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //CE2
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //CE2
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //CE2
                    
                    // CG_CD1_CE1_CZ_diangle = 0.0
                    xyz2double(chain.residues[r].atoms[5],c1); //CG
                    xyz2double(chain.residues[r].atoms[6],c2); //CD1
                    xyz2double(chain.residues[r].atoms[8],c3); //CE1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,0); //CZ
                    chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //CZ
                    chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //CZ
                    chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //CZ
                    continue;
                }
                else if (resn=="TYR")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.82,tor+180); //CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2

                    // CB_CG_CD1_CE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[6],c3); //CD1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE1
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE1

                    // CB_CG_CD2_CE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[7],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,180); //CE2
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //CE2
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //CE2
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //CE2
                    
                    // CG_CD1_CE1_CZ_diangle = 0.0
                    xyz2double(chain.residues[r].atoms[5],c1); //CG
                    xyz2double(chain.residues[r].atoms[6],c2); //CD1
                    xyz2double(chain.residues[r].atoms[8],c3); //CE1
                    calculateCoordinates(c4,c1,c2,c3,1.39,120.0,0); //CZ
                    chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //CZ
                    chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //CZ
                    chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //CZ

                    // CD1_CE1_CZ_OH_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[6],c1); //CD1
                    xyz2double(chain.residues[r].atoms[8],c2); //CE1
                    xyz2double(chain.residues[r].atoms[10],c3); //CZ
                    calculateCoordinates(c4,c1,c2,c3,1.39,119.78,180); //OH
                    chain.residues[r].atoms[11].xyz[0]=1000.*c4[0]; //OH
                    chain.residues[r].atoms[11].xyz[1]=1000.*c4[1]; //OH
                    chain.residues[r].atoms[11].xyz[2]=1000.*c4[2]; //OH
                    continue;
                }
                else if (resn=="TRP")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.43,126.66,tor+180); //CD2
                    chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //CD2
                    chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //CD2
                    chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //CD2

                    // CB_CG_CD1_NE1_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[3],c1); //CB
                    xyz2double(chain.residues[r].atoms[5],c2); //CG
                    xyz2double(chain.residues[r].atoms[6],c3); //CD1
                    calculateCoordinates(c4,c1,c2,c3,1.38,108.5,180); //NE1
                    chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //NE1
                    chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //NE1
                    chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //NE1

                    // CB_CG_CD2_CE2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[7],c3); //CD2
                    calculateCoordinates(c4,c1,c2,c3,1.40,108.5,180); //CE2
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //CE2
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //CE2
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //CE2

                    // CB_CG_CD2_CE3_diangle = 0.0
                    calculateCoordinates(c4,c1,c2,c3,1.40,133.83,0); //CE3
                    chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //CE3
                    chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //CE3
                    chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //CE3

                    // CG_CD2_CE2_CZ2_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[5],c1); //CG
                    xyz2double(chain.residues[r].atoms[7],c2); //CD2
                    xyz2double(chain.residues[r].atoms[8],c3); //CE2
                    calculateCoordinates(c4,c1,c2,c3,1.40,120.0,180); //CZ2
                    chain.residues[r].atoms[12].xyz[0]=1000.*c4[0]; //CZ2
                    chain.residues[r].atoms[12].xyz[1]=1000.*c4[1]; //CZ2
                    chain.residues[r].atoms[12].xyz[2]=1000.*c4[2]; //CZ2

                    // CG_CD2_CE3_CZ3_diangle = 180.0
                    xyz2double(chain.residues[r].atoms[9],c3); //CE3
                    calculateCoordinates(c4,c1,c2,c3,1.40,120.0,180); //CZ3
                    chain.residues[r].atoms[13].xyz[0]=1000.*c4[0]; //CZ3
                    chain.residues[r].atoms[13].xyz[1]=1000.*c4[1]; //CZ3
                    chain.residues[r].atoms[13].xyz[2]=1000.*c4[2]; //CZ3

                    // CD2_CE2_CZ2_CH2_diangle = 0.0
                    xyz2double(chain.residues[r].atoms[7],c1); //CD2
                    xyz2double(chain.residues[r].atoms[8],c2); //CE2
                    xyz2double(chain.residues[r].atoms[12],c3); //CZ2
                    calculateCoordinates(c4,c1,c2,c3,1.40,120.0,0); //CH2
                    chain.residues[r].atoms[11].xyz[0]=1000.*c4[0]; //CH2
                    chain.residues[r].atoms[11].xyz[1]=1000.*c4[1]; //CH2
                    chain.residues[r].atoms[11].xyz[2]=1000.*c4[2]; //CH2
                    continue;
                }
                else if (resn=="ILE" || resn=="PRO") continue; // 2+4+9=15

                /* chi-3: CB-CG-CD-XE: ARG LYS GLN GLU MET */
                if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                else                     fp.read((char *)&tor8,sizeof(int8_t));
                tor=tor8*180./INT8_MAX;
                xyz2double(chain.residues[r].atoms[3],c1); //CB
                xyz2double(chain.residues[r].atoms[5],c2); //CG
                xyz2double(chain.residues[r].atoms[6],c3); //CD
                if (resn=="ARG"||resn=="LYS")
                    calculateCoordinates(c4,c1,c2,c3,1.46,111.68,tor);
                else if (resn=="GLU")
                    calculateCoordinates(c4,c1,c2,c3,1.25,119.02,tor);
                else if (resn=="GLN")
                    calculateCoordinates(c4,c1,c2,c3,1.24,120.86,tor);
                else if (resn=="MET")
                    calculateCoordinates(c4,c1,c2,c3,1.79,100.61,tor);
                
                if (resn=="GLN")
                {
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //OE1
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //OE1
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //OE1
                    calculateCoordinates(c4,c1,c2,c3,1.33,116.50,tor+180);
                }
                chain.residues[r].atoms[7].xyz[0]=1000.*c4[0]; //XE
                chain.residues[r].atoms[7].xyz[1]=1000.*c4[1]; //XE
                chain.residues[r].atoms[7].xyz[2]=1000.*c4[2]; //XE
                if (resn=="GLU")
                {
                    calculateCoordinates(c4,c1,c2,c3,1.25,119.02,tor+180);
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //OE2
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //OE2
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //OE2
                    continue;
                }
                if (resn=="GLN" || resn=="MET") continue; // 2+4+9+3=18
                
                /* chi-4: CG-CD-XE-XZ: ARG LYS */
                if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                else                     fp.read((char *)&tor8,sizeof(int8_t));
                tor=tor8*180./INT8_MAX;
                xyz2double(chain.residues[r].atoms[5],c1); //CG
                xyz2double(chain.residues[r].atoms[6],c2); //CD
                xyz2double(chain.residues[r].atoms[7],c3); //XE
                calculateCoordinates(c4,c1,c2,c3,1.33,124.79,tor);
                if (resn=="LYS")
                {
                    chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //NZ
                    chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //NZ
                    chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //NZ
                    continue; // 2+4+9+3+1=19
                }

                chain.residues[r].atoms[10].xyz[0]=1000.*c4[0]; //CZ
                chain.residues[r].atoms[10].xyz[1]=1000.*c4[1]; //CZ
                chain.residues[r].atoms[10].xyz[2]=1000.*c4[2]; //CZ
                /* chi-5: CD-NE-CZ-NH1: ARG */
                if (use_stdin)          cin.read((char *)&tor8,sizeof(int8_t));
                else if (use_pstream) fp_gz.read((char *)&tor8,sizeof(int8_t));
                else                     fp.read((char *)&tor8,sizeof(int8_t));
                tor=tor8*180./INT8_MAX;
                xyz2double(chain.residues[r].atoms[6],c1); //CG
                xyz2double(chain.residues[r].atoms[7],c2); //NE
                xyz2double(chain.residues[r].atoms[10],c3); //CZ
                calculateCoordinates(c4,c1,c2,c3,1.33,120.64,tor);
                chain.residues[r].atoms[8].xyz[0]=1000.*c4[0]; //NH1
                chain.residues[r].atoms[8].xyz[1]=1000.*c4[1]; //NH1
                chain.residues[r].atoms[8].xyz[2]=1000.*c4[2]; //NH1
                calculateCoordinates(c4,c1,c2,c3,1.33,120.64,tor+180);
                chain.residues[r].atoms[9].xyz[0]=1000.*c4[0]; //NH2
                chain.residues[r].atoms[9].xyz[1]=1000.*c4[1]; //NH2
                chain.residues[r].atoms[9].xyz[2]=1000.*c4[2]; //NH2
            }
            resn.clear();
            DeleteArray(&allCA_arr,L);
            DeleteArray(&allNCAC_arr,L*5+1);
            DeleteArray(&target_stru,3);
            DeleteArray(&refrence_stru,3);
            if (lossy==0)
            {
                vector<int16_t> d16_vec;
                vector<int8_t> d8_vec;
                int16_t d16;
                int8_t d8;
                int i,j;
                size_t d16_count=0;
                if (use_stdin)        cin.read((char *)&d16_count,sizeof(size_t));
                else if (use_pstream) fp_gz.read((char *)&d16_count,sizeof(size_t));
                else                  fp.read((char *)&d16_count,sizeof(size_t));
                for (i=0;i<d16_count;i++)
                {
                    if (use_stdin)        cin.read((char *)&d16,sizeof(int16_t));
                    else if (use_pstream) fp_gz.read((char *)&d16,sizeof(int16_t));
                    else                  fp.read((char *)&d16,sizeof(int16_t));
                    d16_vec.push_back(d16);
                }
                j=0;
                for (r=0;r<L;r++)
                {
                    for (a=0;a<chain.residues[r].atoms.size();a++)
                    {
                        if (a==1) continue;
                        for (i=0;i<3;i++)
                        {
                            if (use_stdin)        cin.read((char *)&d8,sizeof(int16_t));
                            else if (use_pstream) fp_gz.read((char *)&d8,sizeof(int16_t));
                            else                  fp.read((char *)&d8,sizeof(int16_t));
                            if (d8==0)
                            { 
                                pep.chains[c].residues[r].atoms[a].xyz[i]+=d16_vec[j];
                                j++;
                            }
                            else
                                pep.chains[c].residues[r].atoms[a].xyz[i]+=d8;
                        }
                    }
                }
            }
            if (use_stdin)        cin.read((char *)&bfactor,sizeof(int16_t));
            else if (use_pstream) fp_gz.read((char *)&bfactor,sizeof(int16_t));
            else                  fp.read((char *)&bfactor,sizeof(int16_t));
            if (bfactor_mode==0)
            {
                for (r=0;r<L;r++)
                    for (a=0;a<chain.residues[r].atoms.size();a++)
                        chain.residues[r].atoms[a].bfactor=bfactor;
            }
            else
            {
                for (r=0;r<L;r++)
                {
                    if (use_stdin)        cin.read((char *)&db8,sizeof(int8_t));
                    else if (use_pstream) fp_gz.read((char *)&db8,sizeof(int8_t));
                    else                  fp.read((char *)&db8,sizeof(int8_t));
                    bfactor+=10*db8;
                    for (a=0;a<chain.residues[r].atoms.size();a++)
                        chain.residues[r].atoms[a].bfactor=bfactor;
                }
                if (bfactor_mode==2)
                {
                    bfactor=chain.residues[r].atoms[1].bfactor;
                    for (r=0;r<chain.residues.size();r++)
                    {
                        for (a=0;a<chain.residues[r].atoms.size();a++)
                        {
                            if (use_stdin)        cin.read((char *)&db8,sizeof(int8_t));
                            else if (use_pstream) fp_gz.read((char *)&db8,sizeof(int8_t));
                            else                  fp.read((char *)&db8,sizeof(int8_t));
                            chain.residues[r].atoms[a].bfactor=bfactor+10*db8;
                        }
                    }
                }
            }
        }
        else
        {
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
    return lossy;
}

string write_cif_structure(ModelUnit &pep,string &header)
{
    string txt="data_AF\n#\n_entry.id AF\n#\n";
    pdb2fasta(pep);
    stringstream buf;
    int a,r,c,s;
    vector<string> line_vec;
    Split(header,line_vec,'\n');
    string model_group_name="";
    string db_accession="";
    string db_code="";
    string db_name="";
    string ncbi_taxonomy_id="";
    string organism_scientific="";
    string seq_db_align_begin="";
    string seq_db_align_end="";
    string pdbx_description="";
    string asym_id="";
    string::size_type n;
    for (s=0;s<line_vec.size();s++)
    {
        if (StartsWith(line_vec[s],"DBREF "))
        {
            vector<string> dbref_line_vec;
            Split(line_vec[s],dbref_line_vec);
            db_name=dbref_line_vec[5];
            db_accession=dbref_line_vec[6];
            db_code=dbref_line_vec[7];
            seq_db_align_begin=dbref_line_vec[8];
            seq_db_align_end  =dbref_line_vec[9];
            r=atoi(seq_db_align_begin.c_str());
            r/=200;
            r++;
            buf<<"data_AF-"<<db_accession<<"-F"<<r
                <<"\n#\n_entry.id AF-"<<db_accession<<"-F"<<r<<"\n#\n";
            txt=buf.str();
            buf.str(string());
            vector<string>().swap(dbref_line_vec);
        }
        if (StartsWith(line_vec[s],"TITLE "))
        {
            n=line_vec[s].find(" FOR ");
            if (n!=string::npos)
                model_group_name=line_vec[s].substr(10,n-10);
        }
        if (StartsWith(line_vec[s],"COMPND"))
        {
            if (asym_id.size()==0)
            {
                n=line_vec[s].find(" CHAIN: ");
                if (n!=string::npos) asym_id=line_vec[s].substr(n+8);
            }
            if (pdbx_description.size()==0)
            {
                n=line_vec[s].find(" MOLECULE: ");
                if (n!=string::npos) pdbx_description=
                    rstrip(Trim(line_vec[s].substr(n+11)),";");
            }
        }
        if (StartsWith(line_vec[s],"SOURCE"))
        {
            if (organism_scientific.size()==0)
            {
                n=line_vec[s].find(" ORGANISM_SCIENTIFIC: ");
                if (n!=string::npos) organism_scientific=
                    rstrip(Trim(line_vec[s].substr(n+22)),";");
            }
            if (ncbi_taxonomy_id.size()==0)
            {
                n=line_vec[s].find(" ORGANISM_TAXID: ");
                if (n!=string::npos) ncbi_taxonomy_id=
                    rstrip(Trim(line_vec[s].substr(n+17)),";");
            }
        }
    }

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
                <<left<<setw(4)<<pep.chains[c].residues[r].resi<<endl;
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
"#\n";
    if (model_group_name.size()) txt+=""
"_ma_model_list.data_id          1\n"
"_ma_model_list.model_group_id   1\n"
"_ma_model_list.model_group_name \""+model_group_name+"\"\n"
"_ma_model_list.model_id         1\n"
"_ma_model_list.model_name       \"Top ranked model\"\n"
"_ma_model_list.model_type       \"Ab initio model\"\n"
"_ma_model_list.ordinal_id       1\n"
"#\n";

    txt+=""
"_ma_software_group.group_id    1\n"
"_ma_software_group.ordinal_id  1\n"
"_ma_software_group.software_id 1\n"
"#\n"
"_ma_target_entity.data_id   1\n"
"_ma_target_entity.entity_id 1\n"
"_ma_target_entity.origin    \"reference database\"\n"
"#\n";
    if (asym_id.size()) txt+=""
"_ma_target_entity_instance.asym_id   "+asym_id+"\n"
"_ma_target_entity_instance.details   .\n"
"_ma_target_entity_instance.entity_id 1\n"
"#\n";
    bool pdbx_sifts=(db_name.size() && db_accession.size() && seq_db_align_begin.size() && seq_db_align_end.size());
    if (pdbx_sifts && db_code.size() && ncbi_taxonomy_id.size() &&
        organism_scientific.size()) txt+=""
"_ma_target_ref_db_details.db_accession                 "+db_accession+"\n"
"_ma_target_ref_db_details.db_code                      "+db_code+"\n"
"_ma_target_ref_db_details.db_name                      "+db_name+"\n"
"_ma_target_ref_db_details.ncbi_taxonomy_id             "+ncbi_taxonomy_id+"\n"
"_ma_target_ref_db_details.organism_scientific          \""+organism_scientific+"\"\n"
"_ma_target_ref_db_details.seq_db_align_begin           "+seq_db_align_begin+"\n"
"_ma_target_ref_db_details.seq_db_align_end             "+seq_db_align_end+"\n"
"_ma_target_ref_db_details.target_entity_id             1\n"
"#\n";
    txt+=""
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
                <<left<<setw(4)<<pep.chains[c].residues[r].resi<<endl;
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
    if (pdbx_sifts) txt+=""
"_atom_site.pdbx_sifts_xref_db_acc\n"
"_atom_site.pdbx_sifts_xref_db_name\n"
"_atom_site.pdbx_sifts_xref_db_num\n"
"_atom_site.pdbx_sifts_xref_db_res\n";
    int pdbx_sifts_xref_db_num=atoi(seq_db_align_begin.c_str());
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
                    <<left<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)
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
                    <<Trim(pep.chains[c].residues[r].atoms[a].name)<<" 1 ";
                if (pdbx_sifts) buf<<db_accession<<' '<<db_name<<' '
                    <<setw(4)<<pdbx_sifts_xref_db_num<<' '
                    <<aa3to1(pep.chains[c].residues[r].resn)<<' ';
                buf<<endl;
                txt+=buf.str();
                buf.str(string());
            }
            pdbx_sifts_xref_db_num++;
        }
    }
    txt+="#\n";
    vector<string> ().swap(line_vec);
    return txt;
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_cif_structure(const string &filename,ModelUnit &pep,string &header)
{
    if (filename=="-")
        cout<<write_cif_structure(pep,header)<<flush;
    else
    {
        string filename_str=filename;
        if (EndsWith(filename,".gz")) filename_str=filename.substr(0,filename.size()-3);
        ofstream fp(filename_str.c_str());
        fp<<write_cif_structure(pep,header)<<flush;
        fp.close();
        if (EndsWith(filename,".gz"))
            int r=system(((string)("gzip -f "+filename_str)).c_str());
        filename_str.clear();
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
