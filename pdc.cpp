const char* docstring=""
"pdc input.pdb output.pdc [options]\n"
"    convert PDB format input to PDC format output\n"
"\n"
"option:\n"
"    -a={2,1,0} atomic details to preserve\n"
"               0: CA atoms\n"
"               1: N, CA, C atoms\n"
"               2: (default) all atoms\n"
"    -l={0,1,2,3,4} lossless/lossy compression\n"
"               0: (default) lossless compression\n"
"               1: lossless compression of CA; lossy compression of other atoms\n"
"               2: lossy compression of all atoms\n"
"               3: lossless compression, CA only\n"
"               4: lossy compression, CA only\n"
"    -f={0,1,2} input format\n"
"               0: (default) determined by filename\n"
"               1: PDB\n"
"               2: mmCIF/PDBx\n"
;

#include "PDBParser.hpp"

int main(int argc,char **argv)
{
    string infile ="";
    string outfile="";
    int atomic_detail=2;
    int lossy        =0;
    int infmt        =0;

    for (int a=1;a<argc;a++)
    {
        if (StartsWith(argv[a],"-f="))
            infmt=atoi(((string)(argv[a])).substr(3).c_str());
        else if (StartsWith(argv[a],"-a="))
            atomic_detail=atoi(((string)(argv[a])).substr(3).c_str());
        else if (StartsWith(argv[a],"-l="))
            lossy=atoi(((string)(argv[a])).substr(3).c_str());
        else if (infile.size()==0)
            infile=argv[a];
        else if (outfile.size()==0)
            outfile=argv[a];
        else
        {
            cerr<<"ERROR: unknown option "<<argv[a]<<endl;
            return 1;
        }
    }

    if (outfile.size()==0)
    {
        cerr<<docstring;
        return 1;
    }
    
    if (infmt==0)
    {
        if      (EndsWith(infile,".cif") || EndsWith(infile,".cif.gz")) infmt=2;
        else if (EndsWith(infile,".pdb") || EndsWith(infile,".pdb.gz") ||
                 EndsWith(infile,".ent") || EndsWith(infile,".ent.gz")) infmt=1;
        else 
        {
            cerr<<"WARNING! input PDB because format cannot be determined by output filename.\n"
                <<"use -f=1 or -f=2 to specify an output format"<<endl;
            infmt=1;
        }
    }

    string header;
    ModelUnit pdb_entry;
    if (infmt==2) pdb_entry=read_cif_structure(infile.c_str(),header,atomic_detail);
    else          pdb_entry=read_pdb_structure(infile.c_str(),header,atomic_detail);

    if (lossy<=2)
    {
        map<string, map<string,int> >ordMap;
        initialize_atom_order_map(ordMap);
        remove_nonter_OXT(pdb_entry);
        standardize_pdb_order(pdb_entry, ordMap);
        map<string, map<string,int> >().swap(ordMap);
    }
    else // CA only
    {
        standardize_pdb_ca(pdb_entry);
    }
    if (lossy==0) write_pdc_structure(outfile,pdb_entry,header);
    else write_pdc_lossy_structure(outfile,pdb_entry,header,lossy);

    /* clean up */
    string ().swap(infile);
    string ().swap(outfile);
    string ().swap(header);
    deepClean(pdb_entry);
    return 0;
}
