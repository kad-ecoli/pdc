const char* docstring=""
"pdc input.pdb output.pdc [options]\n"
"    convert PDB format input to PDC format output\n"
"\n"
"option:\n"
"    -a={2,1,0} atomic details to preserve\n"
"               0: CA atoms\n"
"               1: N, CA, C atoms\n"
"               2: (default) all atoms\n"
"    -l={0,1}   lossless/lossy compression\n"
"               0: (default) lossless compression\n"
"               1: loosy compression\n"
;

#include "PDBParser.hpp"

int main(int argc,char **argv)
{
    string infile ="";
    string outfile="";
    int atomic_detail=2;
    int loosy        =0;

    for (int a=1;a<argc;a++)
    {
        if (StartsWith(argv[a],"-a="))
            atomic_detail=atoi(((string)(argv[a])).substr(3).c_str());
        else if (StartsWith(argv[a],"-l="))
            loosy=atoi(((string)(argv[a])).substr(3).c_str());
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
    
    string header;
    ModelUnit pdb_entry=read_pdb_structure(infile.c_str(),header,atomic_detail);

    map<string, map<string,int> >ordMap;
    initialize_atom_order_map(ordMap);
    standardize_pdb_order(pdb_entry, ordMap);
    //write_pdb_structure(outfile.c_str(),pdb_entry,header);
    write_pdc_structure(outfile.c_str(),pdb_entry,header);

    /* clean up */
    string ().swap(infile);
    string ().swap(outfile);
    string ().swap(header);
    deepClean(pdb_entry);
    map<string, map<string,int> >().swap(ordMap);
    return 0;
}
