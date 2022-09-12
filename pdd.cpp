const char* docstring=""
"pdd input.pdc output.pdb [options]\n"
"    convert PDC format input to PDB format output\n"
;

#include "PDBParser.hpp"

int main(int argc,char **argv)
{
    string infile ="";
    string outfile="";

    for (int a=1;a<argc;a++)
    {
        if (infile.size()==0)
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
    
    map<string, vector<string> > ordMapR;
    initialize_reverse_atom_order_map(ordMapR);
    string header;
    ModelUnit pdb_entry=read_pdc_structure(infile.c_str(),header,ordMapR);
    write_pdb_structure(outfile.c_str(),pdb_entry,header);

    /* clean up */
    string ().swap(infile);
    string ().swap(outfile);
    string ().swap(header);
    deepClean(pdb_entry);
    map<string, vector<string> >().swap(ordMapR);
    return 0;
}
