// function: insert snps into the fasta file
// written by Hequan Sun, Max Planck Institude for Plant Breeding Research.
#include <iostream>
#include  <fstream>
#include  <sstream>
#include      <map>
#include   <string>
#include <stdlib.h>
#include "split_string.h"
using namespace std;
//
bool read_marker(string markerfile, map<string, char>* markerSet);
//
bool verbose = false;
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "\nFunction: insert given mutations to fasta";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\nUsage: insertMarkerToFasta seq.fasta markerfile.txt outfile_prefix" << endl;
        cout << "\nPls keep order of input." << endl;
        cout << "\nOutput will be outfile_prefix.fasta." << endl << endl;
        return 1;
    }
    cout << "\nInput DNA sequence is from file: " << argv[1] << "." << endl;
    //
    string markerfile = (string)argv[2];
    map<string, char> markerSet;
    if(!read_marker(markerfile, &markerSet))
    {
        cout << "   Error: reading marker set failed. " << endl;
        return 1;
    }
    // write markers (inserted in fasta; some may not)
    string omkrfile = "inserted_" + markerfile; // in case there are repeated info in the original file.
    ofstream mkrfp;
    mkrfp.open(omkrfile.c_str(), ios::out);
    if(!mkrfp.good())
    {
        cout << "   Error: cannot open output file." << endl;
        return false; 
    }    
    // open seq file
    std::ifstream fp(argv[1]);
    if(!fp.is_open())
    {
        cout << "Cannot open dna seq file " << argv[1] << ". Exited." << endl;
        return 1;
    }
    // initialize a file where mutated sequences will be recorded.
    fstream outfp;
    std::stringstream ss("");
    ss << argv[3] << ".fasta"; 
    outfp.open((ss.str()).c_str(), ios::out);
    if(!outfp.is_open())
    {
        cout << "Cannot open file " << ss.str() << " to write mutated dna info. Exited.\n";
        return 1;
    }
    cout << "\nMutated DNA sequence will be collected in file: " << ss.str() << "." << endl;
    // read sequence and mutate
    std::string line("");
    getline(fp, line);
    while(fp.good())
    {
        if(line.size()==0 || line.find("#")!=std::string::npos)
        {
            getline(fp, line);
            continue;
        }        
        int num = 0;        
        if(line[0]=='>')
        {
            string seq("");
            string chrseq("");   
            while(fp.good())
            {             
                getline(fp, chrseq);
                if(chrseq.find(">")!=std::string::npos) break; // next seq
                seq    += chrseq;
            }            
            for(size_t i=0; i<seq.length(); i++)
            {
                std::stringstream key;
                key.str("");
                key << line.substr(1) << "#" << i+1;                
                if(markerSet.find(key.str()) != markerSet.end())
                {
                    mkrfp << "ins\t" << line.substr(1) << "\t" << i+1 << "\t" << seq[i] << "\t" << markerSet[key.str()] << endl;                     
                    seq[i] = markerSet[key.str()];
                    num ++;
                }
                //cout << seq[oripos] << endl;
            }
            outfp << line << "\n";
            for(int pi = 0; pi < seq.size(); pi ++)
            {
                outfp << seq.substr(pi, 1);
                if((pi+1)%80==0 || pi==seq.size()-1)
                {
                    outfp << "\n";
                }
            }  
            cout << "Length of original sequence " << line.substr(1) << " is " << seq.length() << ", of which "
                 << num << " have been point-mutated with given markers."    << endl;
            line.clear();
            line += chrseq;
        }
        else
        {
            getline(fp, line);
        }
    }    
    fp.close();
    outfp.close();
    mkrfp.close();
}
//
bool read_marker(string markerfile, map<string, char>* markerSet)
{
    ifstream fp;
    fp.open(markerfile.c_str(), ios::in);
    if(!fp.good())
    {
        cout << "   Error: cannot open marker file " << markerfile << endl;
        return false;
    }
    while(fp.good())
    {
        string line("");
        getline(fp, line);
        if(line.size()==0 || line[0]=='#') 
        {
            line.clear(); 
            continue;
        }
        vector<string> lineinfo = split_string(line, '\t'); // simP1	Chr1	33	C	T        
        string key = lineinfo[1] +"#"+ lineinfo[2]; // note: here in file chr.pos is 1 based.
        char   val = (lineinfo[4].c_str())[0];        
        if((*markerSet).find(key) == (*markerSet).end())
        {
            (*markerSet).insert(std::pair<string, char>(key, val));
        }
        else
        {
            cout << "   Warning: repeated chr pos info: " << key << endl;
        }
    }
    fp.close();
    return true;
}





















