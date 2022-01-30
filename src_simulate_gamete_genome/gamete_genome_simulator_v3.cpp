/* this function, given 

	inputs:
		1. haplotype genome A
		2. haplotype genome B
		3. crossover landscape along chromosomes	
			
	Output haploid gamete genomes with 0~2 crossovers.
	
   Note, we simulate 0 co:29.3%, 1 co:58.9%, 2 co:11.8% => 0.83 CO per chr (following Beth Rowan et al 2019)
	
   Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
*/
#include       <map>
#include    <string>
#include    <vector>
#include   <fstream>
#include   <sstream>
#include  <iostream>
#include <algorithm>
#include  <string.h>
#include  <stdlib.h>
#include    <math.h>
#include "split_string.h"
using namespace std;
struct CONODE{
    unsigned long arr[10000]; // co breakpoint array
    int           freq[10000];// co breakpoint frequency
    int            cnt;
};
int findCeil(unsigned long arr[], int r, int l, int h);
int myRand(unsigned long arr[], int freq[], int n);
bool read_co_landscape(string fileCOLand, map<string, CONODE>* comap);
// 
int main(int argc, char* argv[])
{
    if(argc < 7)
    {
        cout << "\nFunction: given two haplotype genomes, generate a number of haploid gamete genomes, each haploid / fasta.\n";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\nUsage: gamete_genome_simulator genome_A.fa genome_B.fa co_landscape.txt N_recombinant output_prefix random_seed pure_CO"  << endl;
        cout << "         genome_A.fa.....: genome sequence of parent A"                           << endl;
        cout << "         genome_B.fa.....: genome sequence of parent B"                           << endl;       
        cout << "         co_landscape.txt: frequency count of crossovers within fixed windows"    << endl;
        cout << "         N_recombinant...: number of expected recombinants in the pool [50]"      << endl;
        cout << "         rand_seed.......: for randomizing numbers/location of COs"               << endl;
        cout << "         pure_CO.........: 0=only simulate CO; otherwise, simulate CO and gametes"<< endl;
        cout << "         NOTE: order of chrs in genome_A.fa and genome_B.fa must be the same. "   << endl;
        return 1;
    }
    // set parameters
    string fileGeno1  = (string)argv[1];
    string fileGeno2  = (string)argv[2];
    string fileCOLand = (string)argv[3];
    int    N          = atoi(argv[4]); 
    map<unsigned long, int> existingBP;
    string output_prefix = (string)argv[5];  
    double iseed         = atof(argv[6])*1000/2; // for randomization    
    bool simulate_genome = (bool)atoi(argv[7]);// 0: only CO
    //
    double CO_0_ratio = 0.293; // no CO
    double CO_1_ratio = 0.589; // 1  CO
    double CO_2_ratio = 0.118; // 2  CO
    // step 1 simulate crossover breakpoints
    // int* arr;  // these are candidate break points
    // int* freq; // these are counts of break points
    // step 1.1 read co landscapes
    map<string, CONODE> comap;
    if(!read_co_landscape(fileCOLand, &comap))
    {
        return false;
    }
    // step 1.2 create recombination breakpoints
    map<string, map<int, vector<unsigned long> > > sim_bp; // <chr, <ri/gametei, vector<position> > > // 0~2 COs
    //
    map<string, CONODE>::iterator chritr;
    map<string, CONODE>::iterator chritr_end;
    chritr     = comap.begin();
    chritr_end = comap.end();
    while(chritr != chritr_end)
    {
        string chr   = (*chritr).first;
        CONODE tmpco = (*chritr).second;
        unsigned long arr[tmpco.cnt];
        int           freq[tmpco.cnt];
        int        n = tmpco.cnt;
        //
        if(sim_bp.find(chr) == sim_bp.end() )
        {
            map<int, vector<unsigned long> > tmp_bplist;
            sim_bp.insert(std::pair<string, map<int, vector<unsigned long> > >(chr, tmp_bplist) );
        }
        //
        for(int i = 0; i < tmpco.cnt; i ++)
        {
            arr[i]  = tmpco.arr[i];
            freq[i] = tmpco.freq[i];
        }
        // Use a different seed value for every run.     
        // srand (time( NULL ));
        srand(iseed);
        // Let us generate N=500 breakpoints accroding to given distribution: win_sta + ( std::rand() % ( win_end - win_sta + 1 ) ) => random numbers in this window        
        int CO_cnt[] = {0, 0, 0};
        double scaling_down_CO_0 = 0.8;        
        for (int i = 0; i < N; i++)
        {
            //
            vector<unsigned long> tmp_co_pos; // for this chr of this gamete i
            sim_bp[chr].insert(std::pair<int, vector<unsigned long> >(i, tmp_co_pos));
            // number of CO 
            int nco = 0; // which recombined chr: 0~2 crossovers
            double randra = rand()%1000000000*1.0 / 1000000000.0;
            if(randra<=CO_0_ratio*scaling_down_CO_0) // [0,0.293]
            {
                nco = 0;
                CO_cnt[0] ++;
            }else
            if(randra>CO_0_ratio*scaling_down_CO_0 && randra<=CO_0_ratio+CO_1_ratio) // (0.293, 0.293+0.589]
            {
                nco = 1;
                CO_cnt[1] ++;                
            }else
            if(randra>CO_0_ratio+CO_1_ratio && randra<=1.0)        // (0.293+0.589, 1]
            {
                nco = 2;
                CO_cnt[2] ++;                
            }else
            {
                cout << "   Warning: you should not reach here!";
            }  
            // cout << "   check: nco = " << nco << endl;
            if(nco > 0)
            {
                for(int ni = 0; ni < nco; ni ++)
                {
                    unsigned long this_bp  = myRand(arr, freq, n);
                    if(ni > 0)
                    {
                        unsigned long last_bp = sim_bp[chr][i][ni-1];
                        while(std::find(sim_bp[chr][i].begin(), sim_bp[chr][i].end(), this_bp) != sim_bp[chr][i].end() ||
                              abs(this_bp-last_bp)<2000000)
                        {
                            this_bp  = myRand(arr, freq, n);
                        } 
                    }                
                    unsigned long this_sta = this_bp - 50000 + 1;
                    unsigned long this_end = this_bp + 50000;
                    this_bp                = this_sta + ( std::rand() % ( this_end - this_sta + 1 ) ); // not exactly mid anymore
                    cout << "bp_breakpoint gamete " << i << "\t" << chr << "\t" << this_sta << "\t" << this_end << "\t" << this_bp << endl;
                    // collect
                    sim_bp[chr][i].push_back(this_bp);
                }
            }else
            {
                // no co
                unsigned long this_bp  = 0; // no CO
                long this_sta = -1;
                long this_end = -1;                
                cout << "bp_breakpoint_noCO gamete " << i << "\t" << chr << "\t" << this_sta << "\t" << this_end << "\t" << 0 << endl;                
                // collect
                sim_bp[chr][i].push_back(this_bp);                
            }
        }
        //
        cout << "   Targeted ratio: 0 CO: "<< CO_0_ratio*scaling_down_CO_0 
             << ", 1 CO: "                 << CO_0_ratio+CO_1_ratio - CO_0_ratio*scaling_down_CO_0
             << ", 2 COs: "                << 1-(CO_0_ratio+CO_1_ratio) << endl; 
        cout << "   Info: Simulated at "   << chr       << ", " 
                                           << CO_cnt[0] << " 0 CO, " 
                                           << CO_cnt[1] << " 1 CO, "
                                           << CO_cnt[2] << " 2 COs, in "
                                           << N         << " gamete genomes. "
                                           << endl;
        //
        chritr ++;
    }
    //
    if(!simulate_genome)
    {
        cout << "   Info: only simulating CO breakpoints asked. Stopping now. " << endl;
        return 0;
    }
    // open files
    std::ifstream fp1(fileGeno1.c_str());
    if(!fp1.is_open())
    {
        cout << "Cannot open genome 1 file " << fileGeno1 << ". Exited." << endl;
        return 1;
    }
    std::ifstream fp2(fileGeno2.c_str());
    if(!fp2.is_open())
    {
        cout << "Cannot open genome 2 file " << fileGeno2 << ". Exited." << endl;
        return 1;
    }
    //
    // srand (time( NULL ));
    srand (iseed);
    // step 2 parse chromosome files and generate gamete genomes
    std::string  line1("");
    getline(fp1, line1);
    std::string  line2("");
    getline(fp2, line2); 
    int ichr = 0;
    while(fp1.good() && fp2.good())
    {
        if(line1.size()==0 || line1.find("#")!=std::string::npos)
        {
            getline(fp1, line1);
            continue;
        }
        if(line2.size()==0 || line2.find("#")!=std::string::npos)
        {
            getline(fp2, line2);
            continue;
        }
        // get chr-i from genome 1
        string seq1("");   
        string chrseq1("");   
        string chr1name("");     
        if(line1[0]=='>')
        {
            chr1name = line1.substr(1);
            while(fp1.good())
            {             
                getline(fp1, chrseq1);
                if(chrseq1.find(">")!=std::string::npos) break; // next seq
                seq1    += chrseq1;
            }
            ichr ++;
        }
        else
        {
            getline(fp1, chrseq1);
        }
        // get chr-i from genome 2        
        string seq2("");  
        string chrseq2(""); 
        string chr2name("");                     
        if(line2[0]=='>')
        {
            chr2name = line2.substr(1);        
            while(fp2.good())
            {             
                getline(fp2, chrseq2);
                if(chrseq2.find(">")!=std::string::npos) break; // next seq
                seq2    += chrseq2;
            }
        }
        else
        {
            getline(fp2, chrseq2);
        }        
        // check consistency in sequence ids
        if(chr1name.compare(chr2name) != 0)
        {
            cout << "   Error: inconsistent sequence ID/order in genome files. Please correct! " << endl;
            return 1;
        }else
        {
            cout << "   Info: simulate gamete sequence with recombination from sequences: " << chr1name
                 << ", sequence 1 size " << seq1.size() << " bp"
                 << ", sequence 2 size " << seq2.size() << " bp." << endl;
        }      
        // check if chr in given co landscape map
        if(sim_bp.find(chr1name) == sim_bp.end())
        {
            cout << "   Error: cannot find chr in landscape list. " << endl;
            return 1;
        }
        // recombine sequence
        for(int ri = 0; ri < N; ri ++)
        {
            // open a file where a haploid genome will be recorded.
            fstream ofp;
            std::stringstream ss;
            ss.str("");
            ss << output_prefix << "_recombinant_gamete_" << ri << ".fa";
            if(ichr==1)
            {
                ofp.open((ss.str()).c_str(), ios::out);                 
            }
            else
            {
                ofp.open((ss.str()).c_str(), ios::out | std::fstream::app);
            }
            if(!ofp.is_open())
            {
                cout << "Cannot open file " << ss.str() << " to write pooled recombinants. Exited.\n";
                return 1;
            }
            //
            vector<unsigned long> all_gamete_co_pos = sim_bp[chr1name][ri];
            if(all_gamete_co_pos.size() == 1 && sim_bp[chr1name][ri][0]==0)
            {
                unsigned long pos  = sim_bp[chr1name][ri][0];
                cout << "   check: no recombination at " << chr1name << endl;
                int nchr = rand()%2+1; // which haplotype chr to take
                if(nchr == 1)
                {
                    ofp << ">" << line1.substr(1) << "." << ri << ".g1" << " non-recomb-from.g1-" << ri << " bp=" << pos << endl;
                    std::stringstream ss; 
                    ss.str("");
                    ss << seq1;
                    string iseq = ss.str();
                    for(int pi = 0; pi < iseq.size(); pi ++)
                    {    
                        ofp << iseq.substr(pi, 1);
                        if((pi+1)%80==0 || pi==iseq.size()-1)
                        {
                            ofp << "\n";
                        }
                    }
                }
                else
                {
                    ofp << ">" << line2.substr(1) << "." << ri << ".g2" << " non-recomb-from.g2-" << ri << " bp=" << pos << endl;
                    ss.str("");
                    ss << seq2;     
                    string iseq = ss.str();
                    for(int pi = 0; pi < iseq.size(); pi ++)
                    {
                        ofp << iseq.substr(pi, 1);
                        if((pi+1)%80==0 || pi==iseq.size()-1)
                        {
                            ofp << "\n";
                        }
                    }
                }                 
            }else            
            if(all_gamete_co_pos.size() == 1 && sim_bp[chr1name][ri][0]!=0)
            {
                // case 1: 1 CO
                unsigned long pos  = sim_bp[chr1name][ri][0];
                cout << "   check: to be recombined at " << chr1name << "\t" << pos << endl;
                int nchr = rand()%2+1; // which recombined chr to start with
                if(nchr == 1)
                {
                    ofp << ">" << line1.substr(1) << "." << ri << ".g1" << " recomb-from.g1-" << ri << " bp=" << pos << endl;
                    std::stringstream ss; 
                    ss.str("");
                    ss << seq1.substr(0, pos) << seq2.substr(pos);
                    string iseq = ss.str();
                    for(int pi = 0; pi < iseq.size(); pi ++)
                    {    
                        ofp << iseq.substr(pi, 1);
                        if((pi+1)%80==0 || pi==iseq.size()-1)
                        {
                            ofp << "\n";
                        }
                    }
                }
                else
                {
                    ofp << ">" << line2.substr(1) << "." << ri << ".g2" << " recomb-from.g2-" << ri << " bp=" << pos << endl;
                    ss.str("");
                    ss << seq2.substr(0, pos) << seq1.substr(pos);     
                    string iseq = ss.str();
                    for(int pi = 0; pi < iseq.size(); pi ++)
                    {
                        ofp << iseq.substr(pi, 1);
                        if((pi+1)%80==0 || pi==iseq.size()-1)
                        {
                            ofp << "\n";
                        }
                    }
                }  
            }else
            if(all_gamete_co_pos.size() == 2)
            {
                // case 2: 2 COs
                unsigned long pos1 = sim_bp[chr1name][ri][0];
                unsigned long pos2 = sim_bp[chr1name][ri][1];
                if(pos1 > pos2)
                {
                    pos1 = sim_bp[chr1name][ri][1];
                    pos2 = sim_bp[chr1name][ri][0];
                }
                int nchr = rand()%2+1; // which recombined chr to start with
                if(nchr == 1)
                {
                    ofp << ">" << line1.substr(1) << "." << ri << ".g1" << " recomb-from.g1-" << ri << " bp1=" << pos1 << " bp2=" << pos2 << "\n";
                    std::stringstream ss; 
                    ss.str("");
                    ss << seq1.substr(0, pos1) << seq2.substr(pos1, pos2-pos1) <<  seq1.substr(pos2);
                    string iseq = ss.str();
                    for(int pi = 0; pi < iseq.size(); pi ++)
                    {
                        ofp << iseq.substr(pi, 1);
                        if((pi+1)%80==0 || pi==iseq.size()-1)
                        {
                            ofp << "\n";
                        }
                    } 
                }
                else
                {                   
                    ofp << ">" << line2.substr(1) << "." << ri << ".g2" << " recomb-from.g2-" << ri << " bp1=" << pos1 << " bp2=" << pos2 << "\n"; 
                    ss.str("");
                    ss << seq2.substr(0, pos1) << seq1.substr(pos1, pos2-pos1) <<  seq2.substr(pos2);  
                    string iseq = ss.str();
                    for(int pi = 0; pi < iseq.size(); pi ++)
                    {
                        ofp << iseq.substr(pi, 1);
                        if((pi+1)%80==0 || pi==iseq.size()-1)
                        {
                            ofp << "\n";
                        }
                    } 
                }                  
            }else 
            {
                ; // not considered yet!
            }
        }                               
        // next seq
        line1.clear();
        line1 += chrseq1;
        line2.clear();
        line2 += chrseq2;
    } 
    // close files
    fp1.close();
    fp2.close();
    return 0;
}
// https://gist.github.com/vitaluha/8188337
// find ceiling of r in arr[l..h]
int findCeil(unsigned long arr[], int r, int l, int h)
{
    unsigned long mid;
    while (l < h)
    {
         mid = l + ((h - l) >> 1);  // Same as mid = (l+h)/2
        (r > arr[mid]) ? (l = mid + 1) : (h = mid);
    }
    return (arr[l] >= r) ? l : -1;
}
// find a random number from arr[] according to distribution array defined by freq[]. n is size of arrays.
int myRand(unsigned long arr[], int freq[], int n)
{
    // Create and fill prefix array
    unsigned long prefix[n], i;
    prefix[0] = freq[0];
    for (i = 1; i < n; ++i)
        prefix[i] = prefix[i - 1] + freq[i]; 
    // prefix[n-1] is sum of all frequencies. Generate a random number with value from 1 to this sum
    unsigned long r = (rand() % prefix[n - 1]) + 1; 
    // Find index of ceiling of r in prefix arrat
    int indexc = findCeil(prefix, r, 0, n - 1);
    return arr[indexc];
}
//
bool read_co_landscape(string fileCOLand, map<string, CONODE>* comap)
{
    ifstream ifp;
    ifp.open(fileCOLand.c_str(), ios::in);
    if(!ifp.good())
    {
        return false;
    }else
    {
        cout << "   Info: read CO landscape from file " << fileCOLand << endl;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        // new_Chr1	100001	200000	0.385	16
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<5)
        {
            cout << "   warning: 5 columns of chr co-sta co-end freq count needed, but insufficient at: " 
                 << line << endl;
            continue;
        }
        string        chr  = lineinfo[0];
        unsigned long sta  = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long end  = strtoul(lineinfo[2].c_str(), NULL, 0);
        unsigned long mid  = (sta+end-1)/2;
        int           freq = atoi(lineinfo[4].c_str());
        if( (*comap).find(chr) == (*comap).end() )
        {
            cout << "         : chr "  << chr    << " not found. " << endl;
            CONODE tmpco;
            tmpco.arr[0]  = mid;
            tmpco.freq[0] = freq;
            tmpco.cnt     = 1;
            (*comap).insert(std::pair<string, CONODE>(chr, tmpco));
        }else
        {
            int co_cnt    = (*comap)[chr].cnt;
            (*comap)[chr].arr[co_cnt] = mid;
            (*comap)[chr].freq[co_cnt] = freq;
            (*comap)[chr].cnt ++;
        }
    }
    ifp.close();
    // check
    bool check=false;
    if(check)
    {
        map<string, CONODE>::iterator chritr;
        map<string, CONODE>::iterator chritr_end;
        chritr     = (*comap).begin();
        chritr_end = (*comap).end();
        while(chritr != chritr_end)
        {
            string chr   = (*chritr).first;
            CONODE tmpco = (*chritr).second;
            //
            for(int i = 0; i < tmpco.cnt; i ++)
            {
                cout << "   check: " << chr << "\t" << tmpco.arr[i] << "\t" << tmpco.freq[i] << endl;
            }
            //
            chritr ++;
        }
    }
    //
    return true;
}
























