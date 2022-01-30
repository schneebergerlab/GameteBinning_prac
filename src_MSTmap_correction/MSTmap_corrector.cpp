// function: correct MSTmap-generated linkage group 
// written by Hequan Sun, Max Planck Institude for Plant Breeding Research.
#include  <iostream>
#include   <fstream>
#include   <sstream>
#include       <map>
#include    <string>
#include <algorithm>
#include  <stdlib.h>
#include  <assert.h>
#include "split_string.h"
using namespace std;
//
bool verbose = false;
bool read_gb_genotype(string gtfile, map<string, string>* marker_gt);
double calculate_recomb_freq(string PMpat1, string PMpat2);
bool phase_other_contig(map<int, map<string, int> >  known_phase, 
                        map<string, int>*            known_phase_dict, 
                        map<string, string>          marker_gt, 
                        bool                         output_phase,
                        double                       max_recomb_freq,                        
                        map<int, map<string, int> >* known_phase_updated);
//
int main(int argc, char* argv[])
{
    if(argc < 5)
    {
        cout << "\nFunction: correct MSTmap-generated linkage group ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")"                             << endl;
        cout << "\nUsage: MSTmap_corrector MSTmap_lg.txt TOP_X_LG s2_genotype_contig_seq.txt MAX_RECOMB_FREQ" 
             << endl << endl;
        cout << "   MSTmap_lg.txt is file generatd by MSTmap. "                   << endl;
        cout << "   TOP_X_LG      is the top X LGs with most contigs to output."  << endl;
        cout << "   s2_genotype_contig_seq.txt is from gamete_binning_dip. "      << endl;
        cout << "   MAX_RECOMB_FREQ is maximum recombined individuals when selecting markers to output in map." 
             << endl << endl;
        return 1;
    }
    //
    ifstream ifp;
    ifp.open(argv[1], ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot find MSTmap-generated linkage group file " << argv[1] << endl;
        return 1;
    }
    int TOP_X_LG  = atoi(argv[2]); // the number of lgs to selected
    string gtfile = (string)argv[3];
    map<string, string> marker_gt;
    if(!read_gb_genotype(gtfile, &marker_gt))
    {
        cout << "   Error: cannot read marker genotypes from file " << gtfile << endl;
        return 1;
    }
    double max_recomb = atof(argv[4]); // allow a maximum of x% recombined individuals    
    //
    int group_num = 0;
    map<int, map<string, double> > gmContig_all;      // <group_id, <marker_id, distance> >
    map<int, vector<string> >      gmContigOrder_all; // <group_id, <marker_id> >
    map<string, int>               marker_pairing;    // <marker_id, 2>
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);                
        if(line.find(";BEGINOFGROUP") != std::string::npos)
        {
            group_num ++;
            map<string, double> gmContig;
            vector<string>      gmContigOrder;            
            cout << "   Info: checking MSTmap linkage group " << group_num << " ... "  << endl;
            //
            getline(ifp, line);
            while(line.find(";ENDOFGROUP")==std::string::npos)
            {
                if(line.find("#") == std::string::npos)
                {
                    // this contig end marker
                    vector<string> lineinfo = split_string(line, '\t');
                    //
                    string marker_id      = lineinfo[0];                             // e.g., new_Chr1_9427397_10008496_le
                    string contig_id      = marker_id.substr(0, marker_id.size()-3); // e.g., new_Chr1_9427397_10008496
                    string leri           = marker_id.substr(marker_id.size()-2);
                    string marker_id_pair = "";
                    if(leri.compare("le")==0)
                    {
                        marker_id_pair    = contig_id + "_ri"; 
                    }else
                    {
                        marker_id_pair    = contig_id + "_le";
                    }
                    double distance  = atof(lineinfo[1].c_str());
                    vector<string>::iterator preitr = std::find(gmContigOrder.begin(), gmContigOrder.end(), marker_id_pair);
                    if(preitr != gmContigOrder.end())
                    {
                        preitr ++;
                        gmContigOrder.insert(preitr, marker_id);
                        // paired ri and le contig-end markers
                        marker_pairing.insert(std::pair<string, int>(marker_id,      2));
                        marker_pairing.insert(std::pair<string, int>(marker_id_pair, 2));
                    }else
                    {
                        gmContigOrder.push_back(marker_id);
                    }
                    gmContig.insert(std::pair<string, double>(marker_id, distance));
                }
                // next line
                getline(ifp, line);
            }
            //
            gmContig_all.insert(std::pair<int, map<string, double> >(group_num, gmContig));
            gmContigOrder_all.insert(std::pair<int, vector<string> >(group_num, gmContigOrder));
            //
            cout << "   Info: checking MSTmap linkage group " << group_num << " done. " << endl;            
            bool check_collect = false;
            if(check_collect)
            {
                cout << "          Details of corrected lg " << group_num << endl;
                vector<string>::iterator clgitr;
                vector<string>::iterator clgitr_end;
                clgitr     = gmContigOrder.begin();
                clgitr_end = gmContigOrder.end();
                while(clgitr != clgitr_end)
                {
                    string marker_id = *clgitr; // e.g., new_Chr1_9427397_10008496_le
                    assert(gmContig.find(marker_id) != gmContig.end() );
                    if(marker_pairing.find(marker_id) != marker_pairing.end() )
                    {
                        cout << "             " << marker_id << "\t" << gmContig[marker_id] << endl;
                    }else
                    {
                        cout << "#            " << marker_id << "\t" << gmContig[marker_id] << " # not paired. " << endl;
                    }
                    //
                    clgitr ++;
                }
            }
        }
    }
    // find the targeted number of linkage groups
    map<int, map<string, double> > gmContig_all_selected;         // <group_id, <marker_id, distance> >
    map<int, vector<string> >      gmContigOrder_all_selected;    // <group_id, <marker_id> >
    for(int i = 0; i < TOP_X_LG; i ++)
    {
        if(gmContig_all.size() > 0)
        {
            map<int, map<string, double> >::iterator lgitr;
            map<int, map<string, double> >::iterator lgitr_end;
            lgitr     = gmContig_all.begin();
            lgitr_end = gmContig_all.end();
            int max_marker_num = 0;
            int max_marker_lg  = 0;
            while(lgitr != lgitr_end)
            {
                int                group_num_tmp = (*lgitr).first;
                map<string, double> gmContig_tmp = (*lgitr).second;
                int               marker_num_tmp = gmContig_tmp.size();
                if(marker_num_tmp > max_marker_num)
                {
                    max_marker_num = marker_num_tmp;
                    max_marker_lg  = group_num_tmp;
                }
                //
                lgitr ++;
            }
            // collect 
            gmContig_all_selected.insert(std::pair<int, map<string, double> >(max_marker_lg, gmContig_all[max_marker_lg]));
            gmContigOrder_all_selected.insert(std::pair<int, vector<string> >(max_marker_lg, gmContigOrder_all[max_marker_lg]));
            // clean
            gmContig_all.erase(max_marker_lg);
            gmContigOrder_all.erase(max_marker_lg);
        }
    }
    //
    bool check_collect = true;    
    map<int, string> gmContig_best_marker_selected; // <group_id, best_marker_id> // well supported as phase 0    
    if(check_collect)
    {
        map<int, map<string, double> >::iterator lgitr;
        map<int, map<string, double> >::iterator lgitr_end;
        lgitr     = gmContig_all_selected.begin();
        lgitr_end = gmContig_all_selected.end();    
        int group_new_id = 0; 
        while(lgitr != lgitr_end)
        {
            map<string, double> gmContig      = (*lgitr).second;     
            vector<string>      gmContigOrder = gmContigOrder_all_selected[ (*lgitr).first ];                         
            //            
            group_new_id ++;
            cout << "   Details of selected corrected-lg (a) " << (*lgitr).first << " with new id " << group_new_id << endl;            
            std::stringstream ss;
            ss.str("");
            ss   << "intermediate_corrected_lg_" << group_new_id << "_map.txt"; // old id: (*lgitr).first
            ofstream ofp;
            ofp.open((ss.str()).c_str(), ios::out);
            if(!ofp.good())
            {
                cout << "   Error: cannot open file " << ss.str() << endl;
            }
            // raw group id = (*lgitr).first - 1
            ofp << "# corrected MSTmap linkage group " << group_new_id << " with old id: " << (*lgitr).first << endl;            
            //
            vector<string>::iterator clgitr;
            vector<string>::iterator clgitr_end;
            clgitr     = gmContigOrder.begin();
            clgitr_end = gmContigOrder.end();
            int numbering = 0;
            while(clgitr != clgitr_end)
            {
                string marker_id = *clgitr; // e.g., new_Chr1_9427397_10008496_le
                assert(gmContig.find(marker_id) != gmContig.end() );
                if(marker_pairing.find(marker_id) != marker_pairing.end() )
                {                    
                    cout << "      " << marker_id << "\t" << gmContig[marker_id] << endl;
                    numbering ++;
                    // format: new_Chr1_9427397_10008496_ri    2.095   ;   2
                    ofp  <<             marker_id << "\t" 
                         << gmContig[marker_id]   << "\t"
                         << ";"                   << "\t"
                         << numbering             << endl;
                    // marker genotype 
                    assert(marker_gt.find(marker_id) != marker_gt.end());
                    if(gmContig_best_marker_selected.find((*lgitr).first) == gmContig_best_marker_selected.end())
                    {
                        string gt_tmp = marker_gt[marker_id];
                        std::replace(gt_tmp.begin(), gt_tmp.end(), 'M', 'a'); // phase 0
                        std::replace(gt_tmp.begin(), gt_tmp.end(), 'P', 'b'); // phase 0
                        std::replace(gt_tmp.begin(), gt_tmp.end(), 'U', '-');                    
                        gmContig_best_marker_selected.insert(std::pair<int, string>((*lgitr).first, gt_tmp));
                    }else
                    {
                        string gt_tmp = marker_gt[marker_id];
                        std::replace(gt_tmp.begin(), gt_tmp.end(), 'M', 'a'); // phase 0
                        std::replace(gt_tmp.begin(), gt_tmp.end(), 'P', 'b'); // phase 0
                        std::replace(gt_tmp.begin(), gt_tmp.end(), 'U', '-');                        
                        string gt_exi = gmContig_best_marker_selected[ (*lgitr).first ];
                        //
                        size_t n_tmp = std::count(gt_tmp.begin(), gt_tmp.end(), '-');
                        size_t n_exi = std::count(gt_exi.begin(), gt_exi.end(), '-');
                        if(n_tmp < n_exi)
                        {
                            gmContig_best_marker_selected[ (*lgitr).first ] = gt_tmp;
                        }
                    }
                }else
                {
                    cout << "#     " << marker_id << "\t" << gmContig[marker_id] << " # not paired / discarded. " << endl;
                }
                //
                clgitr ++;
            }
            //
            ofp.close();
            //            
            lgitr ++;           
        }
    }
    // phase all marker, and check pairing
    map<string, int>  marker_pairing_2;    // <contig_id, count>, if count is 2, properly paired    
    if(check_collect)
    {
        map<int, map<string, double> >::iterator lgitr;
        map<int, map<string, double> >::iterator lgitr_end;
        lgitr     = gmContig_all_selected.begin();
        lgitr_end = gmContig_all_selected.end();     
        while(lgitr != lgitr_end)
        {
            map<string, double> gmContig      = (*lgitr).second;     
            vector<string>      gmContigOrder = gmContigOrder_all_selected[ (*lgitr).first ];        
            string best_marker_gt             = gmContig_best_marker_selected[ (*lgitr).first ];                 
            //
            cout << "   Details of selected corrected-lg (b) " << (*lgitr).first << endl;          
            //
            vector<string>::iterator clgitr;
            vector<string>::iterator clgitr_end;
            clgitr     = gmContigOrder.begin();
            clgitr_end = gmContigOrder.end();
            int numbering = 0;
            while(clgitr != clgitr_end)
            {
                string marker_id = *clgitr; // e.g., new_Chr1_9427397_10008496_le
                assert(gmContig.find(marker_id) != gmContig.end() );
                cout << "      " << marker_id << "\t" << gmContig[marker_id] << endl;
                numbering ++;
                // marker genotype 
                assert(marker_gt.find(marker_id) != marker_gt.end());
                string gt_tmp = marker_gt[marker_id];
                std::replace(gt_tmp.begin(), gt_tmp.end(), 'M', 'a'); // phase 0
                std::replace(gt_tmp.begin(), gt_tmp.end(), 'P', 'b'); // phase 0
                std::replace(gt_tmp.begin(), gt_tmp.end(), 'U', '-');
                string gt_tmp_flip = marker_gt[marker_id];
                std::replace(gt_tmp_flip.begin(), gt_tmp_flip.end(), 'M', 'b'); // phase 1
                std::replace(gt_tmp_flip.begin(), gt_tmp_flip.end(), 'P', 'a'); // phase 1
                std::replace(gt_tmp_flip.begin(), gt_tmp_flip.end(), 'U', '-'); 
                //
                double freq_tmp      = calculate_recomb_freq(gt_tmp,      best_marker_gt);    
                double freq_tmp_flip = calculate_recomb_freq(gt_tmp_flip, best_marker_gt);
                cout << "   " << best_marker_gt << " best_marker_gt \n"
                     << "   " << gt_tmp         << " freq_tmp,score======" << freq_tmp      << "\n"
                     << "   " << gt_tmp_flip    << " freq_tmp_flip,score=" << freq_tmp_flip << endl;
                if(freq_tmp < freq_tmp_flip && freq_tmp < max_recomb)
                {
                    // check pair
                    string contig_id = marker_id.substr(0, marker_id.size()-3); // e.g., new_Chr1_9427397_10008496
                    string leri      = marker_id.substr(marker_id.size()-2);                    
                    if(marker_pairing_2.find(contig_id) == marker_pairing_2.end() )
                    {
                        marker_pairing_2.insert(std::pair<string, int>(contig_id, 1));
                    }else
                    {
                        marker_pairing_2[contig_id] += 1;
                    }                                                            
                }else
                if(freq_tmp >= freq_tmp_flip && freq_tmp_flip < max_recomb)                
                {
                    // check pair
                    string contig_id = marker_id.substr(0, marker_id.size()-3); // e.g., new_Chr1_9427397_10008496
                    string leri      = marker_id.substr(marker_id.size()-2);                    
                    if(marker_pairing_2.find(contig_id) == marker_pairing_2.end() )
                    {
                        marker_pairing_2.insert(std::pair<string, int>(contig_id, 1));
                    }else
                    {
                        marker_pairing_2[contig_id] += 1;
                    }                                      
                }else
                {
                    ; // discarded cases
                }
                //
                clgitr ++;
            }
            //       
            lgitr ++;           
        }
    }    
    // ouput finally paired contig end markers 
    map<string, int> known_phase_files;
    map<int, map<string, int> > known_phase; // <grp_id, <marker_id, phase> >
    map<string, int> known_phase_dict; // <marker_id, phase>
    if(check_collect)
    {
        map<int, map<string, double> >::iterator lgitr;
        map<int, map<string, double> >::iterator lgitr_end;
        lgitr     = gmContig_all_selected.begin();
        lgitr_end = gmContig_all_selected.end();   
        int group_new_id = 0;  
        while(lgitr != lgitr_end)
        {
            map<string, double> gmContig      = (*lgitr).second;     
            vector<string>      gmContigOrder = gmContigOrder_all_selected[ (*lgitr).first ];        
            string best_marker_gt             = gmContig_best_marker_selected[ (*lgitr).first ];                 
            //
            group_new_id ++;
            if(known_phase.find(group_new_id) == known_phase.end())
            {
                map<string, int> tmp_phase_map;
                known_phase.insert(std::pair<int, map<string, int> >(group_new_id, tmp_phase_map));
            }
            cout << "   Details of selected corrected-lg (c) " << (*lgitr).first << " with new id " << group_new_id << endl;                  
            std::stringstream ss1;
            ss1.str("");
            ss1   << "intermediate_corrected_lg_" << group_new_id << "_marker_gt_phase.txt";
            known_phase_files.insert(std::pair<string, int>(ss1.str(), 1)); // collect known phase file
            ofstream ofp1;
            ofp1.open((ss1.str()).c_str(), ios::out);
            if(!ofp1.good())
            {
                cout << "   Error: cannot open file " << ss1.str() << endl;
            }            
            //
            vector<string>::iterator clgitr;
            vector<string>::iterator clgitr_end;
            clgitr     = gmContigOrder.begin();
            clgitr_end = gmContigOrder.end();
            int numbering = 0;
            while(clgitr != clgitr_end)
            {
                string marker_id = *clgitr; // e.g., new_Chr1_9427397_10008496_le
                string contig_id = marker_id.substr(0, marker_id.size()-3); // e.g., new_Chr1_9427397_10008496
                string leri      = marker_id.substr(marker_id.size()-2);     
                if(marker_pairing_2.find(contig_id) != marker_pairing_2.end() && marker_pairing_2[contig_id]==2)
                {                    
                    assert(gmContig.find(marker_id) != gmContig.end() );
                    cout << "      " << marker_id << "\t" << gmContig[marker_id] << endl;
                    numbering ++;
                    // marker genotype 
                    assert(marker_gt.find(marker_id) != marker_gt.end());
                    string gt_tmp = marker_gt[marker_id];
                    std::replace(gt_tmp.begin(), gt_tmp.end(), 'M', 'a'); // phase 0
                    std::replace(gt_tmp.begin(), gt_tmp.end(), 'P', 'b'); // phase 0
                    std::replace(gt_tmp.begin(), gt_tmp.end(), 'U', '-');
                    string gt_tmp_flip = marker_gt[marker_id];
                    std::replace(gt_tmp_flip.begin(), gt_tmp_flip.end(), 'M', 'b'); // phase 1
                    std::replace(gt_tmp_flip.begin(), gt_tmp_flip.end(), 'P', 'a'); // phase 1
                    std::replace(gt_tmp_flip.begin(), gt_tmp_flip.end(), 'U', '-'); 
                    //
                    double freq_tmp      = calculate_recomb_freq(gt_tmp,      best_marker_gt);    
                    double freq_tmp_flip = calculate_recomb_freq(gt_tmp_flip, best_marker_gt);
                    cout << "   " << best_marker_gt << " best_marker_gt \n"
                         << "   " << gt_tmp         << " freq_tmp,score======" << freq_tmp      << "\n"
                         << "   " << gt_tmp_flip    << " freq_tmp_flip,score=" << freq_tmp_flip << endl;
                    if(freq_tmp < freq_tmp_flip && freq_tmp < max_recomb)
                    {
                        ofp1 << numbering << "\t" << marker_id << "\t{0}\t";
                        for(int i = 0; i < gt_tmp.size(); i ++)
                        {
                            ofp1 << gt_tmp.substr(i, 1) << "\t";
                        }
                        ofp1 << endl;   
                        //
                        known_phase[group_new_id].insert(std::pair<string, int>(marker_id, 0));         
                        known_phase_dict.insert(std::pair<string, int>(marker_id, 0));
                    }else
                    if(freq_tmp >= freq_tmp_flip && freq_tmp_flip < max_recomb)                
                    {
                        ofp1 << numbering << "\t" << marker_id << "\t{1}\t";
                        for(int i = 0; i < gt_tmp_flip.size(); i ++)
                        {
                            ofp1 << gt_tmp.substr(i, 1) << "\t"; // keep it as raw, required by following steps
                        } 
                        ofp1 << endl;          
                        known_phase[group_new_id].insert(std::pair<string, int>(marker_id, 1));    
                        known_phase_dict.insert(std::pair<string, int>(marker_id, 1));
                    }else
                    {
                        cout << "discarding " << marker_id <<  endl; // discarded cases
                    }
                }
                //
                clgitr ++;
            }
            //
            ofp1.close();
            //            
            lgitr ++;           
        }
    } 
    // output finally selected genetic map
    if(check_collect)
    {
        map<int, map<string, double> >::iterator lgitr;
        map<int, map<string, double> >::iterator lgitr_end;
        lgitr     = gmContig_all_selected.begin();
        lgitr_end = gmContig_all_selected.end();   
        int group_new_id = 0;  
        while(lgitr != lgitr_end)
        {
            map<string, double> gmContig      = (*lgitr).second;     
            vector<string>      gmContigOrder = gmContigOrder_all_selected[ (*lgitr).first ];                         
            //
            group_new_id ++;
            cout << "   Details of selected corrected-lg (d) " << (*lgitr).first << " with new id " << group_new_id << endl;                  
            //
            std::stringstream ss;
            ss.str("");
            ss   << "final_corrected_lg_" << group_new_id << "_map.txt";
            ofstream ofp;
            ofp.open((ss.str()).c_str(), ios::out);
            if(!ofp.good())
            {
                cout << "   Error: cannot open file " << ss.str() << endl;
            }
             // raw group id = (*lgitr).first - 1
            ofp << "# corrected MSTmap linkage group " << group_new_id << " with old id: " << (*lgitr).first << endl;
            //
            vector<string>::iterator clgitr;
            vector<string>::iterator clgitr_end;
            clgitr     = gmContigOrder.begin();
            clgitr_end = gmContigOrder.end();
            int numbering = 0;
            while(clgitr != clgitr_end)
            {
                string marker_id = *clgitr; // e.g., new_Chr1_9427397_10008496_le
                string contig_id = marker_id.substr(0, marker_id.size()-3); // e.g., new_Chr1_9427397_10008496
                string leri      = marker_id.substr(marker_id.size()-2);  
                if(marker_pairing_2.find(contig_id) != marker_pairing_2.end() && marker_pairing_2[contig_id]==2)
                {                                 
                    assert(gmContig.find(marker_id) != gmContig.end() );
                    if(marker_pairing.find(marker_id) != marker_pairing.end() )
                    {                    
                        cout << "      " << marker_id << "\t" << gmContig[marker_id] << endl;
                        numbering ++;
                        // format: new_Chr1_9427397_10008496_ri    2.095   ;   2
                        ofp  <<             marker_id << "\t" 
                             << gmContig[marker_id]   << "\t"
                             << ";"                   << "\t"
                             << numbering             << endl;
                    }else
                    {
                        cout << "#     " << marker_id << "\t" << gmContig[marker_id] << " # not paired / discarded. " << endl;
                    }
                }
                //
                clgitr ++;
            }
            //
            ofp.close();
            //            
            lgitr ++;           
        }
    }           
    //
    ifp.close();
    cout << "   Info: 1st round - integrate phase of contigs not in final map " << endl;
    bool      output_phase = false;
    double max_recomb_freq = 20;
    map<int, map<string, int> > known_phase_updated_1st;
    if(!phase_other_contig(known_phase, 
                           &known_phase_dict, 
                           marker_gt, 
                           output_phase, 
                           max_recomb_freq,
                           &known_phase_updated_1st))
    {
        return false;
    }
    cout << "   Info: 2nd round: integrate phase of other contigs, considering also added contig markers from 1st round. " << endl;
    output_phase = true;
    max_recomb_freq = 35; // more tolerate    
    map<int, map<string, int> > known_phase_updated_2nd;
    if(!phase_other_contig(known_phase_updated_1st, 
                           &known_phase_dict, 
                           marker_gt, 
                           output_phase, 
                           max_recomb_freq,                           
                           &known_phase_updated_2nd))
    {
        return false;
    }    
    //
    return 0;
}
//
bool phase_other_contig(map<int, map<string, int> >  known_phase, 
                        map<string, int>*            known_phase_dict, 
                        map<string, string>          marker_gt, 
                        bool                         output_phase,
                        double                       max_recomb_freq,
                        map<int, map<string, int> >* known_phase_updated)
{
    /*
        map<string, string> marker_gt ............................ <contig_id + "_le/ri", raw_PM_pattern>
        map<int, map<string, int> > known_phase ... <group_new_id, <contig_id + "_le/ri", 0/1> >
        map<string, int> known_phase_dict......................... <contig_id + "_le/ri", 0/1>
    */
    // map<int, map<string, int> > known_phase_updated = known_phase;
    *known_phase_updated = known_phase;
    //
    map<string, string>::iterator mitr;
    map<string, string>::iterator mitr_end;
    mitr     = marker_gt.begin();
    mitr_end = marker_gt.end();
    while(mitr != mitr_end)
    {
        string marker_id = (*mitr).first;
        //
        if((*known_phase_dict).find(marker_id)!= (*known_phase_dict).end())
        {
            mitr ++;
            continue;
        }
        //
        string marker_pm = (*mitr).second;
        string marker_pm_flipped = marker_pm;
        std::replace(marker_pm_flipped.begin(), marker_pm_flipped.end(), 'M', 'p'); // M->p
        std::replace(marker_pm_flipped.begin(), marker_pm_flipped.end(), 'P', 'm'); // P->m
        std::transform(marker_pm_flipped.begin(), marker_pm_flipped.end(), marker_pm_flipped.begin(), ::toupper);        
        // find best match in existing phasing files.
        int best_matching_group = 1;
        int best_matching_score = 100;
        int best_matching_phase = -1;
        string best_matching_marker    = "";
        int best_matching_marker_phase = -1;        
        map<int, map<string, int> >::iterator gitr;
        map<int, map<string, int> >::iterator gitr_end;
        gitr     = known_phase.begin();
        gitr_end = known_phase.end();
        while(gitr != gitr_end)
        {
            int this_group = (*gitr).first;
            map<string, int> tmp_phase_map = (*gitr).second;
            map<string, int>::iterator pmitr;
            map<string, int>::iterator pmitr_end;
            pmitr     = tmp_phase_map.begin();
            pmitr_end = tmp_phase_map.end();
            while(pmitr != pmitr_end)
            {
                string marker_id_phased = (*pmitr).first;
                int    marker_phase     = (*pmitr).second;
                //
                assert( marker_gt.find(marker_id_phased) != marker_gt.end() );
                string marker_pm_phased = marker_gt[ marker_id_phased ];
                if(marker_phase == 1)
                {
                    std::replace(marker_pm_phased.begin(),marker_pm_phased.end(), 'M', 'p'); // M->p
                    std::replace(marker_pm_phased.begin(),marker_pm_phased.end(), 'P', 'm'); // P->m
                    std::transform(marker_pm_phased.begin(),marker_pm_phased.end(),marker_pm_phased.begin(),::toupper);                        
                }
                //
                double tmp_co_freq = calculate_recomb_freq(marker_pm_phased, marker_pm);
              //cout << "   tmp_co_freq=" << tmp_co_freq << " vs best_matching_score=" << best_matching_score << endl;                
                if(tmp_co_freq < best_matching_score)
                {
                  //cout << "   tmp_co_freq=" << tmp_co_freq << " vs best_matching_score=" << best_matching_score << endl;
                    best_matching_score  = tmp_co_freq;
                    best_matching_group  = this_group;
                    best_matching_phase  = 0;
                    best_matching_marker       = marker_id_phased;
                    best_matching_marker_phase = marker_phase;
                }
                tmp_co_freq = calculate_recomb_freq(marker_pm_phased, marker_pm_flipped);
                if(tmp_co_freq < best_matching_score)
                {
                  //cout << "   tmp_co_freq=" << tmp_co_freq << " vs best_matching_score=" << best_matching_score << endl;                
                    best_matching_score = tmp_co_freq;
                    best_matching_group = this_group;
                    best_matching_phase = 1;
                    best_matching_marker       = marker_id_phased;
                    best_matching_marker_phase = marker_phase;                    
                }                
                //
                pmitr ++;
            }
            //
            gitr ++;
        }
        if(best_matching_score < max_recomb_freq)
        {
            cout << "   check_phase: " << marker_id 
                 << " with phase "     << best_matching_phase
                 << " best matched "   << best_matching_marker 
                 << " with phase "     << best_matching_marker_phase
                 << " in group "       << best_matching_group 
                 << " with score "     << best_matching_score << endl;
            cout << "                " << marker_id << endl
                 << "                  raw....... PM:" << marker_pm << endl
                 << "                  flipped... PM:" << marker_pm_flipped << endl
                 << "                  raw.known. PM:" << marker_gt[best_matching_marker] << endl;
            // update
            assert( (*known_phase_updated).find(best_matching_group) != (*known_phase_updated).end() );
            (*known_phase_updated)[best_matching_group].insert(std::pair<string, int>(marker_id, best_matching_phase));
            //
            (*known_phase_dict).insert(std::pair<string, int>(marker_id, best_matching_phase));
        }
        //
        mitr ++;
    }
    // output phases of markers (note: some markers showing more recombinations might be have excluded)
    if(output_phase)
    {
        map<int, map<string, int> >::iterator gitr;
        map<int, map<string, int> >::iterator gitr_end;
        gitr     = (*known_phase_updated).begin();
        gitr_end = (*known_phase_updated).end();
        while(gitr != gitr_end)
        {
            int this_group = (*gitr).first;
            map<string, int> tmp_phase_map = (*gitr).second;
            //
            std::stringstream ss;
            ss.str("");
            ss   << "final_corrected_lg_" << this_group << "_marker_gt_phase.txt";            
            ofstream ofp;
            ofp.open((ss.str()).c_str(), ios::out);
            if(!ofp.good())
            {
                cout << "   Error: cannot open file " << ss.str() << endl;
            }
             // raw group id = (*lgitr).first - 1
            ofp << "# corrected+integrated MSTmap linkage group " << this_group << " with old id not collected here " << endl;            
            //
            map<string, int>::iterator pmitr;
            map<string, int>::iterator pmitr_end;
            pmitr     = tmp_phase_map.begin();
            pmitr_end = tmp_phase_map.end();
            int numbering = 0;
            while(pmitr != pmitr_end)
            {
                string marker_id_phased = (*pmitr).first;
                int    marker_phase     = (*pmitr).second;   
                assert( marker_gt.find(marker_id_phased) != marker_gt.end() );
                string gt_tmp           = marker_gt[ marker_id_phased ];
                //
                std::replace(gt_tmp.begin(), gt_tmp.end(), 'M', 'a'); // phase 0
                std::replace(gt_tmp.begin(), gt_tmp.end(), 'P', 'b'); // phase 0
                std::replace(gt_tmp.begin(), gt_tmp.end(), 'U', '-');   
                ofp << numbering << "\t" << marker_id_phased << "\t{"<< marker_phase << "}\t";
                for(int i = 0; i < gt_tmp.size(); i ++)
                {
                    ofp << gt_tmp.substr(i, 1) << "\t";
                }
                ofp << endl;                                                                         
                //
                pmitr ++;
            } 
            //
            ofp.close();           
            //
            gitr ++;
        }        
    }
    //
    return true;
}
/* convert MSTmap format into JoinMap format
tig00003853_ri          0.000  ;  132
...
*/
//
double calculate_recomb_freq(string PMpat1, string PMpat2)
{
    double recomb_freq = 0.0;
    string pm1 = PMpat1;
    string pm2 = PMpat2;
    std::transform(pm1.begin(), pm1.end(), pm1.begin(), ::toupper);    
    std::transform(pm2.begin(), pm2.end(), pm2.begin(), ::toupper);      
    // only consider effective exact P/M positions
    unsigned long matchedp   = 0;
    unsigned long effectivep = 0;
    assert(pm1.size() == pm2.size());
    for(unsigned long ii = 0; ii < pm1.size(); ii ++)
    {
        string letter1 = pm1.substr(ii, 1);
        string letter2 = pm2.substr(ii, 1);
        //
        if(letter1.compare("-")!=0 && letter2.compare("-") != 0)
        {
            effectivep ++;
        } 
        //
        if( (letter1.compare("B")==0 && letter2.compare("B") == 0) || 
            (letter1.compare("P")==0 && letter2.compare("P") == 0) )
        {
            matchedp ++;
        }else
        if( (letter1.compare("A")==0 && letter2.compare("A") == 0) ||
            (letter1.compare("M")==0 && letter2.compare("M") == 0) )
        {
            matchedp ++;
        }else ;
    }
    recomb_freq = (effectivep-matchedp)*1.0/effectivep*100;
    /*
    //
    cout << PMpat1 << endl
         << PMpat2 << endl;
    //
    cout << "   Info: effectivep.=" << effectivep              << endl
         << "         recombined.=" << effectivep - matchedp   << endl
         << "         recomb_freq=" << recomb_freq             << endl;
    */
    return recomb_freq;
}
//
bool read_gb_genotype(string gtfile, map<string, string>* marker_gt)
{
    ifstream ifp;
    ifp.open(gtfile, ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << gtfile << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        if(line.find("CO-info") != std::string::npos) continue;
        cout << line << endl;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 3)
        {
            cout << "   Warning: skipping insufficient line info at " << line << endl;
            continue;
        }
        string gt        = lineinfo[0];
        string contig_id = lineinfo[1];
        string leri      = "";
        if(lineinfo[2].find("left") != std::string::npos)
        {
            leri = "le";
        }else
        if(lineinfo[2].find("right") != std::string::npos)
        {
            leri = "ri";
        }else        
        {
            cout << "   Error: unexpected info of " << lineinfo[2] << endl;
            return false;
        }
        string marker_id = contig_id + "_" + leri;
        (*marker_gt).insert(std::pair<string, string>(marker_id, gt));
    }
    ifp.close();
    return true;
}
// std::replace(gt.begin(), gt.end(), 'M', 'a'); // phase 0
// std::replace(gt.begin(), gt.end(), 'P', 'b'); // phase 0
// std::replace(gt.begin(), gt.end(), 'U', '-');
// std::replace(gt.begin(), gt.end(), 'M', 'b'); // phase 1
// std::replace(gt.begin(), gt.end(), 'P', 'a'); // phase 1


















