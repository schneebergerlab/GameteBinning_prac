/*////////////////////////////////////////////////////////////////////////////////////////////////////

function.: convert vcf format into plain text, count reads on reference and alternative alleles.
Writen by: Hequan Sun, MPIPZ, Germany
Email....: sun@mpipz.mpg.de/sunhequan@gmail.com

///////////////////////////////////////////////////////////////////////////////////////////////////*/
#include             <stddef.h>
#include             <stdlib.h>
#include              <stdio.h>
#include               <math.h>
#include               <string>
#include             <string.h>
#include              <ctype.h>
#include               <time.h>
#include            "globals.h"
#include   "allele_counter_extract.h"
#include   "allele_counter_convert.h"

using namespace std;
// declare sub functions
void init_global(int argc, char* argv[]);
//
int main(int argc, char* argv[])
{
    double startT  = clock();
    init_global(argc, argv);
    double finishT = clock();
    
    printf("\nTime consumed %.4f seconds.\n", (finishT-startT)/1000000.0);
    return 0;
}
// define functions
void init_global(int argc, char* argv[])
{
    // this function read command line parameters
    string version = "1.1";
    std::string usage_global = "";
    usage_global += "Usage of allele_counter version " + version + ": allele_counter FUNCTION [OPTIONS]\n\n";
    usage_global += " FUNCTION is one of:\n\n";
    usage_global += "  extract  \textracts  base call info     according to markers \n";
    usage_global += "  convert  \tconverts  VCF4.1   into     acceptable formats; or \n";
    usage_global += "           \tconverts  shore pileup into shore consensus_summary format \n";

    if (argc < 2) 
    {
        printf("\nNOT enough input parameters. Check necessary parameters below: \n\n");
        printf("%s\n", usage_global.c_str());
        exit(1);
    }
    
    char* myTask = argv[1];
    if(!strcmp(myTask, "convert"))
    {
        if(argc > 2)
        printf("\nConverting format: \n\n");
        strcatCMD += "convert ";
        allele_counter_convert(argc, argv);
        printf("\nConvertion done.\n");
    }
    else if(!strcmp(myTask, "extract"))
    {
        if(argc > 2)
        printf("\nExtracting consensus according to markers: \n");
        strcatCMD += "extract ";
        allele_counter_extract(argc, argv);
        printf("\nEXtraction done.\n");
    }
    else
    {
        printf("\nINVALID command: %s. Check COMMANDs below: \n\n", myTask);
        printf("%s\n", usage_global.c_str());
    }
}
