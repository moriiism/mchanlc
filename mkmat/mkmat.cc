#include "mi_iolib.h"
#include "mir_graph2d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_mkmat.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[]){
    int status = kRetNormal;
  
    ArgValMkmat* argval = new ArgValMkmat;
    argval->Init(argc, argv);
    argval->Print(stdout);

    //
    // freq_info
    //
    long nline = 0;
    string* line_arr = NULL;
    MiIolib::GenReadFileSkipComment(argval->GetFreqInfoFile(),
                                    &line_arr, &nline);
    printf("nline = %ld\n", nline);
    if(1 != nline){
        printf("bad freq_info.txt\n");
        abort();
    }
    int nsplit = 0;
    string* split_arr = NULL;
    MiStr::GenSplit(line_arr[0], &nsplit, &split_arr);
    double freq_lo = atof(split_arr[0].c_str());
    double freq_up = atof(split_arr[1].c_str());
    int nfreq   = atoi(split_arr[2].c_str());
    delete [] line_arr;
    delete [] split_arr;
    
    //
    // data file
    //
    GraphDataNerr2d* gd2d = new GraphDataNerr2d;
    gd2d->Load(argval->GetInfile());
    gd2d->Sort();
    int sorted = MirMath::IsSorted(gd2d->GetNdata(),
                                   gd2d->GetXvalArr()->GetVal());
    printf("sorted ? = %d\n", sorted);
    MirQdpTool::MkQdp(gd2d, "temp.qdp", "x,y");

    printf("ndata = %ld\n", gd2d->GetNdata());

    // matrix
    // (nrow, ncol)
    long nrow = gd2d->GetNdata();
    long ncol = nfreq * 2;
    double** mat = new double* [nrow];
    for(long irow = 0; irow < nrow; irow ++){
        mat[irow] = new double [ncol];
        for(long icol = 0; icol < ncol; icol ++){
            mat[irow][icol] = 0.0;
        }
    }
    double delta_freq = (freq_up - freq_lo) / nfreq;
    for(long irow = 0; irow < nrow; irow ++){
        for(long icol = 0; icol < ncol/2; icol ++){
            double freq = freq_lo + delta_freq * (icol - 0.5);
            double time = gd2d->GetXvalElm(irow);
            mat[irow][icol]          = 2 * delta_freq * cos(2 * M_PI * freq * time);
            mat[irow][ncol/2 + icol] = 2 * delta_freq * sin(2 * M_PI * freq * time);
        }
    }

    FILE* fp = fopen(argval->GetOutfile().c_str(), "w");
    for(long irow = 0; irow < nrow; irow ++){
        for(long icol = 0; icol < ncol; icol ++){
            fprintf(fp, "%e ", mat[irow][icol]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    
    delete argval;
    
    return status;
}

