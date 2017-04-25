#include "mi_iolib.h"
#include "mir_graph2d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_mkmat2.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[]){
    int status = kRetNormal;
  
    ArgValMkmat2* argval = new ArgValMkmat2;
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
    GraphDataNerr2d* gd2d_1 = new GraphDataNerr2d;
    gd2d_1->Load(argval->GetInfile1());
    gd2d_1->Sort();
    int sorted = MirMath::IsSorted(gd2d_1->GetNdata(),
                                   gd2d_1->GetXvalArr()->GetVal());
    printf("sorted ? = %d\n", sorted);
    MirQdpTool::MkQdp(gd2d_1, "temp1.qdp", "x,y");

    printf("ndata1 = %ld\n", gd2d_1->GetNdata());


    GraphDataNerr2d* gd2d_2 = new GraphDataNerr2d;
    gd2d_2->Load(argval->GetInfile2());
    gd2d_2->Sort();
    sorted = MirMath::IsSorted(gd2d_2->GetNdata(),
                               gd2d_2->GetXvalArr()->GetVal());
    printf("sorted ? = %d\n", sorted);
    MirQdpTool::MkQdp(gd2d_2, "temp2.qdp", "x,y");

    printf("ndata2 = %ld\n", gd2d_2->GetNdata());
    

    // matrix for data1
    // (nrow_1, ncol_1)
    long nrow_1 = gd2d_1->GetNdata();
    long ncol_1 = nfreq * 2;
    double** mat_1 = new double* [nrow_1];
    for(long irow = 0; irow < nrow_1; irow ++){
        mat_1[irow] = new double [ncol_1];
        for(long icol = 0; icol < ncol_1; icol ++){
            mat_1[irow][icol] = 0.0;
        }
    }
    double delta_freq = (freq_up - freq_lo) / nfreq;
    for(long irow = 0; irow < nrow_1; irow ++){
        for(long icol = 0; icol < ncol_1/2; icol ++){
            double freq = freq_lo + delta_freq * (icol - 0.5);
            double time = gd2d_1->GetXvalElm(irow);
            mat_1[irow][icol]            = 2 * delta_freq * cos(2 * M_PI * freq * time);
            mat_1[irow][ncol_1/2 + icol] = 2 * delta_freq * sin(2 * M_PI * freq * time);
        }
    }

    // matrix for data2
    // (nrow_2, ncol_2)
    long nrow_2 = gd2d_2->GetNdata();
    long ncol_2 = nfreq * 2;
    double** mat_2 = new double* [nrow_2];
    for(long irow = 0; irow < nrow_2; irow ++){
        mat_2[irow] = new double [ncol_2];
        for(long icol = 0; icol < ncol_2; icol ++){
            mat_2[irow][icol] = 0.0;
        }
    }
    delta_freq = (freq_up - freq_lo) / nfreq;
    for(long irow = 0; irow < nrow_2; irow ++){
        for(long icol = 0; icol < ncol_2/2; icol ++){
            double freq = freq_lo + delta_freq * (icol - 0.5);
            double time = gd2d_2->GetXvalElm(irow);
            mat_2[irow][icol]            = 2 * delta_freq * cos(2 * M_PI * freq * time);
            mat_2[irow][ncol_2/2 + icol] = 2 * delta_freq * sin(2 * M_PI * freq * time);
        }
    }

    // matrix
    // (nrow, ncol)
    long nrow = gd2d_1->GetNdata() + gd2d_2->GetNdata();
    long ncol = nfreq * 2 * 2;
    double** mat = new double* [nrow];
    for(long irow = 0; irow < nrow; irow ++){
        mat[irow] = new double [ncol];
        for(long icol = 0; icol < ncol; icol ++){
            mat[irow][icol] = 0.0;
        }
    }

    for(long irow = 0; irow < nrow_1; irow ++){
        for(long icol = 0; icol < ncol_1; icol ++){
            mat[irow][icol] = mat_1[irow][icol];
        }
    }
    for(long irow = 0; irow < nrow_2; irow ++){
        for(long icol = 0; icol < ncol_2; icol ++){
            mat[nrow_1 + irow][ncol_1 + icol] = mat_2[irow][icol];
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

