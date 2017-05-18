#ifndef MORIIISM_MCHANLC_MKMAT_ARG_MKMAT2_H_
#define MORIIISM_MCHANLC_MKMAT_ARG_MKMAT2_H_

#include "mi_base.h"

class ArgValMkmat2 : public MiArgBase{
public:
    ArgValMkmat2() :
        MiArgBase(),
        progname_(""),
        infile1_(""),
        infile2_(""),
        freq_info_file_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkmat2(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile1()   const {return infile1_;};
    string GetInfile2()   const {return infile2_;};    
    string GetFreqInfoFile() const {return freq_info_file_;};
    string GetOutdir()   const {return outdir_;};
    string GetOutfileHead()   const {return outfile_head_;};

private:
    string progname_;
    string infile1_;
    string infile2_;
    string freq_info_file_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);    
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MCHANLC_MKMAT_ARG_MKMAT2_H_
