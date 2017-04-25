#ifndef MORIIISM_MCHANLC_MKMAT_ARG_MKMAT_H_
#define MORIIISM_MCHANLC_MKMAT_ARG_MKMAT_H_

#include "mi_base.h"

class ArgValMkmat : public MiArgBase{
public:
    ArgValMkmat() :
        MiArgBase(),
        progname_(""),
        infile_(""),
        freq_info_file_(""),
        outfile_("") {}
    ~ArgValMkmat(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile()   const {return infile_;};
    string GetFreqInfoFile() const {return freq_info_file_;};
    string GetOutfile()   const {return outfile_;};

private:
    string progname_;
    string infile_;
    string freq_info_file_;
    string outfile_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);    
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MCHANLC_MKMAT_ARG_MKMAT_H_
