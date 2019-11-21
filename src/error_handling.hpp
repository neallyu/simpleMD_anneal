#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include <stdexcept>

using namespace std;

struct ReadingFile_Other: public exception
{
    const char *what() const throw () {
        return "other problems";
    }
};


struct ReadingFile_Open: public exception
{
    const char *what() const throw () {
        return "cannot open input file";
    }
};


struct CreatingOutputPath: public exception
{
    const char *what() const throw () {
        return "creating output path";
    }
};

#endif