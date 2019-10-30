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
