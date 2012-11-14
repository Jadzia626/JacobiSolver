/*
**  Class : Jacobi
** ~~~~~~~~~~~~~~~~
*/

namespace eigensolvers {

class Jacobi {

    public:

    void Simple();
    void Cyclic();
    void Parallel();

    private:

    double fFrobenius();
    void   fSearch();

};

} // End namespace eigensolvers
