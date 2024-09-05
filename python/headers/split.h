

#ifndef bbww_SPLIT
#define bbww_SPLIT

#include <stdlib.h>
#include <math.h> 

namespace split {
    bool MET(double met){
        met = abs(met);
        met *= 10000;
        met -= floor(met);
        int met_int = int((met*10));
        return met_int%2 == 0;
    }

}

#endif
