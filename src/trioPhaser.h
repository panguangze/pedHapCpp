//
// Created by caronkey on 10/23/21.
//

#ifndef PEDHAP_TRIOPHASER_H
#define PEDHAP_TRIOPHASER_H

#include "type.h"
class TrioPhaser {
public:
    ChromoPhaser *chromoPhaser;
    TrioPhaser() = default;
    ~TrioPhaser() = default;
    inline void setChromoPhase(ChromoPhaser* chromoPhaser){
        this->chromoPhaser = chromoPhaser;
    }
};


#endif //PEDHAP_TRIOPHASER_H
