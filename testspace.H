#ifndef TESTSPACE_H
#define TESTSPACE_H

#include "className.H"

namespace testspace
{
    NamespaceName("Test");

    int answer;
    void deepthought(const int num);
    void ask();

}

#ifdef NoRepository
#   include "testspace.C"
#endif

#endif
