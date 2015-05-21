#include "testspace.H"

namespace testspace
{
    void deepthought(const int num)
    {
        answer = num;
    }

    void ask()
    {
        std::cout<< "What is the meaning of Life, the Universe, and Everything? " << answer << "\n";
    }

}
