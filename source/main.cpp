
#include "gui.hpp"
#include <seqan3/argument_parser/all.hpp> 

// grap most std::cout and seqan3::debug_stream with  if constexpr (debugging) {}


int main (int argc, char ** argv)
{
    seqan3::argument_parser arg_parser{"Split-CoV-fasta", argc, argv}; 
    
    GUI fm;
    fm.show();
    nana::exec();
    return 0;
}