#include <string>
 
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>

struct dna_deg : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna15; // instead of dna5
    template <typename alph>
    using sequence_container = seqan3::dna15_vector;
};
 

int main ()
{
    std::filesystem::path base_dir{"D:/PMS/CoV"}; // get the directory
    seqan3::debug_stream <<1<< base_dir<<"\n";
 
    // Initialise a file input object with a FASTA file.
    seqan3::sequence_file_input<dna_deg> file_in{base_dir/"allnuc0209.fasta", seqan3::format_fasta{}};
    seqan3::sequence_file_output<dna_deg> file_E{base_dir/"E.20220209.fasta"};
    seqan3::sequence_file_output<dna_deg> file_N{base_dir/"N.20220209.fasta"};
    seqan3::sequence_file_output<dna_deg> file_S{base_dir/"S.20220209.fasta"};

    long e{0L}, n{0L}, s{0L}, t{0L}, m{(1L<<15)-1};
    seqan3::debug_stream << "m= " << m << "\n" ; 
    
    for (auto & record : file_in)
    {
             if (record.id().starts_with("E|"))        {file_E.push_back(std::move(record));++e;}
        else if (record.id().starts_with("N|"))        {file_N.push_back(std::move(record));++n;}
        else if (record.id().starts_with("Spike|"))    {file_S.push_back(std::move(record));++s;}
        if (!(++t & m)) seqan3::debug_stream << "N= " << n << ", E= " << e << ", Spike= " << s << "\n" ; 
    }
    seqan3::debug_stream << "T= " << t  << ". N= " << n << ", E= " << e << ", Spike= " << s << "\n" ; 
 
    return 0;
}