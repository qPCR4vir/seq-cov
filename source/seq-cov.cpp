#include <string>
#include <filesystem>
 
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>
#include <seqan3/argument_parser/all.hpp> 

#include <nana/gui.hpp>
#include <nana/gui/widgets/button.hpp>
#include <nana/gui/widgets/group.hpp>
#include <nana/gui/widgets/label.hpp>
#include <nana/gui/widgets/progress.hpp>
#include <nana/gui/widgets/textbox.hpp>
#include <nana/gui/widgets/combox.hpp>
#include <nana/gui/filebox.hpp>
#include <nana/gui/msgbox.hpp>
#include <nana/gui/tooltip.hpp>
#include <nana/gui/place.hpp>

struct dna_deg : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna15; // instead of dna5
    template <typename alph>
    using sequence_container = seqan3::dna15_vector;
};
 
void split_fasta(const std::filesystem::path& fasta, std::string gene)
{
    std::cout << "Going to split: " << fasta.string() << ", gene: " << gene ;
    return ;
    /*
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
    */
}

int main (int argc, char ** argv)
{
    seqan3::argument_parser arg_parser{"Split-CoV-fasta", argc, argv}; 
    std::filesystem::path base_dir{"D:/PMS/CoV"}; // get the directory
    // seqan3::debug_stream <<1<< base_dir<<"\n";
    
    nana::form fm{nana::api::make_center(1000, 200)};
    fm.caption("Split-CoV-fasta");
    nana::label input_file_label{fm, "Original fasta file:"};
    nana::textbox input_file{fm};
    input_file.tip_string("Original fasta file:");
    nana::combox gene{fm};
    gene.push_back("E").push_back("N").push_back("Spike");
    nana::button set{fm, "&Set"}, 
                 run_split{fm, "&Split"};
    run_split.events().click([&](){split_fasta(input_file.text(), gene.caption());});
    set.events().click([&]()
    {
        nana::filebox fb{fm, true};
        fb.title("Select the original GISAID CoV fasta file");
        fb.add_filter("fasta file", "*.fasta");
        const auto&files = fb.show();
        if (!files.empty())
           input_file.reset(files[0].string());
    });
    auto& p = fm.get_place();
    p.div(R"(vertical  <height=30 input arrange=[120,variable,40]> 
                       <height=10 >
                       <height=30 split arrange=[100,80,40] gap=10> 
                       margin=10)");
    p["input"] << input_file_label << input_file << set ;
    p["split"] << "Select the gene: " << gene << run_split ;
    p.collocate();



    fm.show();
    nana::exec();
    return 0;
}