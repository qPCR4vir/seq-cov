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
 
void split_fasta(const std::filesystem::path& fasta, bool E, bool N, bool S)
{
    std::cout << "\nGoing to split: " << fasta.string() << ", gene E: " << E  
              << ", gene N: " << N << ", gene Spike: " << S   ;
    
    // Initialise a file input object with a FASTA file.
    seqan3::sequence_file_input file_in{fasta};
    auto dir=fasta.parent_path();
    auto name=fasta.filename().string();
    std::filesystem::path e_fasta = E ? dir / ("E." + name) : std::filesystem::path( "E.fasta");
    std::filesystem::path n_fasta = N ? dir / ("N." + name) : std::filesystem::path( "N.fasta");
    std::filesystem::path s_fasta = S ? dir / ("Spike." + name) : std::filesystem::path( "S.fasta");
    seqan3::sequence_file_output file_E{e_fasta};
    seqan3::sequence_file_output file_N{n_fasta};
    seqan3::sequence_file_output file_S{s_fasta};

    long e{0L}, n{0L}, s{0L}, t{0L}, m{(1L<<15)-1};
    seqan3::debug_stream << "m= " << m << "\n" ; 
    
    for (auto & record : file_in)
    {
             if (E && record.id().starts_with("E|"))        {file_E.push_back(std::move(record));++e;}
        else if (N && record.id().starts_with("N|"))        {file_N.push_back(std::move(record));++n;}
        else if (S && record.id().starts_with("Spike|"))    {file_S.push_back(std::move(record));++s;}
        if (!(++t & m)) seqan3::debug_stream << "N= " << n << ", E= " << e << ", Spike= " << s << "\n" ; 
    }
    seqan3::debug_stream << "T= " << t  << ". N= " << n << ", E= " << e << ", Spike= " << s << "\n" ; 
    
}

class GUI: public nana::form
{
    nana::label    input_file_label{*this, "Original fasta file:"};
    nana::textbox  input_file      {*this};
    nana::button   set             {*this, "&Select"}, 
                   run_split       {*this, "S&plit" };
    nana::group    gene            {*this, "Gene"   };
public:
    GUI() : nana::form{nana::api::make_center(1000, 150)}
    {
        caption("Split-CoV-fasta. v0.01.00");

        input_file.tip_string("Original fasta file:").multi_lines(false);
        
        auto& E = gene.add_option("E");
        auto& N = gene.add_option("N");
        auto& S = gene.add_option("Spike");

        run_split.events().click([&]()
        {
            split_fasta(input_file.text(), E.checked(), N.checked(), S.checked());
        });

        set.events().click([&]()
        {
            nana::filebox fb{*this, true};
            fb.title("Select the original GISAID CoV fasta file");
            fb.add_filter("fasta file", "*.fasta");
            const auto&files = fb.show();
            if (!files.empty())
            input_file.reset(files[0].string());
        });

        auto& p = get_place();
        p.div(R"(<width=100 gene><vertical  <height=30 input arrange=[variable,60]> 
                                            <height=30 file >
                                            <height=30 split arrange=[40] gap=10> 
                                            margin=10
                                            min=300>   
            )");
        p["gene"]  << gene ;
        p["input"] << input_file_label << set ;
        p["file"]  << input_file ;
        p["split"] << run_split ;
        p.collocate();
    };
};

int main (int argc, char ** argv)
{
    seqan3::argument_parser arg_parser{"Split-CoV-fasta", argc, argv}; 
    
    GUI fm;
    fm.show();
    nana::exec();
    return 0;
}