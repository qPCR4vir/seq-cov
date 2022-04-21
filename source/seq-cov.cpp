#include <string>
#include <filesystem>
#include <unordered_map>
//#include <ranges>
 
#include <seqan3/std/ranges>                    // include all of the standard library's views
#include <seqan3/alphabet/views/all.hpp>        // include all of SeqAn's views 
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>
#include <seqan3/argument_parser/all.hpp> 
//#include <seqan3/alphabet/range/hash.hpp>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/std/ranges>

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
namespace std
{
/*!\brief Struct for hashing a range of characters.
 * \ingroup alphabet_range
 * \tparam urng_t The type of the range; Must model std::ranges::input_range and the reference type of the range of the
 *                range must model seqan3::semialphabet.
 * \details
 * \experimentalapi{Experimental since version 3.1.}
 */
template <ranges::input_range urng_t>
//!\cond
    requires seqan3::semialphabet<std::ranges::range_reference_t<urng_t>>
//!\endcond
struct hash<urng_t>
{
    /*!\brief Compute the hash for a range of characters.
     * \tparam urng2_t  The same as `urng_t` (+- cvref); used to get forwarding reference in the interface.
     * \param[in] range The input range to process. Must model std::ranges::input_range and the reference type of the
     *                  range of the range must model seqan3::semialphabet.
     * \returns size_t.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <ranges::input_range urng2_t>
    //!\cond
        requires seqan3::semialphabet<std::ranges::range_reference_t<urng2_t>>
    //!\endcond
    size_t operator()(const urng2_t & range) const noexcept
    {
        using alphabet_t = std::ranges::range_reference_t<urng_t>;
        size_t result{0};
        hash<alphabet_t> h{};
        for (auto character : range)
        {
            result *= seqan3::alphabet_size<alphabet_t>;
            result += h(character);
        }
        return result;
    }
};

} // namespace std
class SplitCoVfasta
{
    std::filesystem::path fasta, dir;
    std::string fasta_name;
    bool split_E, split_N, split_S;
    seqan3::dna5_vector e_forw, e_rev, n_forw, n_rev, s_forw, s_rev;    

public:
    SplitCoVfasta(const std::filesystem::path& fasta, bool split_E, bool split_N, bool split_S)
    : fasta{fasta}, split_E{split_E}, split_N{split_N}, split_S{split_S}, dir{fasta.parent_path()}, fasta_name{fasta.filename().string()}
    {
    std::cout << "\nGoing to split: " << fasta.string() << ", gene E: " << split_E  
              << ", gene N: " << split_N << ", gene Spike: " << split_S   ;
    
    }
    void split_fasta( )
    {

        // Initialise a file input object with a FASTA file.
        seqan3::sequence_file_input file_in{fasta};

        long e{0L}, n{0L}, s{0L}, t{0L}, m{(1L<<15)-1};
        seqan3::debug_stream << "m= " << m << "\n" ; 

        using seq_t = decltype(file_in)::sequence_type;
        using id_t = decltype(file_in)::id_type;
        std::unordered_multimap<seq_t, id_t> e_grouped;
        {
            std::filesystem::path e_fasta = split_E ? dir / ("E." + fasta_name) : std::filesystem::path( "E.fasta");
            std::filesystem::path n_fasta = split_N ? dir / ("N." + fasta_name) : std::filesystem::path( "N.fasta");
            std::filesystem::path s_fasta = split_S ? dir / ("Spike." + fasta_name) : std::filesystem::path( "S.fasta");
            seqan3::sequence_file_output file_E{e_fasta};
            seqan3::sequence_file_output file_N{n_fasta};
            seqan3::sequence_file_output file_S{s_fasta};

            for (auto & record : file_in)
            {
                    if (split_E && record.id().starts_with("E|"))        
                    {   
                        e_grouped.emplace(std::make_pair(record.sequence(), record.id()));
                        //seqan3::debug_stream << record.sequence() << "\n" ;
                        file_E.push_back(std::move(record));++e;
                            
                    }
                else if (split_N && record.id().starts_with("N|"))        {file_N.push_back(std::move(record));++n;}
                else if (split_S && record.id().starts_with("Spike|"))    {file_S.push_back(std::move(record));++s;}
                if (!(++t & m)) 
                    seqan3::debug_stream << "T= " << t << ", N= " << n << ", E= " << e << ", Spike= " << s << "\n" ; 
                if (t>100000) break;
            }
        }
        seqan3::debug_stream << "Total= " << t  << ". N= " << n << ", E= " << e << ", Spike= " << s << ". Grouped: " << e_grouped.size() << " sequences into " << e_grouped.bucket_count() << " groups\n" ; 
        std::filesystem::path e_gr = split_E ? dir / ("E.gr" + fasta_name) : std::filesystem::path( "E.fasta");
        seqan3::sequence_file_output file_e_gr{e_gr};

    }
    void set_e_forw(const std::string& pr)
    {
         e_forw.clear();
         auto e = pr | seqan3::views::char_to<seqan3::dna5> ;
         e_forw = seqan3::dna5_vector{e.begin(), e.end()}; 
         seqan3::debug_stream << "\n from " << pr << " e_forw: " << e_forw;
    }
    void set_e_rev(const std::string& pr)
    {
         auto e = pr | seqan3::views::char_to<seqan3::dna5> | std::views::reverse | seqan3::views::complement ; 
         e_rev = seqan3::dna5_vector{e.begin(), e.end()}; 
         seqan3::debug_stream  << " from " << pr << ", e_rev: " << e_rev;
    }
};
 

class GUI: public nana::form
{
    nana::label    input_file_label{*this, "Original fasta file:"};
    nana::textbox  input_file      {*this};
    nana::button   set             {*this, "&Select"}, 
                   run_split       {*this, "S&plit" };
    nana::group    gene            {*this, "Gene"   };
    nana::textbox  e_forw_tb       {*this};
    nana::textbox  e_rev_tb        {*this};
    nana::textbox  n_forw_tb       {*this};
    nana::textbox  n_rev_tb        {*this};    
    nana::textbox  s_forw_tb       {*this};
    nana::textbox  s_rev_tb        {*this};
public:
    GUI() : nana::form{nana::api::make_center(1000, 200)}
    {
        caption("Split-CoV-fasta. v0.01.00");

        input_file.tip_string("Original fasta file:").multi_lines(false);
        e_forw_tb.tip_string("Seq. forward primer").multi_lines(false).reset("ACAggTACgTTAATAgTTAATAgCgT");
        e_rev_tb.tip_string("Seq. reverse primer").multi_lines(false).reset("CAATATTgCAgCAgTACgCACA");
        n_forw_tb.tip_string("Seq. forward primer").multi_lines(false);
        n_rev_tb.tip_string("Seq. reverse primer").multi_lines(false);
        s_forw_tb.tip_string("Seq. forward primer").multi_lines(false);
        s_rev_tb.tip_string("Seq. reverse primer").multi_lines(false); 

        auto& E = gene.add_option("E");
        auto& N = gene.add_option("N");
        auto& S = gene.add_option("Spike");

        run_split.events().click([&]()
        {
            SplitCoVfasta sp{input_file.text(), E.checked(), N.checked(), S.checked()};
            sp.set_e_forw(e_forw_tb.text());
            sp.set_e_rev(e_rev_tb.text());
            sp.split_fasta();
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
        p.div(R"(<width=100 gene><vertical  margin=10 gap=10 min=300
                                            <height=30 input arrange=[variable,60, 60] gap=10> 
                                            <height=30 file >
                                            <height=30 e_pr arrange=[40,180,180] gap=10> 
                                            <height=30 n_pr arrange=[40,180,180] gap=10> 
                                            <height=30 s_pr arrange=[40,180,180] gap=10> 
                                            >   
            )");
        p["gene"]  << gene ;
        p["input"] << input_file_label << set  << run_split ;
        p["file"]  << input_file ;
        p["e_pr"] << "E-gene" << e_forw_tb <<  e_rev_tb;
        p["n_pr"] << "N-gene" << n_forw_tb <<  n_rev_tb;
        p["s_pr"] << "Spike" << s_forw_tb <<  s_rev_tb;
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