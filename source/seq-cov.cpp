#include <string>
#include <iostream>
#include <vector>
#include <filesystem>
//#include <algorithm>
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
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>

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

struct SeqGr
{
    int count{0};
    int beg, end;  // position of the actual target IN THIS grouped seq.
    std::string id;
};

class SplitCoVfasta
{
public:    
    using sequence_type = seqan3::dna5_vector;
private:
    struct SplitGene
    {
        SplitCoVfasta const &parent;
        bool          split;
        std::sting    gene;
        
        sequence_type forw, rev, target; 
        int           beg{0}, end{0}, len{0}, count{0}; 
        std::sting    start;
        std::unordered_map<sequence_type, SeqGr> grouped; 
        seqan3::sequence_file_output file_E;   // todo implement conditional split
        
        SplitGene(SplitCoVfasta const *parent, std::string gene, bool split)
        : parent{*parent}, 
          gene{gene}, 
          split{split}, 
          start{gene+"|"},
          file_fasta_split{parent.dir / (gene + "." + parent.fasta_name)};
        {
            //if (split) file_E.open( );  // todo implement conditional split
            std::cout << std::boolalpha  
                << "\nGoing to split: " << (parent.dir / parent.fasta_name).string() 
                << ", gene " << gene <<": " << split << " into "
                << (parent.dir / (gene + "." + parent.fasta_name)).string();
        }
        
        bool check(auto& record)
        {
            if (!split_e.check_rec(record)) 
                return false;

            file_fasta_split.push_back(std::move(record));
            return true;
        }

        bool check_rec(auto& record)
        {
            if (!split_e.split) return false;
            if (!record.id().starts_with(start)) return false;
            const seqan3::dna5_vector &sq = record.sequence();
            count++;
            
            if (!len)
            {   
                // assume the first sequence is OK
                // todo don't assume the first sequence is OK
                SeqGr sg=set_seq_pos(sq); 
                sg.id = record.id();
                sg.count = 1;
                len = sg.end - sg.beg;

                if (sg.end + 20 > sq.size())  // todo set this 20 as program parameter
                    end  = sq.size();     // seq not too long
                else 
                    end = sg.end + 20 ;

                if (sg.beg < 20)  // todo set this 20 as program parameter
                    beg = 0;
                else 
                {
                    beg = beg - 20 ;
                    sg.beg = 20;
                }
                sg.end = sg.end - sg.beg;

                bg = sq.begin()+beg;
                en = sq.begin()+end;
                target = seqan3::dna5_vector{bg, en};

                grouped[target] = sg;
                // todo what if incomplete ??
                // todo make sure this first target-sequence is correct !!!!
                // this will be the "standart" target-sequence
                
                seqan3::debug_stream << "\nTarget sequence: " << start << target << "\n" ; 
                return true; 
            }

            int lend= end > sq.size() ? sq.size() : end;
            int lbeg= lend - len < beg ? (lend - len < 0 ? 0 : lend - len) : beg ;

            auto bg = sq.begin()+lbeg;
            auto en = sq.begin()+lend;
            SeqGr& sg = grouped[seqan3::dna5_vector{bg, en}];  // first try
            if (sg.count++) return true; // duplicate seq. More than 99% of cases.
            
            // new, unknown seq.
            seqan3::debug_stream << "\n"  <<  seqan3::dna5_vector{bg, en} << "\n" ;

            sg=set_seq_pos(sq);   // true align to correct position 
            sg.id = record.id();
            
            if (sg.end + 20 > sq.size())  // todo set this 20 as program parameter
                lend = sg.end + 20 ;

            if (sg.beg < 20)  // todo set this 20 as program parameter
                lbeg = 0;
            else 
            {
                lbeg = sg.beg - 20 ;
                sg.beg = 20;
            }
            sg.end = sg.end - sg.beg;

            if (beg <= lbeg && lend <= end) return true; 

            grouped.erase(seqan3::dna5_vector{bg, en});

            SeqGr& sgr = grouped[seqan3::dna5_vector{sq.begin()+lbeg, sq.begin()+lend}];
            if (sgr.count++) return true; // duplicate seq. More than 99% of cases.

            sgr = sg;
            return true;
        }
        SeqGr set_seq_pos(const seqan3::dna5_vector& s)
        {
            // new, unknown seq. We need to find the right position of the target sequence
            // less than 1% of the seq. May be around 50k?
            SeqGr sg;
            
            auto output_config = seqan3::align_cfg::output_score{}        | seqan3::align_cfg::output_begin_position{} |
                                 seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};
            auto method = seqan3::align_cfg::method_local{};
            seqan3::align_cfg::scoring_scheme scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};
            seqan3::align_cfg::gap_cost_affine gap_costs{seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-1}};
            auto config = method | scheme | gap_costs | output_config;

            /*
            auto config = seqan3::align_cfg::method_global{
                    seqan3::align_cfg::free_end_gaps_sequence1_leading{false},    // target seq
                    seqan3::align_cfg::free_end_gaps_sequence2_leading{true},     // primer
                    seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},   // target seq
                    seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};   // primer
            */
            
            auto results = seqan3::align_pairwise(std::tie(s, forw), config);
            auto & res = *results.begin();
            seqan3::debug_stream << "Score: "     << res.score() ;
            seqan3::debug_stream << ", Target: (" << res.sequence1_begin_position() << "," << res.sequence1_end_position() << ")";
            seqan3::debug_stream << ", Primer: (" << res.sequence2_begin_position() << "," << res.sequence2_end_position() << ")\n"
                                 << "Alignment: " << res.alignment();
            
            if (res.score() > 15)  // forw primer found. todo set this 15 as program parameter
            {
                sg.beg = res.sequence1_begin_position();   // todo check beg primer align
                if (len)
                {
                    sg.end = sg.beg + len;
                    return sg; 
                }
            }
            else 
            {
                if (!len)
                    throw std::exception{"First " + gene + " sequence don't contain the forward primer"};
                else
                    sg.beg = -1;  // todo ??
            }
            // we need to find rev primer location

            auto res_rev = seqan3::align_pairwise(std::tie(s, rev), config);
            auto & res_r = *res_rev.begin();
            seqan3::debug_stream << "Score: " << res_r.score() ;
            seqan3::debug_stream << ", Target: ("     << res_r.sequence1_begin_position() << "," << res_r.sequence1_end_position() << ")";
            seqan3::debug_stream << ", rev Primer: (" << res_r.sequence2_begin_position() << "," << res_r.sequence2_end_position() << ")\n"
                                 << "Alignment: " << res_r.alignment();

            if (res_r.score() > 15)  // rev primer found. todo set this 15 as program parameter
            {   
                sg.end = res_r.sequence1_end_position() ;  // todo check end primer align
                if (len)
                    sg.beg = sg.end - len;
                return sg;
            }

            if (!len)
                    throw std::exception{"First " + gene + " sequence don't contain the reverse primer"};

            throw std::exception{gene + " sequence don't contain the forward and reverse primers"};
        }

        void write_grouped (std::filesystem::path& const dir, std::string& const fasta_name)
        {
            if (!split) return;
            std::filesystem::path gr = dir / ("E_gr." + fasta_name);
            seqan3::sequence_file_output file_e_gr{gr};

        }
        
        SplitGene& set_forw(const std::string& pr)
        {
            if (pr.empty()) return *this;
            auto e = pr | seqan3::views::char_to<seqan3::dna5> ;
            forw = seqan3::dna5_vector{e.begin(), e.end()}; 
            seqan3::debug_stream << "\n from " << pr << ", forw: " << forw;
            return *this;
        }
        void set_rev(const std::string& pr)
        {
            if (pr.empty()) return *this;
            auto e = pr | seqan3::views::char_to<seqan3::dna5> | std::views::reverse | seqan3::views::complement ; 
            rev = seqan3::dna5_vector{e.begin(), e.end()}; 
            seqan3::debug_stream  << " from " << pr << ", rev: " << rev;
            return *this;
        }
    };    
public:
    std::filesystem::path fasta;
    std::filesystem::path dir       {fasta.parent_path()};
    std::string           fasta_name{fasta.filename().string()};
private:
    std::vector<SplitGene> genes;

public:
    SplitCoVfasta(const std::filesystem::path& fasta)
    : fasta{fasta}
    {}

    void add_gene(std::string gene, bool split, std::string fw="", std::string rv="")  // todo implement conditional split
    {
        genes.emplace_back(this, gene, split).set_forw(fw)
                                             .set_rev (rv);
    }
    void split_fasta( )
    {
        // Initialise a file input object with a FASTA file.
        seqan3::sequence_file_input file_in{fasta};

        long t{0L}, m{(1L<<15)-1};
        seqan3::debug_stream << "\nchunk - m= " << m << "\n" ; 

        for (auto & record : file_in)
        {
            for (auto & gene : genes)
                    if (gene.check(record)) break;

            if (!(++t & m)) 
            {
                seqan3::debug_stream << "\tT= " << t  << "\n" ;
                for (auto & gene : genes)
                   seqan3::debug_stream << gene.gene <<"= " << gene.count 
                                        << ". Grouped: "    << gene.grouped.size() << "\n" ; 
            }
        }
        seqan3::debug_stream << "\nTotal= " << t  << "\n" ;
        for (auto & gene : genes)
            seqan3::debug_stream << gene.gene <<"= " << gene.count 
                                 << ". Grouped: "    << gene.grouped.size() << "\n" ; 
        

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
        n_forw_tb.tip_string("Seq. forward primer").multi_lines(false).reset("CCAAAAggCTTCTACgCAgA");
        n_rev_tb.tip_string("Seq. reverse primer").multi_lines(false).reset("TgCCTggAgTTgAATTTCTTgA");
        s_forw_tb.tip_string("Seq. forward primer").multi_lines(false);
        s_rev_tb.tip_string("Seq. reverse primer").multi_lines(false); 

        auto& E = gene.add_option("E"); E.check(true);
        auto& N = gene.add_option("N"); N.check(true);
        auto& S = gene.add_option("Spike");

        run_split.events().click([&]()
        {
            std::filesystem::path fasta{input_file.text()};
            if (!fasta.is_regular()) return;  // todo msg

            SplitCoVfasta sp{};

            // todo implement conditional split
            if (E.checked()) 
                sp.add_gene("E", E.checked(), e_forw_tb.text(), e_rev_tb.text())
            if (N.checked())  
                sp.add_gene("N", N.checked(), n_forw_tb.text(), n_rev_tb.text())
            if (S.checked()) 
                sp.add_gene("Spike", S.checked(), s_forw_tb.text(), s_rev_tb.text())

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