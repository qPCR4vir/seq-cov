#include <string>
#include <iostream>
#include <vector>
#include <filesystem>
//#include <algorithm>
#include <unordered_map>
#include <ranges>
#include <stdexcept>
	
// #include <seqan3/std/ranges>                    // include all of the standard library's views
#include <seqan3/alphabet/views/all.hpp>        // include all of SeqAn's views 
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>
//#include <seqan3/alphabet/range/hash.hpp>
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/argument_parser/all.hpp> 
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
#include <nana/gui/widgets/spinbox.hpp>
#include <nana/gui/widgets/panel.hpp>
#include <nana/gui/filebox.hpp>
#include <nana/gui/msgbox.hpp>
#include <nana/gui/tooltip.hpp>
#include <nana/gui/place.hpp>

/*
struct dna_deg : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna15; // instead of dna5
    template <typename alph>
    using sequence_container = seqan3::dna15_vector;
};
*/

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
    int beg{0}, end{0};  // position of the actual target IN THIS grouped seq.
    std::string id;
};

class SplitCoVfasta
{
public:    
    using sequence_type = seqan3::dna5_vector;
    using sequence_file_output = decltype(seqan3::sequence_file_output{"name.fasta"});
    
private:
    struct SplitGene
    {
        SplitCoVfasta const &parent;
        bool                split, 
                            group,  // we are asked to split, group or ignore this?
                            ignore{!(split || group)};
        const std::string   gene;

        static constexpr int notfound = std::numeric_limits<int>::lowest(); 
        
        sequence_type       forw, rev, target; 
        int                 fw_match, rv_match;
        int                 beg{0}, end{0}, len{0}, count{0}; 
        const std::string   start{gene+"|"};
        using grouped_by_seq = std::unordered_map<sequence_type, SeqGr>;
        grouped_by_seq                                  grouped; 
        std::unordered_map<std::string, grouped_by_seq> daily,  
                                                        monthly;
        sequence_file_output file_fasta_split{(parent.dir / (gene + "." + parent.fasta_name)
                                              ).replace_extension("fasta")};   
        
        SplitGene(SplitCoVfasta const &parent, std::string gene, bool split, bool group)
        : parent{parent}, 
          gene{gene}, 
          split{split}, group{group}
        {
            //if (split) file_E.open( );  // todo implement conditional split
            std::cout << std::boolalpha  
                << "\nGoing to split: " << (parent.dir / parent.fasta_name).string() 
                << ", gene " << start << gene <<": " << split << " into "
                << (parent.dir / (gene + "." + parent.fasta_name)).string();
            
        }
        
        /// record identified and ...?
        bool check(auto& record)
        {
            // seqan3::debug_stream << '\n' << record.id();
            if (ignore) return false;
            if (!record.id().starts_with(start)) return false;
            
            if (group) check_rec(record);

            if (split) file_fasta_split.push_back(std::move(record));
            
            return true;
        }

        bool check_rec(auto& record)
        {
            sequence_type &sq = record.sequence();
            count++;
            int flank = parent.flank;
            // >Gene name|Isolate name|YYYY-MM-DD|Isolate ID|Passage details/history|Type^^
            //           1            212345678910
            std::string id = std::move(record.id());

            std::size_t date_beg = id.find('|', 1 + id.find('|')) + 1;       // pos next to 2nd |
            std::string date     = id.substr(date_beg, 10);
            std::string month    = id.substr(date_beg,  7);

            if (!len)
            {   
                // assume the first sequence is OK
                // todo don't assume the first sequence is OK

                SeqGr sg = set_seq_pos(sq); 
                sg.id    = std::move(id);
                sg.count = 1;
                len = sg.end - sg.beg;

                end = std::min<int>(sg.end + flank , sq.size()); // seq not too long
                beg = std::max<int>(0, sg.beg - flank);
                    
                sg.beg = sg.beg - beg;
                sg.end = sg.end - beg;

                auto bg = sq.begin()+beg;
                auto en = sq.begin()+end;
                target = sequence_type{bg, en};

                grouped                  [target] = sg ;
                daily   [std::move(date)][target] = sg ;  // end of date
                monthly[std::move(month)][target] = sg ;  // end of month 

                // todo what if incomplete ??
                // todo make sure this first target-sequence is correct !!!!
                // this will be the "standart" target-sequence
                //seqan3::debug_stream << " Target sequence: " << start << target 
                //                     << " (" << beg << ", "  << end  <<")\n" ;
                return true; 
            }

            int lend= std::min<int>(end, sq.size());
            int lbeg= std::max<int>(  0, std::min(lend - len, beg)) ;

            auto bg = sq.begin()+lbeg;
            auto en = sq.begin()+lend;
            auto cur_seq = seqan3::dna5_vector{bg, en};

            SeqGr& sg1 = grouped[cur_seq];  // first try
            if (sg1.count++)  // duplicate seq. More than 99% of cases.
            {
                SeqGr& sgd = daily   [date][cur_seq];            
                if (!sgd.count++) // but it is new in daily
                {
                    sgd.id = record.id();
                    sgd.beg = sg1.beg ;
                    sgd.end = sg1.end ;
                }
                SeqGr& sgm = monthly[month][cur_seq];            
                if (!sgm.count++)  // but it is new in monthly
                {
                    sgm.id = record.id();
                    sgm.beg = sg1.beg ;
                    sgm.end = sg1.end ;
                }
                return true;
            }

            // new, unknown seq if limited as fisrt. sg1 is a new, blank group already in grouped

            SeqGr sg=set_seq_pos(sq);   // true align to correct position - find limits
            
            if (sg.beg == notfound || sg.end == notfound)
            {
                // failed, is new in grouped, but we report it for completeneds of staistic ??!
                sg1.id = record.id();
                sg1.beg = sg.beg ;
                sg1.end = sg.end ;  // sg1.count already set to 1

                SeqGr& sgd = daily   [date][cur_seq];            
                if (!sg1.count++) sgd = sg1; 

                SeqGr& sgm = monthly[month][cur_seq];            
                if (!sgm.count++) sgm = sg1; 

                //seqan3::debug_stream  << start <<  seqan3::dna5_vector{bg, en} << " -- Failed! " 
                //                      << " (" << lbeg << ", " << lend  <<")\n" ;
                return true; // todo ???????????????
            }
            //seqan3::debug_stream << start <<  seqan3::dna5_vector{bg, en} 
            //                     << " (" << lbeg << ", " << lend  <<")\n" ;

            int nlend= std::min<int>(sg.end + flank, sq.size() );
            int nlbeg= std::max<int>(sg.beg - flank, 0);
            
            if (lbeg == nlbeg && lend == nlbeg)  // new seq in the same pos
            {
                sg1.id = id;
                sg1.beg = sg.beg - nlbeg;
                sg1.end = sg.end - nlbeg;

                daily   [date][cur_seq] = sg1; 
                monthly[month][cur_seq] = sg1; 

                return true;
            }
                
            // new seq in a new pos
            grouped.erase(cur_seq);
            //seqan3::debug_stream << " New try ----- (" << lbeg << ", " << lend  <<")\n" ;

            auto new_seq = sequence_type{sq.begin()+nlbeg, sq.begin()+nlend};
            SeqGr&   sgr = grouped[new_seq];

            if (sgr.count++)                                // duplicate seq. More than 99% of cases.
            {
                SeqGr& sgd = daily   [date][cur_seq];            
                if (!sgd.count++) // but it is new in daily
                {
                    sgd.id = id;
                    sgd.beg = sgr.beg ;
                    sgd.end = sgr.end ;
                }
                SeqGr& sgm = monthly[month][cur_seq];            
                if (!sgm.count++)  // but it is new in monthly
                {
                    sgm.id = id;
                    sgm.beg = sgr.beg ;
                    sgm.end = sgr.end ;
                }

                return true;
            }
            //seqan3::debug_stream << " It was new !!!!!!!\n" ;
            sgr.beg = sg.beg - nlbeg;
            sgr.end = sg.end - nlbeg;
            sgr.id = id;
            daily   [date][new_seq] = sgr; 
            monthly[month][new_seq] = sgr; 
            return true;
        }
        SeqGr set_seq_pos(const sequence_type& s)
        {
            // new, unknown seq. We need to find the right position of the target sequence
            // less than 1% of the seq. May be around 50k?
            SeqGr sg;
            
            auto output_config = seqan3::align_cfg::output_score{}          | 
                                 seqan3::align_cfg::output_begin_position{} |
                                 seqan3::align_cfg::output_end_position{}   | 
                                 seqan3::align_cfg::output_alignment{};
            auto method = seqan3::align_cfg::method_local{};
            seqan3::align_cfg::scoring_scheme     scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, 
                                                         seqan3::mismatch_score{-1}}};
            seqan3::align_cfg::gap_cost_affine gap_costs{seqan3::align_cfg::open_score{0}, 
                                                         seqan3::align_cfg::extension_score{-1}};
            auto config = method | scheme | gap_costs | output_config;

            // try the fw primer
            auto results = seqan3::align_pairwise(std::tie(s, forw), config);
            auto & res = *results.begin();
            //seqan3::debug_stream << "Alignment: " << res.alignment() << " Score: "     << res.score() ;
            //seqan3::debug_stream << ", Target: (" << res.sequence1_begin_position() << "," << res.sequence1_end_position() << ")";
            //seqan3::debug_stream << ", Primer: (" << res.sequence2_begin_position() << "," << res.sequence2_end_position() << "). " ;
            
            if (res.score() > fw_match)  // forw primer found. 
            {
                sg.beg = res.sequence1_begin_position() - res.sequence2_begin_position() ;
                if (len)  // not the first time/seq - the "ref." seq was already set
                {
                    sg.end = sg.beg + len;
                    return sg; // only find fw pr pos. Don't count for insertions/deletions or bad rev.
                }
                else if (sg.beg < 0)
                 throw std::runtime_error{"First " + gene + " sequence don't contain "
                   "the full forward primer. Score: " + std::to_string(res.score()) +
                   " that begin at position "  + std::to_string(res.sequence2_begin_position())  };
            }
            else 
            {
                if (!len)
                    throw std::runtime_error{"First " + gene + " sequence don't contain "
                                             "the forward primer. Score: " + std::to_string(res.score() )};
                else
                    sg.beg = notfound;  // mark beg pos as not found !!
            }
            // we need to find rev primer location

            auto res_rev = seqan3::align_pairwise(std::tie(s, rev), config);
            auto & res_r = *res_rev.begin();
            //seqan3::debug_stream << "\nAlignment: " << res_r.alignment() << " Score: " << res_r.score() ;
            //seqan3::debug_stream << ", Target: ("     << res_r.sequence1_begin_position() << "," << res_r.sequence1_end_position() << ")";
            //seqan3::debug_stream << ", rev Primer: (" << res_r.sequence2_begin_position() << "," << res_r.sequence2_end_position() << "). ";

            if (res_r.score() > rv_match)  // rev primer found. 
            {   
                sg.end = res_r.sequence1_end_position() + (rev.size() - res_r.sequence2_end_position());
                if (len)
                    sg.beg = sg.end - len;
                else if (sg.end > s.size())
                 throw std::runtime_error{"First " + gene + " sequence don't contain "
                   "the full reverse primer. Score: " + std::to_string(res.score()) +
                   " that end at position "  + std::to_string(res.sequence2_end_position())  };
                return sg;
            }

            if (!len)
                    throw std::runtime_error{"First " + gene + " sequence don't contain "
                                              "the reverse primer. Score: " + std::to_string(res.score() )};
            sg.end = notfound;  // mark both beg and end pos as not found !! Better try to align whole seq?? not sure..
            return sg;
        }

        void write_grouped ()
        {
            
            using types = seqan3::type_list<std::vector<seqan3::dna5>, std::string>;
            using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
            using record_t = seqan3::sequence_record<types, fields>;
            using sgr_t = decltype(grouped)::value_type;
            using sgr_p = sgr_t*;

            if (!group) return;  // todo ?????

            std::filesystem::path gr = parent.dir / (gene + ".grouped-" + parent.fasta_name);
            std::filesystem::path grd= parent.dir / (gene + ".daily-"   + parent.fasta_name);
            std::filesystem::path grm= parent.dir / (gene + ".monthly-" + parent.fasta_name);

            seqan3::debug_stream << "Going to write: " << gr.string() << "\n" ;
            seqan3::sequence_file_output file_e_gr {gr },
                                         file_e_grd{grd},
                                         file_e_grm{grm};
            
            std::vector<sgr_p> gr_v;
            gr_v.reserve(grouped.size());
            for (sgr_t& sgr : grouped) gr_v.push_back(&sgr);

            std::sort(gr_v.begin(), gr_v.end(), []( sgr_p &a,  sgr_p &b)
            {return a->second.count > b->second.count;});

            for(sgr_p sg:gr_v)
            {
                auto id = "x_" + std::to_string(sg->second.count) 
                         + "_" + std::to_string(sg->second.beg) 
                         + "_" + std::to_string(sg->second.end)
                         + "_" + sg->second.id ;
                file_e_gr.push_back(record_t{std::move(sg->first), std::move(id)});
            }

            for(auto& [date, group]: daily)
                for(auto &[seq, gr] : group)
                {
                    auto id = "d_" + date + 
                         +"_x_"+ std::to_string(gr.count) 
                         + "_" + std::to_string(gr.beg) 
                         + "_" + std::to_string(gr.end)
                         + "_" + gr.id ;
                    file_e_grd.push_back(record_t{std::move(seq), std::move(id)});
                }

            for(auto& [date, group]: monthly)
                for(auto &[seq, gr] : group)
                {
                    auto id = "d_" + date + 
                         +"_x_"+ std::to_string(gr.count) 
                         + "_" + std::to_string(gr.beg) 
                         + "_" + std::to_string(gr.end)
                         + "_" + gr.id ;
                    file_e_grm.push_back(record_t{std::move(seq), std::move(id)});
                }

        }
        
        SplitGene& set_forw(const std::string& pr)
        {
            if (pr.empty()) return *this;
            auto e = pr | seqan3::views::char_to<seqan3::dna5> ;
            forw = seqan3::dna5_vector{e.begin(), e.end()}; 
            fw_match = std::round(parent.match * double(forw.size()) / 100.0);
            seqan3::debug_stream << "\n From forward " << pr 
                                 << " = [ " << forw << " ] with minimum match:" << fw_match;
            return *this;
        }
        SplitGene& set_forw(seqan3::dna5_vector forw)
        {
            if (forw.empty()) return *this;
            fw_match = std::round(parent.match * double(forw.size()) / 100.0);
            seqan3::debug_stream << "\n From forward = [ " << forw << " ] with minimum match:" << fw_match;
            return *this;
        }
        SplitGene& set_rev(const std::string& pr)
        {
            if (pr.empty()) return *this;
            auto e = pr | seqan3::views::char_to<seqan3::dna5> | std::views::reverse | seqan3::views::complement ; 
            rev = seqan3::dna5_vector{e.begin(), e.end()}; 
            rv_match = std::round(parent.match * double(rev.size()) / 100.0);
            seqan3::debug_stream  << "\n From reverse " << pr << " = [ "  << rev << " ] with minimum match:" << rv_match;
            return *this;
        }
        SplitGene& set_rev(seqan3::dna5_vector rev)
        {
            if (rev.empty()) return *this;
            rv_match = std::round(parent.match * double(rev.size()) / 100.0);
            seqan3::debug_stream  << " From reverse = [ "  << rev << " ] with minimum match:" << rv_match;
            return *this;
        }        
    };    
public:
    std::filesystem::path fasta;
    std::filesystem::path dir       {fasta.parent_path()};
    std::string           fasta_name{fasta.filename().string()},
                          from, to;
    int                   flank;
    double                match;
    bool                  check_date = !(from.empty() && to.empty());

private:
    std::vector<SplitGene> genes;

public:
    SplitCoVfasta(const std::filesystem::path& fasta, int flank, double match, std::string from, std::string to)
    : fasta{fasta}, flank{flank}, match{match}, from{from}, to{to}
    {}

    void add_gene(std::string gene, bool split, bool group, std::string fw="", std::string rv="")  // todo implement conditional split
    {
        genes.emplace_back(*this, gene, split, group).set_forw(fw)
                                                     .set_rev (rv);
    }

    void split_fasta( )
    {
        
        // Initialise a file input object with a FASTA file.
        seqan3::sequence_file_input file_in{fasta};

        long t{0L}, m{(1L<<18)-1};
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
            //if (t>10000000) break;
        }
        seqan3::debug_stream << "\nTotal= " << t  << "\n" ; 
        for (auto & gene : genes)
        {
            seqan3::debug_stream << gene.gene <<"= " << gene.count 
                                 << ". Grouped: "    << gene.grouped.size() << ". " ; 
            gene.write_grouped();
        }

    }
};
 
class GeneGUI: public nana::panel<false>
{
    public:
    GeneGUI(nana::window parent, std::string gene_name="", std::string fw="", std::string rv="")
    : nana::panel<false>(parent)
    {
        gene.multi_lines(false).reset(gene_name);
        forw.tip_string("Seq. forward primer").multi_lines(false).reset(fw);
        rev.tip_string("Seq. reverse primer").multi_lines(false).reset(rv);
        input_file.tip_string("fasta file with primers").multi_lines(false);
        set.events().click([&]()
        {
            nana::filebox fb{*this, true};
            fb.title("Select a fasta file with primer sequences for gene " + this->gene.text());
            fb.add_filter("fasta file", "*.fasta");
            const auto&files = fb.show();
            if (!files.empty())
            {
                 this->input_file.reset(files[0].string());
                 seqan3::debug_stream << "\nGoing to load: " << files[0].string();
                 seqan3::sequence_file_input file_in{files[0]};
                 std::string fw, rv;
                 for (auto & primer : file_in)
                 {
                    if (fw.empty())
                    {
                        seqan3::debug_stream << "\nGoing to convert fw: \n>" 
                                             << primer.id() << "\n" << primer.sequence();
                        auto e = primer.sequence() | seqan3::views::to_char;
                        fw = std::string{e.begin(), e.end()}; 
                        this->forw.reset(fw);
                        continue;
                    }
                    seqan3::debug_stream << "\nGoing to convert rv: \n>" 
                                         << primer.id() << "\n" << primer.sequence();
                    auto e = primer.sequence() | seqan3::views::to_char;
                    rv = std::string{e.begin(), e.end()}; 
                    this->rev.reset(rv);
                    break;

                 }
            }  
        });

        plc.div(R"( <min=500 all arrange=[35, 60,80,100,180,180,30,200] gap=10> )");
        plc["all"] << "Gene" << gene << split << group << forw << rev << set << input_file;
        plc.collocate();
    }
    nana::textbox  gene     {*this};  // todo use a combox 
    nana::checkbox split    {*this, "split full"},
                   group    {*this, "group target"};
    nana::textbox  forw     {*this};
    nana::textbox  rev      {*this};
    nana::button   set      {*this, "..."}; 
    nana::textbox  input_file{*this};
    nana::place    plc      {*this};
};
class GUI: public nana::form
{
    nana::label    input_file_label{*this, "Original fasta file:"};
    nana::textbox  input_file      {*this},
                   from            {*this},
                   to              {*this};
    nana::spinbox  flank           {*this},
                   match           {*this},
                   period          {*this};
    nana::button   set             {*this, "&Select"}, 
                   run_split       {*this, "S&plit" };
    GeneGUI        E               {*this, "E"},
                   N               {*this, "N"},
                   S               {*this, "Spike"},
                   R               {*this, "NSP12"};
    
    
public:
    GUI() : nana::form{nana::api::make_center(1100, 350)}
    {
        caption("Split-CoV-fasta. v1.00.00");

        input_file.tip_string("Original fasta file:").multi_lines(false);
        flank.range(0, 100, 1);
        flank.value("20");
        match.range(50.0, 100.0, 1.0);
        match.value("70.0");
        period.range(0, 12, 1);
        period.value("3");

        E.split.check(false);  E.group.check(true);
        N.split.check(false);  N.group.check(true);
        S.split.check(false);  S.group.check(false);  // R

        run_split.events().click([&]()
        {
            std::filesystem::path fasta{input_file.text()};
            if (!std::filesystem::is_regular_file(fasta)) return;  // todo msg

            SplitCoVfasta sp{fasta, flank.to_int(), match.to_double(), from.text(), to.text()};

            auto Add_Gene = [&sp](auto &g)
            {
            if (g.split.checked() || g.group.checked()) 
                sp.add_gene(g.gene.text(),     
                            g.split.checked(), g.group.checked(), 
                            g.forw.text(),     g.rev.text());
            };

            // todo implement conditional split
            Add_Gene(E);
            Add_Gene(N);
            Add_Gene(S);  
            Add_Gene(R);

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
        p.div(R"(<vertical  margin=10 gap=10 min=350
                    <height=25 input arrange=[variable, 35,50,  75,50,  60, 60] gap=10> 
                    <height=10  >
                    <height=30 file >
                    <height=10  >
                    <height=133 vertical genes gap=10> 
                    <height=10  >
                    <height=30 no_dates arrange=[ 35,80,  35,80,  45,40,180] gap=10> 
                  >   
            )");
        p["input"] << input_file_label << "Flank:"      << flank 
                                       << "Min match %" << match << set  << run_split ;
        p["file"]  << input_file ;
        p["genes"] << E << N << S << R;
        p["dates"] << "From date: " << from << " To date: " << to
                   << "Separate by: " << period << " months \n(0 - all time in 1 period).";
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