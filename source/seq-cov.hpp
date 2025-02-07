
// file seq_cov.hpp

#pragma once

#include <string>
#include <vector>
#include <filesystem>
#include <unordered_map>
#include <ranges>

// // Define a fallback value for hardware_destructive_interference_size
// #ifndef __cpp_lib_hardware_interference_size
// constexpr std::size_t hardware_destructive_interference_size = 64; // Typical cache line size
// #else
// constexpr std::size_t hardware_destructive_interference_size = std::hardware_destructive_interference_size;
// #endif

#include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp> // for nucleotide_scoring_scheme
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
//#include <seqan3/alignment/decorator/gap_decorator.hpp>

enum {  debugging_NOT_USED= 0, 
        debugging_NO_DEBUG= 1, 
        debugging_ERROR   = 10, 
        debugging_WARNING = 20, 
        debugging_INFO    = 30, 
        debugging_VERBOSE = 35,
        debugging_DEBUG   = 40, 
        debugging_TRACE   = 80,
        debugging_ALL     = 1000};
constexpr auto debugging = debugging_ALL;  // How much debug info to print
// todo grap most std::cout and seqan3::debug_stream with  if constexpr (debugging) { }

namespace cov
{

/*!\brief Struct for hashing a range of characters.
 * \ingroup alphabet_range
 * \tparam urng_t The type of the range; Must model std::ranges::input_range and the reference type of the range of the
 *                range must model seqan3::semialphabet.
 * \details
 * \experimentalapi{Experimental since version 3.1.}
 */
template <std::ranges::input_range urng_t>
//!\cond
    requires seqan3::semialphabet<std::ranges::range_reference_t<urng_t>>
//!\endcond
struct hash
{
    /*!\brief Compute the hash for a range of characters.
     * \tparam urng2_t  The same as `urng_t` (+- cvref); used to get forwarding reference in the interface.
     * \param[in] range The input range to process. Must model std::ranges::input_range and the reference type of the
     *                  range of the range must model seqan3::semialphabet.
     * \returns size_t.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <std::ranges::input_range urng2_t>
    //!\cond
        requires seqan3::semialphabet<std::ranges::range_reference_t<urng2_t>>
    //!\endcond
    size_t operator()(const urng2_t & range) const noexcept
    {
        using alphabet_t = std::ranges::range_reference_t<urng_t>;
        size_t result{0};
        std::hash<alphabet_t> h{};
        for (auto character : range)
        {
            result *= seqan3::alphabet_size<alphabet_t>;
            result += h(character);
        }
        return result;
    }
};    
using namespace seqan3::literals;

using oligo_seq_alph = seqan3::dna15;
using oligo_seq_t    = seqan3::dna15_vector;
using msa_seq_alph   = seqan3::gapped<oligo_seq_alph>;
using msa_seq_t      = std::vector<msa_seq_alph>;
using sequence_file_output = decltype(seqan3::sequence_file_output{"name.fasta"});
struct MSA : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = msa_seq_alph; // instead of dna5
    using sequence_legal_alphabet = msa_seq_alph;
};
struct OLIGO : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = oligo_seq_alph; // instead of dna5
    
};
struct SeqPos {long beg{0}, end{0};};
struct oligo
{
    oligo_seq_t   seq, ref_seq;
    msa_seq_t     msa_seq, msa_ref;
    std::string   name, code;
    long          ref_beg{0}, ref_end{0}, 
                  msa_beg{0}, msa_end{0};
    int           match{0}, len{0}, msa_len{0};
    bool          reverse{false}, is_probe{false};
};

struct pattern_q
{
    oligo&      primer;
    std::string pattern;
    int         Q{0}, mm{0}, N{0}, gap{0}, crit{0};
};
struct target_q
{
    std::vector<pattern_q> patterns;
    std::string target_pattern;
};
struct parsed_id
{
    std::string country;      // 
    std::string isolate;      // 
    std::string EPI_ISL;      // 
    int         year{0}, month{0}, day{0}; // 

    std::string continent;
    std::string region;
    std::string clade;
    std::string pango;
    std::string pango_version;
};

using Metadata = std::unordered_map<std::string, parsed_id>;  // isolate -> metadata
struct clade_count
{
    parsed_id   *id = nullptr;  // pointer to Metadata parsed_id in SplitCoVfasta::metadata
    int         count = 0;
};
struct country_count
{
    std::unordered_map<std::string, clade_count> clades;
    int         count = 0;
};
struct continent_count
{
    std::unordered_map<std::string, country_count> countries;
    int         count = 0;
};
struct day_count
{
    std::unordered_map<std::string, continent_count> continents;
    int         count = 0;
};
struct month_count
{
    std::unordered_map<int, day_count> days;
    int         count = 0;
};
struct year_count
{
    std::unordered_map<int, month_count> months;
    int count = 0;
};
struct target_count
{
    std::unordered_map<int, cov::year_count> years;
    target_q target;
    int count = 0;
};
enum class GISAID_format {fasta, gene_msa, msa, allprot, allaa, allcodon, allnucprot, allnucaa, allnucprotcodon, allnucprotcodonaa};
class SplitCoVfasta;
struct SplitGene
{
    seqan3::nucleotide_scoring_scheme<int8_t> mismatch; // hamming distance is default

    SplitCoVfasta      &parent;
    const std::string   gene;
    std::vector<oligo>  f_primers, r_primers, probes_s, probes_a;
    std::vector<oligo>  all_oligos;
    int                 extern_forw_idx{}, extern_rev_idx{};  // index of the external primer in f_primers and r_primers respectively

    long          ref_beg{0}, ref_end{0},  // amplicon positions in the reference sequence
                  msa_beg{0}, msa_end{0};  // amplicon positions in the MSA
    int           msa_len{0}, ref_len{0};  // length of the amplicon

    bool read_oligos(const std::filesystem::path& oligos);

    // References and pointers to either key or data stored in the container 
    // are only invalidated by erasing that element, even when the corresponding 
    // iterator is invalidated.
    using grouped_by_msa_seq = std::unordered_map<msa_seq_t  , target_count, hash<msa_seq_t>>;
    using grouped_by_seq     = std::unordered_map<oligo_seq_t, target_count, hash<oligo_seq_t>>;
    grouped_by_msa_seq  msa_grouped; 
    grouped_by_seq      grouped;
    int                 count{0}; 

    SplitGene( SplitCoVfasta &parent, std::string gene);
    
    /// record identified and ...?
    target_count& check_msa_rec     (auto& record);
    target_count& check_rec         (auto& record);

    SeqPos find_ampl_pos            (const oligo_seq_t& target);
    void target_pattern             (target_q  &tq, const oligo_seq_t &full_target, long ampl_beg);
    void evaluate_target            (target_q  &tq, const oligo_seq_t &full_target, long ampl_beg);
    void evaluate_target_primer     (pattern_q &pq, const oligo_seq_t &full_target, long offset);
    void evaluate_msa_target        (target_q  &tq, const msa_seq_t &full_target);
    void target_msa_pattern         (target_q  &tq, const msa_seq_t &full_target);
    void evaluate_msa_target_primer (pattern_q &pq, const msa_seq_t &full_target);
    void align_to_msa   (pattern_q &oligo_pattern_quality, const msa_seq_t &full_target   );  // not used
    void re_align       (pattern_q &oligo_pattern_quality, oligo_seq_t &oligo_target);  // to exact target

    bool reconstruct_msa_seq(const msa_seq_t &s, oligo_seq_t &seq,
                         long msa_beg, long msa_end, int tent_len);
    void write_msa_grouped();
    void write_grouped();
    bool quick_check(const oligo_seq_t &target, const oligo &primer, long offset); // check if the primer matches the target at the expected position
};    

/// Given a GISAID fasta file, split it into detally evaluated PCR targets defined by a set of oligos from another small fasta file.
///
/// There are pver 17 000 000 SARS-CoV-2 sequences in GISAID database, and the number is growing.
/// Our goal is to make this information more manejable reducing it to a few hundred thousand of fasta sequences 
/// where each item reprecent a group of identical sequences of the amplicon for that PCR tardet. 
/// Additional splitting of that grouping will be saved into different files.each grouping sequences additionaly by
/// year, month, day, continent, country and clade. 
/// The id or description of each resulting fasta item specify in a defined format the number of original isolates 
/// grouped in that item followed representative metadata includin an EPI_ISL isoloate name, country, region, clade, pango lineage, etc.
/// It also include a detailed evaluation of the PCR target sequence, including the number of mismatches, Ns, gaps, and critical mismatches
/// for each primer 
///  
/// Objects 
class SplitCoVfasta
{
 public:
    std::filesystem::path fasta;
    GISAID_format         format    {GISAID_format::fasta};
    std::filesystem::path dir       {fasta.parent_path()};
    std::string           fasta_name{fasta.stem().string()},
                          from, to;
    int                   flank{5},
                          crit_term_nt{4},
                          crit_mm{3},
                          crit_N{3},
                          crit_gap{2};  // todo how to check partial seqs?
    double                match;
    bool                  check_date = !(from.empty() && to.empty());
    msa_seq_t             msa_ref;
    oligo_seq_t           ref_seq;      
    std::vector<int>      msa_pos;


 private:
    std::vector<SplitGene> genes;
    Metadata               metadata;
    std::ifstream          open_metadata() const;
    std::unordered_map<std::string, size_t>  parse_metadata_header(std::ifstream &metadata_file);
    
    /// extract isolate name from virus_name like "BTC-4694" from "hCoV-19/United Arab Emirates/BTC-4694/2021"
    static std::string_view isolate(const std::string_view virus_name);

 public:
    SplitCoVfasta(const std::filesystem::path& fasta,
                                           int flank, 
                                        double match, 
                                    std::string from, 
                                    std::string to)
    : fasta{fasta}, flank{flank}, match{match}, from{from}, to{to}
    {}

    void add_gene(const std::filesystem::path& oligos, std::string gene, 
                     std::string fw="", std::string rv="")  ///< todo: implement conditional split
    {
        genes.emplace_back(*this, gene).read_oligos(oligos);
    }

    GISAID_format check_format();
    void split();
    bool extract_msa_seq(const msa_seq_t &full_msa_seq, 
                               msa_seq_t &msa_fragment, 
                             oligo_seq_t &reconstructed_seq,
                        long msa_beg, long msa_end, int tent_len = 0) ;
 private:
    void update_target_count(target_count& tc, parsed_id& pid);
    void split_msa( );
    void split_fasta( );
    void set_ref_pos();
    void set_msa_ref_pos();

    void parse_id_allnuc(const std::string& id, parsed_id& pid);  ///< parse the id of the original fasta record, only if not found in metadata
    void parse_id       (const std::string& id, parsed_id& pid);  ///< parse the id of the original fasta record, only if not found in metadata

    void read_metadata();
    
};

} // namespace cov