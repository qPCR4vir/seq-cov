/**
* Evaluate the quality of a PCR across time, regions and CoV-2 clades.
*
* @file seq_cov.hpp
*
*/
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

/** Hashing a range of characters.
 * 
 * \tparam urng_t The type of the range; Must model std::ranges::input_range and the reference type of the range of the
 *                range must model seqan3::semialphabet.
 *
 */
template <std::ranges::input_range urng_t>
    requires seqan3::semialphabet<std::ranges::range_reference_t<urng_t>>
struct hash
{
    /** Compute the hash for a range of characters.
     * 
     * \tparam urng2_t  The same as `urng_t` (+- cvref); used to get forwarding reference in the interface.
     * \param[in] range The input range to process. Must model std::ranges::input_range and the reference type of the
     *                  range of the range must model seqan3::semialphabet.
     * \returns size_t.
     */
    template <std::ranges::input_range urng2_t>
        requires seqan3::semialphabet<std::ranges::range_reference_t<urng2_t>>
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


/**
 * Contains information and methods for processing PCR targets.
 * 
 * Evaluated PCR targets defined by a set of oligos from an small fasta file.
 * It include a detailed evaluation of the PCR target sequence, including the patterns and number of mismatches, Ns, gaps, and critical mismatches
 * for each oligo and the match pattern for the amplicon.  
 * Some functions are now split into two versions: one for the old MSA format and a new version for the new 'raw' sequence format. 
 * As the sequences are now not pre-aligned we don't know the exact position 
 * of the amplicon and of each oligo/primer. 
 * Thus, for each new sequence we have the following cases: 
 * 1-exact known position and sequence: just update the counters. 
 * 2-right position but new sequence: 
 *         simple check that one of external primers matches the expected position: 
 *          reuse the newly retrieved parsed_id pid to hold all the patterns, that have to be generated first. 
 * 3-New position: 
 *     first discard wrong seq/pid 
 *     them we need to align one of the external primers to the expected target first 
 *         (if not found expand by 1000 that region, 
 *             if fails use the whole sequence, 
 *             and try inverted too) 
 *     and readjust the coordinates of the amplicon and repeat to check if we are now in 1- or 2-. 
 */
struct PCRSplitter
{
    seqan3::nucleotide_scoring_scheme<int8_t> mismatch; ///< Scoring scheme for mismatches (default is hamming distance)

    SplitCoVfasta &   parent;             ///< Parent SplitCoVfasta instance
    const std::string pcr_name;           ///< PCR name

    PCRSplitter( SplitCoVfasta &parent,   ///< Parent SplitCoVfasta instance
                std::string   pcr_name    ///< PCR name
             );
    
    std::vector<oligo>   f_primers,   ///< Forward primers for this PCR
                         r_primers,   ///< Reverse primers for this PCR
                         probes_s,    ///< Probes sense direction
                         probes_a;    ///< Probes antisense direction

    std::vector<oligo> all_oligos;    ///< All oligos (primers and probes)

    int extern_forw_idx{},      ///< Index of the external primer in f_primers 
        extern_rev_idx{};       ///< Index of the external primer in r_primers

    long ref_beg{0},     ///< Amplicon start position in the reference sequence
         ref_end{0};     ///< Amplicon end position in the reference sequence

    long msa_beg{0},     ///< Amplicon start position in the MSA
         msa_end{0};     ///< Amplicon end position in the MSA

    int msa_len{0},      ///< Length of the amplicon in MSA
        ref_len{0};      ///< Length of the amplicon in the reference

    // Grouping types for targets
    using grouped_by_msa_seq = std::unordered_map<msa_seq_t, target_count, hash<msa_seq_t>>;
    using grouped_by_seq     = std::unordered_map<oligo_seq_t, target_count, hash<oligo_seq_t>>;

    grouped_by_msa_seq msa_grouped;  ///< Grouped targets based on MSA sequences
    grouped_by_seq     grouped;      ///< Grouped targets based on raw nucleotide sequences

    int                count{0};     ///< Total count of identified targets

    /// Reads oligo definitions from a fasta file including positions and sequences
    bool read_oligos(const std::filesystem::path& oligos  ///< Path to the oligos file
                    );

    /// Processes an MSA record and returns a reference to a target_count keep in the msa_grouped map
    target_count & check_msa_rec(auto & record        ///< MSA record to process
                                );

    /// Processes a raw record and returns a reference to a target_count keep in the grouped map
    target_count & check_rec   (auto & record         ///< non-MSA record to process
                               );

    /// Finds the amplicon positions in the target sequence
    SeqPos find_ampl_pos(const oligo_seq_t & target   ///< The nucleotide sequence of the target
                        );

    /// Determines the target pattern from the full sequence
    void target_pattern(target_q & tq,                    ///< Target pattern output structure
                        const oligo_seq_t & full_target,  ///< Full target sequence
                        long ampl_beg                     ///< Amplicon beginning position
                       );

    /// Fully evaluates the target sequence for match quality
    void evaluate_target(target_q & tq,                    ///< Target quality output structure
                         const oligo_seq_t & full_target,  ///< Full target sequence
                         long ampl_beg                     ///< Amplicon beginning position
                        );

    /// Evaluates the alignment of a primer against the target sequence
    void evaluate_target_primer(pattern_q & pq,                   ///< primer input & pattern quality output structure
                                const oligo_seq_t & full_target,  ///< Full target sequence
                                long offset                       ///< Offset applied to the amplicon position respective to the position in the reference sequence
                               );

    /// Evaluates an MSA target sequence
    void evaluate_msa_target(target_q & tq,                 ///< Target quality output structure
                             const msa_seq_t & full_target  ///< Full MSA sequence
                            );

    /// Determines the target pattern from an MSA sequence
    void target_msa_pattern(target_q & tq,                 ///< Target pattern output structure
                            const msa_seq_t & full_target  ///< Full MSA sequence
                           );

    /// Evaluates primer alignment against an MSA target
    void evaluate_msa_target_primer(pattern_q & pq,                ///< Pattern quality output structure
                                    const msa_seq_t & full_target  ///< Full MSA sequence
                                   );

    /// Aligns the target to the MSA sequence (not currently used)
    void align_to_msa(pattern_q & oligo_pattern_quality,  ///< Alignment quality output structure
                      const msa_seq_t & full_target       ///< Full MSA sequence
                     );

    /// Re-aligns the primer to the exact target
    void re_align(pattern_q & oligo_pattern_quality,  ///< Alignment quality output structure
                  oligo_seq_t & oligo_target          ///< Nucleotide sequence of the exact oligo target
                 );

    /// Reconstructs an original sequence from an MSA by eliminating gaps
    bool reconstruct_msa_seq(const msa_seq_t & s,        ///< Full MSA sequence
                             oligo_seq_t     & seq,      ///< Output reconstructed nucleotide sequence
                             long              msa_beg,  ///< Beginning position in MSA
                             long              msa_end,  ///< Ending position in MSA
                             int               tent_len  ///< Tentative length for reservation
                            );

    void write_msa_grouped();     ///< Writes grouped of originally MSA target data to output
    void write_grouped    ();     ///< Writes grouped of originally raw target data to output

    /// Checks whether a primer matches the target at the expected position
    bool quick_check(const oligo_seq_t & target,         ///< Target sequence
                     const oligo       & primer,         ///< Oligo (primer) to check
                     long                offset          ///< Offset applied to the amplicon position respective to the position in the reference sequence
                    );
};    


/**
 * Manages splitting a GISAID FASTA file into evaluated PCR targets using defined oligos.
 * 
 * Given a GISAID fasta file, split it into detally evaluated PCR targets defined by a set of oligos from another small fasta file.
 * There are over 17 000 000 SARS-CoV-2 sequences in the GISAID database, and the number is growing.
 * Our goal is to make this information more manejable reducing it to a few hundred thousand of short fasta sequences 
 * where each item reprecent a group of identical sequences of the amplicon for a PCR target. 
 * Additional splitting of that grouping will be saved into different files, each grouping sequences additionaly by
 * year, month, day, continent, country and clade. 
 * The id or description of each resulting fasta item specify in a defined format the number of original isolates 
 * grouped in that item followed by representative metadata including an EPI_ISL, isoloate name, country, region, clade, pango lineage, etc.
 * It also include a detailed evaluation of the PCR target sequence, including the patterns and number of mismatches, Ns, gaps, and critical mismatches
 * for each oligo and the match pattern for the amplicon.  
 * This information may be used by external code to evaluate the quality of the PCR across time, regions and CoV-2 clades.
 * 
 * Some functions are now split into two versions: one for the old MSA format and a new version for the new 'raw' sequence format. 
 * As the sequences are now not pre-aligned we don't know the exact position 
 * of the amplicon and of each oligo/primer. 
 * We are saving and grouping the possible matching patterns by equal target sequences member map grouped. 
 * Thus, for each new sequence we have the following cases: 
 * 1-exact known position and sequence: just update the counters. 
 * 2-right position but new sequence: 
 *         simple check that one of external primers matches the expected position: 
 *          reuse the newly retrieved parsed_id pid to hold all the patterns, that have to be generated first. 
 * 3-New position: 
 *     first discard wrong seq/pid 
 *     them we need to align one of the external primers to the expected target first 
 *         (if not found expand by 1000 that region, 
 *             if fails use the whole sequence, 
 *             and try inverted too) 
 *     and readjust the coordinates of the amplicon and repeat to check if we are now in 1- or 2-. 
 */
class SplitCoVfasta
{
 public:
    std::filesystem::path fasta;                              ///< Path to the input GISAID FASTA file
    GISAID_format         format {GISAID_format::fasta};      ///< Format of the input file
    std::filesystem::path dir    {fasta.parent_path()};       ///< Directory of the input file
    std::string           fasta_name {fasta.stem().string()}; ///< Name of the file without extension

    std::string           from, to;                                   ///< Date range  
    bool                  check_date = !(from.empty() && to.empty()); ///< Flag to determine whether to check dates

    int                   flank {5};          ///< Flanking region size
    int                   crit_term_nt {4};   ///< Critical terminal nucleotide threshold
    int                   crit_mm {3};        ///< Critical mismatch threshold
    int                   crit_N {3};         ///< Critical ambiguous nucleotide threshold
    int                   crit_gap {2};       ///< Critical gap threshold
    double                match;              ///< Minimum match percentage for primers

    msa_seq_t             msa_ref;   ///< Reference MSA sequence
    oligo_seq_t           ref_seq;   ///< Reference nucleotide sequence
    std::vector<int>      msa_pos;   ///< MSA positions

 private:
    std::vector<PCRSplitter> pcrs;   ///< Container for PCR-specific information

    Metadata               metadata;                  ///< Mapping of isolate names to parsed metadata
    std::ifstream          open_metadata() const;     ///< Opens the metadata file
    std::unordered_map<std::string, size_t> 
                           parse_metadata_header(std::ifstream &metadata_file);  ///< Parses the header of the metadata file

    /// Extract the isolate name from virus_name like "BTC-4694" from "hCoV-19/United Arab Emirates/BTC-4694/2021"
    static std::string_view isolate(const std::string_view virus_name ///< @param virus_name: Full virus name (e.g. "hCoV-19/United Arab Emirates/BTC-4694/2021")
                                   );     
    
 public:
    SplitCoVfasta(const std::filesystem::path& fasta,    ///< Path to the main giant input GISAID FASTA file
                                           int flank,    ///< Flanking region size (currently set to 0? 5?)
                                        double match,    ///< Minimum match percentage for primers
                                    std::string from,    ///< Date range start
                                    std::string to       ///< Date range end
                                )
    : fasta{fasta}, flank{flank}, match{match}, from{from}, to{to}
    {}

    /// Adds a PCR entry and reads its oligo definitions
    void add_gene(const std::filesystem::path& oligos,    ///< Path to the oligos file
                                   std::string pcr_name,  ///< PCR name
                                   std::string fw = "",   ///< Optional forward primer override
                                   std::string rv = ""    ///< Optional reverse primer override
                                )
    {
        pcrs.emplace_back(*this, pcr_name).read_oligos(oligos);
    }

    GISAID_format check_format();  ///< Checks the format of the input FASTA file: MSA or raw sequences

    void split();                  ///< Splits the input file into grouped PCR targets  

    /// Extracts a fragment from an MSA and reconstructs a nucleotide sequence
    bool extract_msa_seq(const msa_seq_t &full_msa_seq,              ///< Full MSA sequence, including the gaps that make the alignment
                               msa_seq_t &msa_fragment,              ///< Output MSA fragment
                             oligo_seq_t &reconstructed_seq,         ///< Output reconstructed nucleotide sequence
                                    long msa_beg,                    ///< MSA beginning position
                                    long msa_end,                    ///< MSA ending position
                                    int tent_len = 0                 ///< Tentative length for sequence reservation
                       ) ;
    
 private:
    /// Updates a target count using parsed metadata
    void update_target_count(target_count& tc,       ///< Target count to update
                             parsed_id& pid          ///< Parsed metadata structure
                            );  

    void split_msa( );         ///< Splits sequences based on MSA data
    void split_fasta( );       ///< Splits the FASTA sequences into groups
    void set_ref_pos();        ///< Sets reference positions for the amplicon
    void set_msa_ref_pos();    ///< Sets MSA reference positions for the amplicon

    void parse_id_allnuc(const std::string& id, parsed_id& pid);  ///< parse the id of the original fasta record, only if not found in metadata
    void parse_id       (const std::string& id, parsed_id& pid);  ///< parse the id of the original fasta record, only if not found in metadata

    void read_metadata();     ///< Reads metadata from the corresponding file
    
};

} // namespace cov