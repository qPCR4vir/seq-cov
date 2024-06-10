#pragma once

#include <string>
#include <vector>
#include <filesystem>
#include <unordered_map>

#include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/core/debug_stream.hpp>         // for debug_stream

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

namespace cov
{

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

struct oligo
{
    oligo_seq_t   seq;
    std::string   name, code;
    int           beg{0}, end{0}, match{0};
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
};
struct country_count
{
    std::string country;
    int         count = 0;
};
struct day_count
{
    std::unordered_map<std::string, country_count> countries;
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
    std::unordered_map<std::string, cov::year_count> years;
    target_q target;
    int count = 0;
};
struct parsed_id
{
    std::string country, isolate;
    int         year, month, day;
    long        EPI_ISL;
};

class SplitCoVfasta;
struct SplitGene
{
    SplitCoVfasta const &parent;
    const std::string   gene;
    int                 forw_idx{}, rev_idx{};    
    std::vector<oligo>  f_primers, r_primers, probes_s, probes_a;

    bool read_oligos(const std::filesystem::path& oligos);

    bool                split, 
                        group,  // we are asked to split, group or ignore this?
                        ignore{!(split || group)};
    

    static constexpr int notfound = std::numeric_limits<int>::lowest(); 
    

    int                 fw_match, rv_match;     // deprecate
    msa_seq_t           target;      
    std::vector<int>    msa_target_pos;
    int                 beg{0}, end{0}, len{0}, count{0}; 
    const std::string   start{gene+"|"};

    using grouped_by_seq = std::unordered_map<msa_seq_t, target_count>;
    grouped_by_seq      grouped; 

    SplitGene( SplitCoVfasta const &parent, std::string gene);
    
    /// record identified and ...?
    target_count& check_rec(auto& record);
    void evaluate_target(target_q  & tq, msa_seq_t& target);

    bool reconstruct_seq(const msa_seq_t& s, oligo_seq_t& seq, 
                            int& beg, int& end, std::vector<int>& msa_pos, int tent_len);

    void write_grouped ();
};    


class SplitCoVfasta
{
 public:    


 public:
    std::filesystem::path fasta;
    std::filesystem::path dir       {fasta.parent_path()};
    std::string           fasta_name{fasta.filename().string()},
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

 public:
    SplitCoVfasta(const std::filesystem::path& fasta,
                                           int flank, 
                                        double match, 
                                    std::string from, 
                                    std::string to)
    : fasta{fasta}, flank{flank}, match{match}, from{from}, to{to}
    {}

    void add_gene(const std::filesystem::path& oligos, std::string gene, 
                     std::string fw="", std::string rv="")  // todo implement conditional split
    {
        genes.emplace_back(*this, gene).read_oligos(oligos);
    }

    void split_fasta( );
    void set_ref_pos();
    void parse_id(const std::string& id, parsed_id& pid);
};

} // namespace cov