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
        
        SplitGene(SplitCoVfasta const &parent, std::string gene, bool split, bool group);
        
        /// record identified and ...?
        bool check(auto& record);

        bool check_rec(auto& record);

        SeqGr set_seq_pos(const sequence_type& s);

        void write_grouped ();
        
        SplitGene& set_forw(const  std::string&   pr  );
        SplitGene& set_forw(seqan3::dna5_vector   forw);
        SplitGene& set_rev (const  std::string&   pr  );
        SplitGene& set_rev (seqan3::dna5_vector   rev );
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

    void split_fasta( );
};
 