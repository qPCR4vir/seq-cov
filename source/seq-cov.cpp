#include <string>
#include <charconv>
#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <unordered_map>
#include <ranges>
#include <stdexcept>
#include <execution>

	
// #include <seqan3/std/ranges>                    // include all of the standard library's views
#include <seqan3/alphabet/views/all.hpp>        // include all of SeqAn's views 
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>
//#include <seqan3/alphabet/range/hash.hpp>
#include <seqan3/core/debug_stream.hpp>         
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>

#include "seq-cov.hpp"

namespace cov
{
/*
struct dna_deg : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna15; // instead of dna5
    template <typename alph>
    using sequence_container = seqan3::dna15_vector;
};
*/


// struct SplitGene
 
SplitGene::SplitGene(SplitCoVfasta &parent, std::string gene)
        : parent{parent}, gene{gene}
{
    
    //if (split) file_E.open( );  // todo implement conditional split
    std::cout << std::boolalpha  
        << "\nGoing to split: " << (parent.dir / parent.fasta_name).string() 
        << ", gene " << gene << " into "
        << (parent.dir / (gene + "." + parent.fasta_name)).string();
    
}

bool SplitGene::read_oligos(const std::filesystem::path& path_oligos)
{
        if (path_oligos.empty()) return false;   // todo: more checks?

        seqan3::debug_stream << "\nGoing to load: " << path_oligos.string();
        seqan3::sequence_file_input<OLIGO> file_in{path_oligos};
        std::string fw, rv;
        ref_beg = ref_end = 0;
        f_primers.clear();
        r_primers.clear();
        probes_s.clear();
        probes_a.clear();
        int forw_idx{0}, rev_idx{0};
        all_oligos.clear();

        for (auto & primer : file_in)
        {
            std::string id = std::move(primer.id());
            seqan3::debug_stream << "\nGoing to check: " << id << "\n" << primer.sequence();

            // parse beg, end from id = >SARS_NF+A -13900 MN908947.3: Seq pos: 28775-28794
            oligo pr;
            std::stringstream ss{id};
            ss >> pr.name >> pr.code ;

            std::string seq_pos = id.substr(id.find("Seq pos: ") + 9);
            size_t dash_pos = seq_pos.find("-");
            pr.ref_beg = std::stoi(seq_pos.substr(0, dash_pos));
            pr.ref_end = std::stoi(seq_pos.substr(dash_pos + 1));
            // todo check if beg, end are valid
            seqan3::debug_stream << " with beg: " << pr.ref_beg << " and end: " << pr.ref_end;
            pr.seq  = std::move(primer.sequence());
            pr.len  = pr.seq.size();
            pr.match = std::round(parent.match * double(pr.len) / 100.0);
            seqan3::debug_stream << " with minimum matches:" << pr.match;

            if (pr.ref_beg < pr.ref_end)   
            {
                pr.reverse = false;

                // throw std::runtime_error if incorrect lenth of the primer
                if (pr.len != (pr.ref_end - pr.ref_beg + 1)) 
                    throw std::runtime_error{"Forward primer " + pr.name + " has incorrect length: " + std::to_string(pr.len) +
                    " instead of " + std::to_string(pr.ref_end - pr.ref_beg + 1) + " at position " + std::to_string(pr.ref_beg) 
                    + " to " + std::to_string(pr.ref_end)};

                if ( !ref_beg || ref_beg > pr.ref_beg)  // extern forward primer
                {
                    ref_beg = pr.ref_beg;
                    forw_idx = f_primers.size();
                }   
                f_primers.push_back(pr);   
                all_oligos.push_back(pr);                  
                continue;
            }
            else   
            {
                pr.reverse = true;

                // throw std::runtime_error if incorrect lenth of the primer
                if (pr.len != (pr.ref_beg - pr.ref_end + 1)) 
                    throw std::runtime_error{"Reverse primer " + pr.name + " has incorrect length: " + std::to_string(pr.len) +
                    " instead of " + std::to_string(pr.ref_beg - pr.ref_end + 1) + " at position " + std::to_string(pr.ref_beg) 
                    + " to " + std::to_string(pr.ref_end)};

                if ( !ref_end || ref_end < pr.ref_beg)  // extern reverse primer
                {
                    ref_end = pr.ref_beg;
                    rev_idx = r_primers.size();
                }
                r_primers.push_back(pr);
                all_oligos.push_back(pr);    
                continue;
            }
        }
        ref_len = ref_end - ref_beg + 1;

        /*for (auto & primer : f_primers) all_oligos.push_back(&primer);
        for (auto & primer : r_primers) all_oligos.push_back(&primer);
        for (auto & primer : probes_s) all_oligos.push_back(&primer);
        for (auto & primer : probes_a) all_oligos.push_back(&primer);*/

        return true;
}

bool SplitGene::reconstruct_seq(const msa_seq_t& full_msa_seq, oligo_seq_t& seq, long msa_beg, long msa_end, int tent_len)
{
    //seqan3::debug_stream << "\nReconstructing sequence from " << msa_beg << " to " << msa_end << '\n';
    
    seq.clear();
    seq.reserve(tent_len);
    bool reverse = (msa_beg > msa_end);
    // go through the sequence and eliminate gaps to reconstruct the original sequence
    if (reverse) 
    {
        // original sequence is in reverse order

        for ( int i = msa_beg; i >= msa_end; --i)
        {
            if (full_msa_seq[i].holds_alternative<seqan3::gap>() ) continue;
            seq.push_back(full_msa_seq[i].convert_unsafely_to<oligo_seq_alph>().complement());
        }
        return true;
    }
    // original sequence is in the same order
    for ( int i = msa_beg; i <= msa_end; ++i)
    {
        if (full_msa_seq[i].holds_alternative<seqan3::gap>() ) continue;
        seq.push_back(full_msa_seq[i].convert_unsafely_to<oligo_seq_alph>());
    }
    return true;
}
void SplitGene::re_align(pattern_q &pq, const oligo_seq_t &oligo_target)
{
    
    seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                        << pq.primer.seq <<" - oligo seq\n" 
                        << oligo_target << " - target seq\n" ;

    auto output_config = seqan3::align_cfg::output_score{}          | 
                         seqan3::align_cfg::output_begin_position{} |
                         seqan3::align_cfg::output_end_position{}   | 
                         seqan3::align_cfg::output_alignment{};
    auto method = seqan3::align_cfg::method_global{};
    seqan3::align_cfg::scoring_scheme     scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, 
                                                    seqan3::mismatch_score{-1}}};
    /*seqan3::align_cfg::gap_cost_affine gap_costs{seqan3::align_cfg::open_score{0}, 
                                                    seqan3::align_cfg::extension_score{-1}};*/
    auto config = method |  scheme |output_config; // gap_costs |
    //seqan3::debug_stream << "\nTarget: " << oligo_target << "\nPrimer: " << pq.primer.seq << '\n';
    
    for (auto const & res : seqan3::align_pairwise(std::tie(oligo_target, pq.primer.seq), config))
    //auto results = seqan3::align_pairwise(std::tie(target, pq.primer.seq), config);
    //seqan3::debug_stream << " Going to check results\n"; 
    //for (auto & res : results)
    // if (res.score() > pq.primer.match)  // primer found. 
    {
        //seqan3::debug_stream << "Alignment: " << res.alignment() << " Score: "     << res.score() ;
        //seqan3::debug_stream << ", Target: (" << res.sequence1_begin_position() << "," << res.sequence1_end_position() << ")";
        //seqan3::debug_stream << ", Primer: (" << res.sequence2_begin_position() << "," << res.sequence2_end_position() << "). " ;
    
            //brief Helper function to print elements of a tuple separately.
            auto&[a_t, a_p] = res.alignment();
            pq.pattern.clear();
            
            int len_a_pr = a_p.size();
            int len_a_tg = a_t.size(); 
            int len = std::max(len_a_pr, len_a_tg);

            //seqan3::debug_stream << " cont...\n" ;
            int o_p = res.sequence2_begin_position();
            if (o_p) pq.pattern = std::string(o_p, '-');

            for (int i = 0; i < len; ++i)
            {
                if (i >= len_a_pr)  // extra nt in aligned target ?? ignore ?
                {
                    if (o_p >= pq.primer.len) continue  ;  // ??
                    o_p++;
                    pq.pattern.push_back(a_t[i].to_char());   
                    if  (a_t[i] == 'N'_dna15)  
                    {
                        pq.N++;
                        continue;
                    }
                    pq.mm++;
                    if (o_p >= pq.primer.len - parent.crit_term_nt)      pq.crit++;  // ?????
                    continue;
                }
                if (i >= len_a_tg)  // no more aligned target?
                {
                    if (o_p >= pq.primer.len) continue  ;  // ??
                    o_p++;
                    pq.pattern.push_back('-');
                    pq.mm++;
                    if (o_p >= pq.primer.len - parent.crit_term_nt)                   pq.crit++;
                    continue;
                }
                if (a_p[i].holds_alternative<seqan3::gap>() &&
                    a_t[i].holds_alternative<seqan3::gap>()    )       continue;  // ??

                if (a_t[i].holds_alternative<seqan3::gap>()    )     // deletion in target (or insertion in ref?)
                {
                    if (o_p >= pq.primer.len) continue  ;  // ??
                    o_p++;
                    pq.pattern.push_back('-');  
                    continue;
                }
                auto t = a_t[i].convert_unsafely_to<oligo_seq_alph>();        
                if (a_p[i].holds_alternative<seqan3::gap>())  // insertion in target 
                {
                    if (!o_p) continue;  // ignore insertions before primer
                    // assume ???? that the inserted base is not in the primer
                    pq.pattern.push_back(t.to_char());  // inserted in the pattern but not advance the primer
                    pq.mm++;
                    if (pq.primer.len - o_p <= parent.crit_term_nt)                   pq.crit++;  // ??
                    continue;
                }
                auto p = a_p[i].convert_unsafely_to<oligo_seq_alph>();   
                o_p++;
                if (!mismatch.score(t, p))     
                { 
                    if  (t == 'N'_dna15)  
                    {
                        pq.N++;
                        pq.pattern.push_back('N');  
                        continue;
                    }                    
                    pq.pattern.push_back('.');  
                    continue;
                }   
                pq.pattern.push_back(t.to_char());   
                pq.mm++;
                if (o_p >= pq.primer.len - parent.crit_term_nt)                   pq.crit++;
            }
            if (res.score() > pq.primer.match)  
                pq.Q = pq.mm + pq.crit * 4;  // primer found. 
            else pq.Q = 1000;

            seqan3::debug_stream << a_p << " - aligned oligo\n"
                                 << pq.pattern << " - pattern\n"
                                 << a_t << " - aligned target\n" ;
            return;
    }
}
void SplitGene::align(pattern_q &pq,  ///< oligo_pattern_quality
                    const msa_seq_t &full_target)
{
    oligo_seq_t oligo_target;
    msa_seq_t   msa_fragment;
    cov::oligo &primer = pq.primer;
    parent.extract_seq(full_target, msa_fragment, oligo_target, primer.msa_beg, primer.msa_end, primer.seq.size());

    pq.pattern.clear();
    pq.pattern.reserve(primer.seq.size()+1);
    
    // go through the sequence and eliminate gaps to reconstruct the original sequence

    for ( int i = 0, p=0; i < primer.msa_len; ++i)
    {
        if (primer.msa_ref[i].holds_alternative<seqan3::gap>() &&
              msa_fragment[i].holds_alternative<seqan3::gap>()    )           continue;

        if (msa_fragment[i].holds_alternative<seqan3::gap>()    )     // deletion in target (or insertion in ref?)
        {
            pq.pattern.push_back('-');  
            continue;
        }
        auto s = msa_fragment[i].convert_unsafely_to<oligo_seq_alph>();        
        if (primer.msa_ref[i].holds_alternative<seqan3::gap>())  // insertion in target (or deletion in ref?)
        {
            if (!mismatch.score(s, primer.seq[p]))               // the inserted base is in the primer ????
            { 
                pq.pattern.push_back('.'); 
                ++p; 
                continue;
            }
            // assume ???? that the inserted base is not in the primer
            pq.pattern.push_back(msa_fragment[i].to_char());  // inserted in the pattern but not advance the primer
            pq.mm++;
            if (primer.len - p <= parent.crit_term_nt)                   pq.crit++;
            continue;
        }

        if (!mismatch.score(s, primer.seq[p]))     
        { 
            pq.pattern.push_back('.');  
            ++p;
            continue;
        }
        ++p;        
        pq.pattern.push_back(s.to_char());   
        if  (s == 'N'_dna15)  
        {
            pq.N++;
            continue;
        }
        pq.mm++;
        if (i >= primer.len - parent.crit_term_nt)                   pq.crit++;
    }
    pq.Q = pq.mm + pq.crit * 4;
    
     seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<" - oligo \n" 
                             << pq.pattern << '\n'
                             << oligo_target << '\n'
                             << pq.primer.msa_ref << " -ref" << '\n'
                             << msa_fragment << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 
}
void SplitGene::target_pattern(target_q & tq, const msa_seq_t& full_target)
{
    tq.target_pattern.clear();
    tq.target_pattern.reserve(ref_len);

    for ( int i = msa_beg; i <= msa_end; ++i)
    {
        if (parent.msa_ref[i].holds_alternative<seqan3::gap>() &&
               full_target[i].holds_alternative<seqan3::gap>()    )           continue;

        if (parent.msa_ref[i].holds_alternative<seqan3::gap>() ||
               full_target[i].holds_alternative<seqan3::gap>()    )          
        {    
            tq.target_pattern.push_back(full_target[i].to_char());     
            continue;  
        }      
        auto s =    full_target[i].convert_unsafely_to<oligo_seq_alph>();
        auto r = parent.msa_ref[i].convert_unsafely_to<oligo_seq_alph>();
        if (!mismatch.score(s, r))   
        { 
            tq.target_pattern.push_back('.');  
            continue;
        }
        tq.target_pattern.push_back(s.to_char());   
    }
}

void SplitGene::evaluate_target(target_q & tq, const msa_seq_t& full_target)
{
    for (oligo& primer : all_oligos)   
    {
        tq.patterns.emplace_back(primer);  // registering/creating the pattern_q is cheap and fast but difficult to parallelize
    }
    std::for_each(std::execution::par_unseq, tq.patterns.begin(), tq.patterns.end(), [&](pattern_q &pq)
    //for (pattern_q& pq : tq.patterns)   
    {
        evaluate_target_primer(pq, full_target);
    }//
    );
    target_pattern(tq, full_target);

    /*
    for (auto & primer : f_primers) evaluate_target_primer(tq, primer, sq);
    for (auto & primer : r_primers) evaluate_target_primer(tq, primer, sq);
    for (auto & probe  : probes_s ) evaluate_target_primer(tq, probe , sq);
    for (auto & probe  : probes_a ) evaluate_target_primer(tq, probe , sq);
    /* for (auto & pq : tq.patterns)
        seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             //<< target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 */}

void SplitGene::evaluate_target_primer(pattern_q &pq, const msa_seq_t& full_target)
{
    cov::oligo &primer = pq.primer;
    oligo_seq_t oligo_target;
    reconstruct_seq(full_target, oligo_target, primer.msa_beg, primer.msa_end, primer.seq.size());
    if (oligo_target.size() != primer.len)
        //return align(pq, full_target); // todo: try to align the primer with the oligo_target sequence
        return re_align(pq, oligo_target); // todo: try to align the primer with the oligo_target sequence
 
    pq.pattern = std::string(primer.len, '.');
    for (int i = 0; i < primer.len; ++i)
    {
        if (!mismatch.score(primer.seq[i], oligo_target[i]))   continue;

        if  (oligo_target[i] == 'N'_dna15)  
        {
            pq.N++;
            pq.pattern[i] = 'N';
            continue;
        }
        pq.pattern[i] = oligo_target[i].to_char();
        pq.mm++;
        if (i >= primer.len - parent.crit_term_nt)                   pq.crit++;
     }
    pq.Q = pq.mm + pq.crit * 4;
/*     seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             << oligo_target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 */
}

target_count& SplitGene::check_rec(auto& record)
{
    msa_seq_t& full_target = record.sequence();
    count++;

    target_count & target_c = grouped[{full_target.begin()+msa_beg, full_target.begin()+msa_end}]; 
    if (!target_c.count)  // new target sequence
        evaluate_target(target_c.target, full_target);
    
    target_c.count++;
    return target_c;   
}

void SplitGene::write_grouped ()
{
    using types    = seqan3::type_list<msa_seq_t, std::string>;
    using fields   = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using record_t = seqan3::sequence_record<types, fields>;
    using sgr_t    = decltype(grouped)::value_type;  // std::unordered_map<msa_seq_t, target_count>;
    using sgr_p    = sgr_t*;

    // std::filesystem::path gr  = parent.dir / (gene + ".grouped-" + parent.fasta_name);
       std::filesystem::path grdc= parent.dir / (gene + ".daily-coountry-"   + parent.fasta_name);
       std::filesystem::path grd = parent.dir / (gene + ".daily-"            + parent.fasta_name);
    // std::filesystem::path grm = parent.dir / (gene + ".monthly-" + parent.fasta_name);

    seqan3::debug_stream << "Going to write: " << grd.string() << "\n" ;
    seqan3::sequence_file_output<fields, seqan3::type_list<seqan3::format_fasta> >  file_e_grdc{grdc};  // , file_e_gr {gr },  file_e_grm{grm};
    seqan3::sequence_file_output<fields, seqan3::type_list<seqan3::format_fasta> >  file_e_grd {grd };  // , file_e_gr {gr },  file_e_grm{grm};
    
    std::vector<sgr_p> gr_v;
    gr_v.reserve(grouped.size());
    for (sgr_t& sgr : grouped) gr_v.push_back(&sgr);  // vector of pointers to grouped target_count sequences

    std::sort(gr_v.begin(), gr_v.end(), []( sgr_p &a,  sgr_p &b)    {return a->second.count > b->second.count;});

    for (sgr_p           sg :  gr_v            )   // pointers to grouped target_count   seq: target_count
    {
        // discard sequences with too many N.
        if (std::all_of(sg->second.target.patterns.begin(), sg->second.target.patterns.end(), 
                        [&](pattern_q &pq) {return pq.N > parent.crit_N;})) 
            continue;

        for (auto& [year,    yc]:  sg->second.years)
        for (auto& [month,   mc]:  yc.months   )  
        for (auto& [day,     dc]:  mc.days     )
        {    
            for (auto& [country, cc]:  dc.countries)
            {
                auto id = std::format("{} |{:04d}-{:02d}-{:02d}|{}|{}|{}", 
                            cc.id.EPI_ISL, year, month, day, cc.count, country, cc.id.isolate);

                for (auto& pq : sg->second.target.patterns)
                {
                    id += std::format("|{}_Q_{}_mm_{}_N_{}_crit_{}_pat_{}", 
                                    pq.primer.name, pq.Q, pq.mm, pq.N, pq.crit, pq.pattern);
                }
                id += "|" + sg->second.target.target_pattern;    
                file_e_grdc.push_back<record_t>( record_t{sg->first, std::move(id)} );
            }

            auto& [country, cc] =  *dc.countries.begin();
            auto id = std::format("{} |{:04d}-{:02d}-{:02d}|{}|{}|{}",
                        cc.id.EPI_ISL, year, month, day, dc.count, country, cc.id.isolate);


            for (auto& pq : sg->second.target.patterns)
            {
                id += std::format("|{}_Q_{}_mm_{}_N_{}_crit_{}_pat_{}", 
                                pq.primer.name, pq.Q, pq.mm, pq.N, pq.crit, pq.pattern);
            }
            id += "|" + sg->second.target.target_pattern;
                
            file_e_grd.push_back<record_t>( record_t{sg->first, std::move(id)} );
        }
    }
}

// class SplitCoVfasta

bool SplitCoVfasta::extract_seq(const msa_seq_t &full_msa_seq, 
                                          msa_seq_t &msa_fragment, 
                                        oligo_seq_t &reconstructed_seq,
                        long msa_beg, long msa_end, int tent_len )  // = 0
{
    //seqan3::debug_stream << "\nReconstructing sequence from MSA pos " << msa_beg << " to " << msa_end << '\n';
    //seqan3::debug_stream <<  msa_pos.size() << '\t';

    reconstructed_seq.clear();
    if (tent_len) reconstructed_seq.reserve(tent_len);

    bool reverse = (msa_beg > msa_end);
    long len = msa_end - msa_beg + 1;
    // go through the sequence and eliminate gaps to reconstruct the original sequence
    if (reverse) 
    {
        // original sequence is in reverse order
        msa_fragment.reserve(-len);
        for ( int i = msa_beg; i >= msa_end; --i)
        {
            if (full_msa_seq[i].holds_alternative<seqan3::gap>() )
            {
                msa_fragment.push_back(seqan3::gap{});
                continue;    
            }
            auto b = full_msa_seq[i].convert_unsafely_to<oligo_seq_alph>().complement();
            msa_fragment.push_back(b);
            reconstructed_seq.push_back(b);
        }
        return true;
    }
    // original sequence is in the same order
    msa_fragment.reserve(len);
    for ( int i = msa_beg; i <= msa_end; ++i)
    {
            auto b = full_msa_seq[i];
            msa_fragment.push_back(b);
            if (b.holds_alternative<seqan3::gap>() ) continue;
            reconstructed_seq.push_back(b.convert_unsafely_to<oligo_seq_alph>());
    }
    return true;
}
void SplitCoVfasta::set_ref_pos( )
{
    // Initialise a file input object with a FASTA file.
    seqan3::sequence_file_input<MSA> file_in{fasta};

    auto&& ref_rec = *file_in.begin();
    msa_ref = std::move(ref_rec.sequence());
    int msa_len = msa_ref.size();

    ref_seq.reserve(32000);
    msa_pos.clear();
    msa_pos.reserve(ref_seq.capacity());

    // go through the sequence and eliminate gaps to reconstruct the original sequence
    for ( int i = 0; i < msa_len; ++i) 
    {
        if (msa_ref[i].holds_alternative<seqan3::gap>())  continue;
        ref_seq.push_back(msa_ref[i].convert_unsafely_to<oligo_seq_alph>());
        msa_pos.push_back(i);
    }

    seqan3::debug_stream << "\n\nMSA Reference: " << ref_rec.id() << ",\t MSA lenth = " << msa_len << ", reference lenth = " << ref_seq.size() << '\n';
    
    // for (auto & gene : genes) gene.set_ref_pos(); Set/Check correct positions of primers on the reference sequence
    for (auto & gene : genes) 
    {
        // set gene target positions, First check ref_beg and ref_end were alrready set (in read_oligos)
        if (gene.ref_beg == gene.ref_end) // run time error if not set
            throw std::runtime_error{"Primer " + gene.gene + " has no reference positions set"};
        gene.msa_beg = msa_pos[gene.ref_beg-1];
        gene.msa_end = msa_pos[gene.ref_end-1];
        gene.msa_len = gene.msa_end - gene.msa_beg + 1;
        
        // check all primers
        for (auto & primer : gene.all_oligos)
        {
            if (primer.ref_beg == primer.ref_end) // run time error if not set
                throw std::runtime_error{"Primer " + primer.name + " has no reference positions set"};
            primer.msa_beg = msa_pos[primer.ref_beg-1];
            primer.msa_end = msa_pos[primer.ref_end-1];
            primer.msa_len = primer.msa_end - primer.msa_beg + 1;
            extract_seq(msa_ref, primer.msa_ref, primer.ref_seq, primer.msa_beg, primer.msa_end, primer.seq.size());
        }

    }

}

void SplitCoVfasta::parse_id(const std::string& id, parsed_id& pid)
{
    // >hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30|China
    // d_2019-12-30_x_1_c_China - Wuhan_EPI_02124_i_WIV04_p_E_Sarbeco_F2_Q_0_mm_0_N_0_crit_0_pat_......................
    // >hCoV-19/Sweden/6332562676/2022|EPI_ISL_16076002|2022-11-28|Europe
    // hCoV-19/Wales/PHWC-PRWWBD/2021|EPI_ISL_6145484|2021-10-19|Europe
    // hCoV-19/Indonesia/JK-GSILab-624880/2021|EPI_ISL_3230236|2021-07-25|Asia
    // hCoV-19/Indonesia/JK-GSILab-622248/2021|EPI_ISL_3230241|2021-07-16|Asia
    // hCoV-19/USA/AK-PHL10443/2021|EPI_ISL_3232474|2021-07-13|NorthAmerica
    // hCoV-19/Ireland/D-NVRL-Z21IRL04853/2021|EPI_ISL_8349991|2021-12-13|Europe

    std::size_t country_beg = 8;
    std::size_t country_end = id.find('/',             8) - 1;          
    std::size_t country_len = country_end - country_beg   + 1;
    std::size_t isolate_beg = country_end                 + 2;
    std::size_t isolate_end = id.find('/', isolate_beg)   - 1;
    std::size_t isolate_len = isolate_end - isolate_beg   + 1;
    std::size_t EPI_ISL_beg = isolate_end+1 + 6 + 8 ;               
    std::size_t EPI_ISL_end = id.find('|', EPI_ISL_beg+1) - 1;
    std::size_t EPI_ISL_len = EPI_ISL_end - EPI_ISL_beg   + 1;
    std::size_t year_beg    = EPI_ISL_end                 + 2;   
    std::size_t month_beg   = year_beg  + 5;
    std::size_t day_beg     = month_beg + 3;  
    std::size_t region_beg  = day_beg   + 3;

    pid.isolate = id.substr(isolate_beg, isolate_len);                 
    pid.EPI_ISL = id.substr(EPI_ISL_beg, EPI_ISL_len); 
    pid.country = id.substr(region_beg) + " - " + id.substr(country_beg, country_len) ;  

    auto s = id.data() ;
    auto y = s + year_beg ;
    auto m = s + month_beg;
    auto d = s + day_beg  ;    
    
    std::from_chars(y, m    , pid.year );   
    std::from_chars(m, d    , pid.month); 
    std::from_chars(d, d + 3, pid.day  );

/*     seqan3::debug_stream << "\n" << pid.isolate << " - EPI: " << pid.EPI_ISL << " - " 
                         << pid.year << " - " << pid.month << " - " << pid.day 
                         << " - " << pid.country << '\n' ; */

}

void SplitCoVfasta::split_fasta( )
{
    set_ref_pos();

    // Initialise a file input object with a FASTA file.
    seqan3::sequence_file_input<MSA> file_in{fasta};

    long t{0L}, m{(1L<<18)-1};
    seqan3::debug_stream << "\nchunk - m= " << m << "\n" ; 


    for (auto && record : file_in)             // read each sequence in the file
    {
        // seqan3::debug_stream << "\n" << record.id() << '\n' ;
        parsed_id pid;
        parse_id(record.id(), pid);


        std::for_each(std::execution::par_unseq, genes.begin(), genes.end(), [&](auto& gene) 
        //for (auto & gene : genes)  // todo: parallelize (std::execution::par)
        {
            target_count& tc = gene.check_rec(record);
            year_count& yc = tc.years[pid.year];
            yc.count++;
            month_count& mc = yc.months[pid.month];
            mc.count++;
            day_count& dc = mc.days[pid.day];
            dc.count++;
            country_count& cc = dc.countries[pid.country];
            if (!cc.count) cc.id = pid;
            cc.count++;
        }//
        );

        if (!(++t & m))                      // print a dot every 2^18 sequences for progress indication
        {
            // seqan3::debug_stream << '.' ;
        
            seqan3::debug_stream << "\tT= " << t  << "\n" ;
            for (auto & gene : genes)
                seqan3::debug_stream << gene.gene <<"= " << gene.count 
                                     << ". Grouped: "    << gene.grouped.size() << "\n" ; 
        }
        //if (t>700000) break;
    }
    seqan3::debug_stream << "\nTotal= " << t  << "\n" ; 

    for (auto & gene : genes)                // write grouped sequences for each gene/target
    {
        seqan3::debug_stream << gene.gene <<"= " << gene.count 
                             << ". Grouped: "    << gene.grouped.size() << ". " ; 
        gene.write_grouped();
    }

}

/* https://docs.seqan.de/seqan3/main_user/cookbook.html 
Reading records in chunks
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/chunk.hpp>
 
int main()
{
    seqan3::sequence_file_input fin{std::filesystem::current_path() / "my.fastq"};
 
    // `&&` is important because seqan3::views::chunk returns temporaries!
    for (auto && records : fin | seqan3::views::chunk(10))
    {
        // `records` contains 10 elements (or less at the end)
        seqan3::debug_stream << "Taking the next 10 sequences:\n";
        seqan3::debug_stream << "ID:  " << (*records.begin()).id() << '\n'; // prints first ID in batch
    }
The example above will iterate over the file by reading 10 records at a time. If no 10 records are available anymore, it will just print the remaining records.
*/

}