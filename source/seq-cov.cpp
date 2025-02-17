
// file seq-cov.cpp

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
#include <chrono>
	
// #include <seqan3/std/ranges>                    // include all of the standard library's views
#include <seqan3/alphabet/all.hpp>
// #include <seqan3/alphabet/views/all.hpp>        // include all of SeqAn's views 
// #include <seqan3/alphabet/hash.hpp>
// #include <seqan3/alphabet/nucleotide/dna15.hpp> // seqan3::gapped<seqan3::dna15>
//#include <seqan3/alphabet/range/hash.hpp>
#include <seqan3/core/debug_stream.hpp>         
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/argument_parser/all.hpp> 
// #include <seqan3/alignment/all.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "seq-cov.hpp"

// todo grap most std::cout and seqan3::debug_stream with  if constexpr (debugging) { }


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

 bool PCRSplitter::quick_check(const oligo_seq_t &target, const oligo &primer, long offset) // check if the primer matches the target at the expected position
 {
    long pr_beg = primer.ref_beg + offset;
    long len = primer.seq.size();
    if (pr_beg < 0 || pr_beg + len > target.size()) return false;             // primer is out of the current target... but (pos + ref_len < full_target.size())
    int mm = len - primer.match;                                              // mm - not allowed number of mismatches
    if (mm < 0) mm = 0;                                   // throw std::runtime_error{"Primer " + primer.name + " has incorrect number of permisible mismatches: " 
    for (int i = 0; i < len - mm; ++i)        // dont check the last, allowed mismatches
    {
        if (mismatch.score(primer.seq[i], target[pr_beg + i]))   
            if (! mm-- ) return false;                           // already consumed all allowed mismatches
    }
    return true;
 }

// struct SplitGene
 
PCRSplitter::PCRSplitter(SplitCoVfasta const &parent, std::string pcr_name)
        : parent{parent}, pcr_name{pcr_name}
{
    
    //if (split) file_E.open( );  // todo implement conditional split
    std::cout << std::boolalpha  
        << "\nGoing to split: " << (parent.dir / parent.fasta_name).string() 
        << ", PCR " << pcr_name << " into "
        << (parent.dir / (pcr_name + "." + parent.fasta_name)).string();
    
}

bool PCRSplitter::read_oligos(const std::filesystem::path& path_oligos)
{
    // just once - no need to be super efficient    
    if (path_oligos.empty()) return false;   // todo: more checks?

        if constexpr (debugging >= debugging_INFO) 
        {  seqan3::debug_stream << "\nGoing to load: " << path_oligos.string(); }
        seqan3::sequence_file_input<OLIGO> file_in{path_oligos};
        std::string fw, rv;
        ref_beg = std::numeric_limits<long>::max();        // agresive initialization: but the whole class is create every time the user click 'split'
        ref_end = std::numeric_limits<long>::min();
        f_primers.clear();
        r_primers.clear();
        probes_s.clear();
        probes_a.clear();
        //extern_forw_idx = extern_rev_idx = -1;
        all_oligos.clear();

        for (auto & primer : file_in)
        {
            std::string id{primer.id()};
            if constexpr (debugging >= debugging_VERBOSE)    seqan3::debug_stream << "\nGoing to check: " << id << "\n" << primer.sequence();

            // E_Sarbeco_F2 -13814 MN908947.3: Seq pos: 26308-26329, PosHint1=[26196	26208	26220	26219], PosHint2=[26171	26300], PosHint3=[26162	25995	25920	26168], PosHint4=[25116	26376]
            
            oligo pr;

            std::stringstream ss{id};
            ss >> pr.name >> pr.code ;                                 // Parse the primer name and code (first two tokens)
            
            std::size_t posPos = id.find("Seq pos: ");                 // Find and parse the "Seq pos: " field (e.g. "Seq pos: 26308-26329,")
            if (posPos != std::string::npos)
            {
                std::string seqPos = id.substr(posPos + 9);            // skip "Seq pos: "
                std::size_t commaPos = seqPos.find(',');               // Trim at the first comma in case hints follow immediately.
                if (commaPos != std::string::npos)   seqPos = seqPos.substr(0, commaPos);
                std::size_t dashPos = seqPos.find('-');
                if (dashPos != std::string::npos)
                {
                    pr.ref_beg = std::stoi(seqPos.substr(0, dashPos));
                    pr.ref_end = std::stoi(seqPos.substr(dashPos + 1));
                }
            }                                                          // for check if beg, end are valid see below
            if constexpr (debugging >= debugging_VERBOSE)    seqan3::debug_stream << " with beg: " << pr.ref_beg << " and end: " << pr.ref_end;

            // Parse optional positional hints using a regex. Expected patterns:
            //   PosHint1=[val ...]  -> vector<int>
            //   PosHint2=[beg end]  -> SeqPos (hint2)
            //   PosHint3=[val ...]  -> vector<int>
            //   PosHint4=[beg end]  -> SeqPos (hint4)
            std::regex hintRegex(R"(PosHint(\d)=\[(.*?)\])");                           // ask OpenAI o3
            auto hintBegin = std::sregex_iterator(id.begin(), id.end(), hintRegex);
            auto hintEnd   = std::sregex_iterator();
            for (auto it = hintBegin; it != hintEnd; ++it)
            {
                std::smatch match = *it;
                int hintIndex = std::stoi(match[1].str());
                std::string hintStr = match[2].str();
               
                std::istringstream hintStream(hintStr);        // Tokenize by any whitespace.
                std::vector<int> tokens;
                int num;
                while (hintStream >> num)    tokens.push_back(num);
                
                switch(hintIndex)
                {
                    case 1:  if (hint1.empty())           hint1 = tokens;                        break;

                    case 2:  if (tokens.size() == 2)           // Expect two tokens: begin and end values.
                        {
                            hint2.beg = tokens[0];
                            hint2.end = tokens[1];
                        } 
                        else if (debugging >= debugging_ERROR) seqan3::debug_stream << "Invalid hint2: "  << hintStr <<
                              "(please set to beg \t end of most frecuent amplicon region, like in PosHint2=[26171	26300]):  "<< '\n';
                        break;

                    case 3:  if (hint3.empty())    hint3 = tokens;                             break;
                    
                    case 4:  if (tokens.size() == 2)          // Expect two tokens: begin and end values.
                        {
                            hint4.beg = tokens[0];
                            hint4.end = tokens[1];
                        }
                        else if (debugging >= debugging_ERROR) seqan3::debug_stream << "Invalid hint4: "  << hintStr <<
                              "(please set to beg \t end of region with almost all the rest of the amplicon positions, like in PosHint4=[25116	26376]):  "<< '\n';
                        break;
                    
                    default:                        break;
                }
            }
            pr.seq   = primer.sequence();
            pr.len   = pr.seq.size();
            pr.match = std::round(parent.match * double(pr.len) / 100.0);
            if constexpr (debugging >= debugging_VERBOSE)  seqan3::debug_stream << " with minimum matches:" << pr.match;

            if (pr.ref_beg < pr.ref_end)                    // update this PCR-amplicon
            {
                pr.reverse = false;

                // throw std::runtime_error if incorrect lenth of the primer
                if (pr.len != (pr.ref_end - pr.ref_beg + 1)) 
                    throw std::runtime_error{"Forward primer " + pr.name + " has incorrect length: " + std::to_string(pr.len) +
                    " instead of " + std::to_string(pr.ref_end - pr.ref_beg + 1) + " at position " + std::to_string(pr.ref_beg) 
                    + " to " + std::to_string(pr.ref_end)};

                if ( ref_beg > pr.ref_beg)               // current extern forward primer
                {
                    ref_beg = pr.ref_beg;
                    extern_forw_idx = f_primers.size();
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

                if (ref_end < pr.ref_beg)                          // extern reverse primer
                {
                    ref_end = pr.ref_beg;
                    extern_rev_idx = r_primers.size();
                }
                r_primers.push_back(pr);
                all_oligos.push_back(pr);    
                continue;
            }
        }
 
        if constexpr (debugging >= debugging_INFO)        // logg info about hints
        {
            seqan3::debug_stream << "\nHints for " << pcr_name << " PCR: ";

            seqan3::debug_stream << "\nHint1: ";
            for (auto & h : hint1) seqan3::debug_stream << h << " ";

            seqan3::debug_stream << "\nHint2: " << hint2.beg << " " << hint2.end;

            seqan3::debug_stream << "\nHint3: ";
            for (auto & h : hint3) seqan3::debug_stream << h << " ";

            seqan3::debug_stream << "\nHint4: " << hint4.beg << " " << hint4.end << std::endl;
        }  

        ref_len = ref_end - ref_beg + 1;
        if (ref_len < 40) 
            throw std::runtime_error{"Amplicon length is too short: " + std::to_string(ref_len)};

        /*for (auto & primer : f_primers) all_oligos.push_back(&primer);
        for (auto & primer : r_primers) all_oligos.push_back(&primer);
        for (auto & primer : probes_s) all_oligos.push_back(&primer);
        for (auto & primer : probes_a) all_oligos.push_back(&primer);*/

        return true;
}

// eliminate gaps and put the sequencce into the oligo_seq_t seq
bool PCRSplitter::reconstruct_msa_seq(const msa_seq_t& full_msa_seq, oligo_seq_t& seq, long msa_beg, long msa_end, int tent_len)
{
    if constexpr (debugging >= debugging_TRACE+3) 
                seqan3::debug_stream << "\nReconstructing sequence from MSA positions " << msa_beg << " to " << msa_end << '\n';
    
    seq.clear();
    seq.reserve(tent_len);
    bool reverse = (msa_beg > msa_end);
    
    if (reverse)  // original sequence is in reverse order
    {       
        // go through the sequence and eliminate gaps to reconstruct the original sequence
        for ( int i = msa_beg; i >= msa_end; --i)
        {
            if (full_msa_seq[i].holds_alternative<seqan3::gap>() ) continue;
            seq.push_back(full_msa_seq[i].convert_unsafely_to<oligo_seq_alph>().complement());
        }
        return true;
    }
    // original sequence is in the same order
    // go through the sequence and eliminate gaps to reconstruct the original sequence
    for ( int i = msa_beg; i <= msa_end; ++i)
    {
        if (full_msa_seq[i].holds_alternative<seqan3::gap>() ) continue;
        seq.push_back(full_msa_seq[i].convert_unsafely_to<oligo_seq_alph>());
    }
    return true;
}

long find_primer(const oligo_seq_t& target, const oligo& primer, auto& config)
{
    auto results = seqan3::align_pairwise(std::tie(target, primer.seq), config);
    if constexpr (debugging >= debugging_TRACE) seqan3::debug_stream << "check if there are results" << '\n';
 
    auto res_beg = results.begin();
    if (res_beg == results.end())    return -1;
    
    auto & res = *res_beg;
    if constexpr (debugging >= debugging_TRACE+3) 
    {seqan3::debug_stream << "Alignment: " << res.alignment() << " Score: "     << res.score() ;
    seqan3::debug_stream << ", Target: (" << res.sequence1_begin_position() << "," << res.sequence1_end_position() << ")";
    seqan3::debug_stream << ", Primer: (" << res.sequence2_begin_position() << "," << res.sequence2_end_position() << "). " ;}
    
    if (res.score() < primer.match)    return -1;
    
    return res.sequence1_begin_position() - res.sequence2_begin_position() ;  // target begin position - primer begin position
    
} 


// find one of the external primers in the sequence target and return the position of the amplicon
SeqPos PCRSplitter::find_ampl_pos(const oligo_seq_t& target)
{
    // new, unknown seq. We need to find the right position of the target sequence
    SeqPos sg;
    oligo & fw_pr = f_primers[extern_forw_idx];
    
    auto output_config =    seqan3::align_cfg::output_score{}          | 
                            seqan3::align_cfg::output_begin_position{} |
                            seqan3::align_cfg::output_end_position{}   ;
    auto method = seqan3::align_cfg::method_local{};
    seqan3::align_cfg::scoring_scheme     scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score   { 1}, 
                                                                                   seqan3::mismatch_score{-1}}};
    seqan3::align_cfg::gap_cost_affine gap_costs{seqan3::align_cfg::open_score{0}, 
                                                 seqan3::align_cfg::extension_score{-1}};
    auto config = method | scheme | gap_costs | output_config;
    auto filter_v = std::views::filter( [&](auto && res) { return res.score() >= fw_pr.match; });

    if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "try to align the external fw primer to target: " << '\n';

    auto results = seqan3::align_pairwise(std::tie(target, fw_pr.seq), config);

    for (auto & res : results | filter_v)
    {
        if constexpr (debugging >= debugging_TRACE+3) 
        {
        seqan3::debug_stream << /*"Alignment: " << res.alignment() << */" Score: "     << res.score() ;
        seqan3::debug_stream << ", Target: (" << res.sequence1_begin_position() << "," << res.sequence1_end_position() << ")";
        seqan3::debug_stream << ", Primer: (" << res.sequence2_begin_position() << "," << res.sequence2_end_position() << "). " ;}
               
        sg.beg = res.sequence1_begin_position() - res.sequence2_begin_position() ;  // target begin position - primer begin position
        if (ref_len)  // not the first time/seq - the "ref." seq was already set
        {
            sg.end = sg.beg + ref_len;
            return sg;                     // only find fw pr pos. Don't count for insertions/deletions or bad rev.
        }
        else 
        if (sg.beg < 0)  // the primer match before this "reference" target begin (we are asuming the first seq is the reference)
                         // becouse it is used to set ref_len 
            throw std::runtime_error{"First " + pcr_name + " sequence don't contain "
            "the full forward primer. Score: " + std::to_string(res.score()) +
            " that begin at position "  + std::to_string(res.sequence2_begin_position())  };
        break;
    }
    if (!ref_len && sg.beg == sg.npos) // the fw primer was not found in the target and we don't have the ref_len
        throw std::runtime_error{"First " + pcr_name + " sequence don't contain the forward primer but we want it to be the ref." };

    if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "try to align the external rev primer to target: " << '\n';

    auto & rev_pr = r_primers[extern_rev_idx];
    auto res_rev = seqan3::align_pairwise(std::tie(target, rev_pr.seq), config);
    auto rfilter_v = std::views::filter( [&](auto && res) { return res.score() >= rev_pr.match; });
    for (auto & res_r : res_rev | rfilter_v)
    {
        if constexpr (debugging >= debugging_TRACE+3) 
        {seqan3::debug_stream << /*"\nAlignment: " << res_r.alignment() <<  */" Score: "<< res_r.score() ;
        seqan3::debug_stream << ", Target: ("     << res_r.sequence1_begin_position() << "," << res_r.sequence1_end_position() << ")";
        seqan3::debug_stream << ", rev Primer: (" << res_r.sequence2_begin_position() << "," << res_r.sequence2_end_position() << "). ";}

        sg.end = res_r.sequence1_end_position() + (rev_pr.seq.size() - res_r.sequence2_end_position()); // target begin position + primer length - primer begin position
        if (ref_len)  // not the first time/seq - the "ref." seq was already set
            { 
                if ( sg.beg == sg.npos) sg.beg = sg.end - ref_len; // only find rev pr pos. Don't count for insertions/deletions or bad rev.
                return sg;
            }  
        else if (sg.end > target.size())
            throw std::runtime_error{"First " + pcr_name + " sequence don't contain "
            "the full reverse primer. Score: " + std::to_string(res_r.score()) +
            " that end at position "  + std::to_string(res_r.sequence2_end_position())  };
        break;
    }
    if (sg.beg != sg.npos && sg.end != sg.npos) return sg;    // we have both primers but still not the ref_len
    sg.beg = sg.end = sg.npos;                                // mark as not found !! 
    if constexpr (debugging >= debugging_TRACE) seqan3::debug_stream << pcr_name + " PCR don't match the external primers in the Target region of length: " << target.size() << '\n';
    return sg;
}

void PCRSplitter::re_align(pattern_q &pq, oligo_seq_t &oligo_target)  // re_align the oligo to an exact target and build primer match pattern on oligo_target
{
    if constexpr (debugging >= debugging_TRACE) 
        seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                            << pq.primer.seq <<" - oligo seq\n" 
                            << oligo_target  <<" - target seq\n" ;

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
    if constexpr (debugging >= debugging_TRACE+3) 
                seqan3::debug_stream << "\nTarget: " << oligo_target << "\nPrimer: " << pq.primer.seq << '\n';

    if constexpr (debugging >= debugging_TRACE+3) 
                seqan3::debug_stream << " Going to check results\n";     

    for (auto const & res : seqan3::align_pairwise(std::tie(oligo_target,  // sequence1 = target
                                                            pq.primer.seq  // sequence2 = primer
                                                            ), config))
    //auto results = seqan3::align_pairwise(std::tie(target, pq.primer.seq), config);
    //for (auto & res : results)
    // if (res.score() > pq.primer.match)  // primer found. 
    {
        if constexpr (debugging >= debugging_TRACE) 
                seqan3::debug_stream << "Alignment: " << res.alignment() << " Score: "     << res.score() ;
        if constexpr (debugging >= debugging_TRACE) 
                seqan3::debug_stream << ", Target: (" << res.sequence1_begin_position() << "," << res.sequence1_end_position() << ")";
        if constexpr (debugging >= debugging_TRACE) 
                seqan3::debug_stream << ", Primer: (" << res.sequence2_begin_position() << "," << res.sequence2_end_position() << "). " ;
    
            //brief Helper function to print elements of a tuple separately: aligned_target, aligned_primer
            auto&[aligned_target, aligned_primer] = res.alignment();
            pq.pattern.clear();
            
            int len_a_pr = aligned_primer.size();
            int len_a_tg = aligned_target.size(); 
            int len = std::max(len_a_pr, len_a_tg);

            if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << " cont...\n" ;

            int o_p = res.sequence2_begin_position();  // add needed '-' to the pattern before aligned region ??
            if (o_p) pq.pattern = std::string(o_p, '-');

            for (int i = 0; i < len; ++i)  // go through the aligned region, len = std::max(len_a_pr, len_a_tg)
            {
                if (i >= len_a_pr)  // extra nt in aligned target ?? ignore ?
                {
                    if (o_p >= pq.primer.len) continue  ;  // ??
                    o_p++;
                    pq.pattern.push_back(aligned_target[i].to_char());   
                    if  (aligned_target[i] == 'N'_dna15)  
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
                if (aligned_primer[i].holds_alternative<seqan3::gap>() &&
                    aligned_target[i].holds_alternative<seqan3::gap>()    )       continue;  // ??

                if (aligned_target[i].holds_alternative<seqan3::gap>()    )     // deletion in target (or insertion in ref?)
                {
                    if (o_p >= pq.primer.len) continue  ;  // ??
                    o_p++;
                    pq.pattern.push_back('-');  
                    continue;
                }
                auto t = aligned_target[i].convert_unsafely_to<oligo_seq_alph>();        
                if (aligned_primer[i].holds_alternative<seqan3::gap>())  // insertion in target 
                {
                    if (!o_p) continue;  // ignore insertions before primer
                    // assume ???? that the inserted base is not in the primer
                    pq.pattern.push_back(t.to_char());  // inserted in the pattern but not advance the primer
                    pq.mm++;
                    if (pq.primer.len - o_p <= parent.crit_term_nt)                   pq.crit++;  // ??
                    continue;
                }
                auto p = aligned_primer[i].convert_unsafely_to<oligo_seq_alph>();   
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

            if constexpr (debugging >= debugging_TRACE)
            seqan3::debug_stream << aligned_primer        << " - aligned oligo\n"
                                 << pq.pattern            << " - pattern\n"
                                 << aligned_target        << " - aligned target\n" ;
            return;
    }
    pq.Q = 1000;
    auto tostring = oligo_target | seqan3::views::to_char;
    pq.pattern = std::string{tostring.begin(), tostring.end()};

    if constexpr (debugging >= debugging_TRACE)
    seqan3::debug_stream << pq.primer.seq << " - unaligned oligo\n"
                         << pq.pattern    << " - unaligned pattern\n"
                         << oligo_target  << " - unaligned target\n" ;
}

void PCRSplitter::align_to_msa(pattern_q &pq,      ///< oligo_pattern_quality  // build primer match pattern on oligo_target on to the expected MSA position. - not used !!
                    const msa_seq_t &full_target)
{
    oligo_seq_t oligo_target;
    msa_seq_t   msa_fragment;
    cov::oligo &primer = pq.primer;
    parent.extract_msa_seq(full_target, msa_fragment, oligo_target, primer.msa_beg, primer.msa_end, primer.seq.size());

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

void PCRSplitter::target_msa_pattern(target_q & tq, const msa_seq_t& full_target) // build target match pattern on ref on to the expected MSA positions msa_beg, msa_end.
{
    // msa_seq_t& full_target is already aligned to the reference sequence: just create the target pattern
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

void PCRSplitter::evaluate_msa_target(target_q & tq, const msa_seq_t& full_target)
{
    for (oligo& primer : all_oligos)   
    {
        tq.patterns.emplace_back(primer);  // registering/creating the pattern_q is cheap and fast but difficult to parallelize
    }
    //for (pattern_q& pq : tq.patterns)   
    std::for_each(std::execution::par_unseq, tq.patterns.begin(), tq.patterns.end(), [&](pattern_q &pq)
    {
        evaluate_msa_target_primer(pq, full_target);
    });
    target_msa_pattern(tq, full_target);

    for (auto & pq : tq.patterns)
        if constexpr (debugging >= debugging_DEBUG)
        seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             //<< target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 }

void PCRSplitter::evaluate_msa_target_primer(pattern_q &pq, const msa_seq_t& full_target) // build primer match pattern on full_target on to the expected MSA positions msa_beg, msa_end of the primer.
{
    cov::oligo &primer = pq.primer;
    oligo_seq_t oligo_target;

    reconstruct_msa_seq(full_target, oligo_target, primer.msa_beg, primer.msa_end, primer.seq.size());
    if (oligo_target.size() != primer.len)
        //return align(pq, full_target); // todo: try to align the primer with the oligo_target sequence
        return re_align(pq, oligo_target); // re align the primer with the oligo_target sequence
 
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
        seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             << oligo_target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 
}

// Some functions are now split into two versions: one for the old MSA format and a new version for the new 'raw' sequence format. 
// As the sequences are now not pre-aligned we don't know the exact position 
// of the amplicon and of each oligo/primer. 
// We are saving and grouping the possible matching patterns by equal target sequences member map grouped. 
// Thus, for each new sequence we have the following cases: 
// 1-exact known position and sequence: just update the counters. 
// 2-right position but new sequence: 
//         simple check that one of external primers matches the expected position: 
//          reuse the newly retrieved parsed_id pid to hold all the patterns, that have to be generated first. 
// 3-New position: 
//     first discard wrong seq/pid 
//     them we need to align one of the external primers to the expected target first 
//         (if not found expand by 1000 that region, 
//             if fails use the whole sequence, 
//             and try inverted too) 
//     and readjust the coordinates of the amplicon and repeat to check if we are now in 1- or 2-. 
std::optional<std::reference_wrapper<target_count>> PCRSplitter::check_rec(const auto& record)
{
    oligo_seq_t const & full_target = record.sequence();
    count++;
    oligo& fw_pr = f_primers[extern_forw_idx];
    // high_resolution_clock clock start
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();


    for (long pos : hint1)                                              // first try hint1: most frecuent single amplicon positions with quick check
    {
        if (pos + ref_len >= full_target.size()) continue;               // out of the sequence
        if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << '\n' << pcr_name << ": quick_check the forward primer at the Hint1 position: " << pos << "\n";
        if (!quick_check(full_target, fw_pr, pos - ref_beg)) continue;
        
        if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "\nFound the forward primer at the Hint1 position: " << pos << "\n";
        if constexpr (debugging >= debugging_NOT_USED)   amplicon_pos_beg[pos]++;
        count_hint1++;

        const auto& [it, is_new_seq] = grouped.try_emplace({full_target.begin()+pos, 
                                                            full_target.begin()+pos + ref_len}, target_count{});  
        target_count & target_c = it->second;
        if (is_new_seq) evaluate_target(target_c.target, full_target, pos); // 2-known position but new sequence
        target_c.count++;
        // add high_resolution duration to used time  
        hint1_time_used += std::chrono::high_resolution_clock::now() - start;
        return std::ref(target_c);          
    }
    // add high_resolution duration to wasted time 
    std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
    hint1_time_wasted += finish - start;
    start = finish;

    // return std::nullopt;    

    if (full_target.size() > hint2.beg + ref_len) try                                         // found the right position inside the Hint2 region
    {
        long end = std::min<long>(full_target.size(), hint2.end);

        SeqPos sg = find_ampl_pos({full_target.begin() + hint2.beg, 
                                   full_target.begin() + end      });
        
        if (sg.beg != sg.npos && sg.end != sg.npos && sg.beg < sg.end)            // found the right position inside the Hint2 region
        {
            sg.beg = std::max<long>(0L                , hint2.beg + sg.beg);
            sg.end = std::min<long>(full_target.size(), hint2.beg + sg.end);
            if (ref_len <= sg.end - sg.beg) // the amplicon is long enough
            {
                if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "\nFound the amplicon in Hint2 region position: " << sg.beg << " to " << sg.end << std::endl;
                if constexpr (debugging >= debugging_NOT_USED) amplicon_pos_beg[sg.beg]++;
                count_hint2++;

                const auto& [it, is_new_seq] = grouped.try_emplace({full_target.begin()+sg.beg,
                                                                    full_target.begin()+sg.end}, target_count{});
                target_count & target_c = it->second;
                if (is_new_seq) evaluate_target(target_c.target, full_target, sg.beg); // New position and new sequence
                target_c.count++;
                hint2_time_used += std::chrono::high_resolution_clock::now() - start;
                return std::ref(target_c);     
            }
        }
        if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "No amplicon found in Hint2 region\n";
    } 
    catch (std::exception & e) 
    {
        if constexpr (debugging >= debugging_DEBUG)  seqan3::debug_stream << "ERROR: No amplicon found, becouse Error: " << e.what() 
          << " checking new target sequence: " << record.id() << " with seq:\n" << full_target <<  "\n";
    }
    finish = std::chrono::high_resolution_clock::now();
    hint2_time_wasted += finish - start;
    start = finish;

    // return std::nullopt;    

    for (long pos : hint3)                                                          // try hint3: less frecuent single amplicon positions with quick check
    {
        if (pos + ref_len >= full_target.size()) continue;  // out of the sequence
        if (!quick_check(full_target, fw_pr, pos - ref_beg)) continue;
        
        if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "\nFound the forward primer at the Hint3 position: " << pos << "\n";
        if constexpr (debugging >= debugging_NOT_USED)     amplicon_pos_beg[pos]++;
        count_hint3++;

        const auto& [it, is_new_seq] = grouped.try_emplace({full_target.begin()+pos, 
                                                            full_target.begin()+pos + ref_len}, target_count{}); 
        target_count & target_c = it->second;
        if (is_new_seq) evaluate_target(target_c.target, full_target, pos);                            // 2-known position but new sequence
        target_c.count++;
        hint3_time_used += std::chrono::high_resolution_clock::now() - start;
        return std::ref(target_c);          
    }
    finish = std::chrono::high_resolution_clock::now();
    hint3_time_wasted += finish - start;
    start = finish;

    // return std::nullopt;    

    if (full_target.size() > hint4.beg + ref_len) try                                                  // found the right position inside the Hint4 region
    {      
        long end = std::min<long>(full_target.size(), hint4.end);
        SeqPos sg = find_ampl_pos({full_target.begin()+hint4.beg, 
                                   full_target.begin()+end});
        
        if (sg.beg != sg.npos && sg.end != sg.npos && sg.beg < sg.end)                           // found the right position inside the Hint4 region
        {
            sg.beg = std::max<long>(0L                , hint4.beg + sg.beg);
            sg.end = std::min<long>(full_target.size(), hint4.beg + sg.end);
            if (ref_len <= sg.end - sg.beg) // the amplicon is long enough
            {
                if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "\nFound the amplicon in Hint4 region position: " << sg.beg << " to " << sg.end << "\n";
                if constexpr (debugging >= debugging_NOT_USED) amplicon_pos_beg[sg.beg]++;
                count_hint4++;

                const auto& [it, is_new_seq] = grouped.try_emplace({full_target.begin()+sg.beg,
                                                                    full_target.begin()+sg.end}, target_count{});
                target_count & target_c = it->second;
                if (is_new_seq) evaluate_target(target_c.target, full_target, sg.beg);               // New position and new sequence
                target_c.count++;
                hint4_time_used += std::chrono::high_resolution_clock::now() - start;
                return std::ref(target_c);   
            }   
        }
        if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "No amplicon found in Hint4 region\n";
    } 
    catch (std::exception & e) 
    {
        if constexpr (debugging >= debugging_DEBUG)  seqan3::debug_stream << "ERROR: No amplicon found, becouse Error: " << e.what() 
            << " checking new target sequence: " << record.id() << " with seq:\n" << full_target <<  "\n";
    }
    finish = std::chrono::high_resolution_clock::now();
    hint4_time_wasted += finish - start;
    start = finish;
    // return std::nullopt;    

    try                                                                                       // search full_target
    {
        SeqPos sg = find_ampl_pos(full_target);
        
        if (sg.beg != sg.npos && sg.end != sg.npos && sg.beg < sg.end)  // found the right position 
        {
            sg.beg = std::max<long>(0L                , sg.beg);
            sg.end = std::min<long>(full_target.size(), sg.end);
            if (ref_len <= sg.end - sg.beg) // the amplicon is long enough, 
            {
                if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "\nFound the amplicon in Hint4 region position: " << sg.beg << " to " << sg.end << "\n";
                if constexpr (debugging >= debugging_INFO) amplicon_pos_beg[sg.beg]++;
                count_full++;

                const auto& [it, is_new_seq] = grouped.try_emplace({full_target.begin()+sg.beg,
                                                                    full_target.begin()+sg.end}, target_count{});
                target_count & target_c = it->second;
                if (is_new_seq) evaluate_target(target_c.target, full_target, sg.beg);     // New position and new sequence
                target_c.count++;
                full_seq_time_used += std::chrono::high_resolution_clock::now() - start;
                return std::ref(target_c);      
            }
        }
    } 
    catch (std::exception & e) 
    {
        if constexpr (debugging >= debugging_DEBUG)  seqan3::debug_stream << "ERROR: No amplicon found, becouse Error: " << e.what() 
            << " checking new target sequence: " << record.id() << " with seq:\n" << full_target <<  "\n";
    }
    finish = std::chrono::high_resolution_clock::now();
    full_seq_time_wasted += finish - start;
    if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "No amplicon found in this target\n";
// todo: try to align the reference amplicon to the full_target sequence
    // create a counter for the not found sequences gruoped every 1000 nt lenth  not_found[full_target.size()/1000]++;
    count_not_found++;
not_found_lenghts[full_target.size()/1000]++;
    return std::nullopt;
}

void PCRSplitter::evaluate_target(target_q  &tq, const oligo_seq_t &full_target, long ampl_beg)
{
    int offset = ampl_beg - ref_beg;
    for (oligo& primer : all_oligos)   
    {
        tq.patterns.emplace_back(primer);  // registering/creating the pattern_q is cheap and fast but difficult to parallelize
    }
    if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "\nGoing to evaluate a new target sequence for " << tq.patterns.size() << " primer pattenrs\n";
    //for (pattern_q& pq : tq.patterns)   
    std::for_each(std::execution::par_unseq, tq.patterns.begin(), tq.patterns.end(), [&](pattern_q &pq)
    {
        evaluate_target_primer(pq, full_target, offset);
    }
    );
    if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "\nGoing to evaluate the target sequence self " << tq.patterns.size() << '\n';
    target_pattern(tq, full_target, ampl_beg);

    for (auto & pq : tq.patterns)
        if constexpr (debugging >= debugging_DEBUG)
        seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             //<< target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 }

void PCRSplitter::evaluate_target_primer(pattern_q &pq, const oligo_seq_t &full_target, long offset) // build primer match pattern on full_target at positions beg
{
    cov::oligo &primer = pq.primer;
    pq.pattern = std::string(primer.len, '.');
    if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ": "  << pq.primer.seq <<'\n'  ;

    int pr_beg = primer.ref_beg + offset ; //   pr_end = primer.ref_end + offset;
    
    for (int i = 0; i < primer.len; ++i)
    { 
        // get the target nt at the position i taking into consideration the primer is forward or reverse
        seqan3::dna15 t{primer.reverse ? full_target[pr_beg - i].complement() 
                                       : full_target[pr_beg + i]};

        if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "Checking position " << i << ": " << primer.seq[i] << " vs " << t << '\n';
        if (!mismatch.score(primer.seq[i], t))   continue;
        if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "Mismatch found " << '\n';

        if  (t == 'N'_dna15)  
        {
            pq.N++;
            pq.pattern[i] = 'N';
            continue;
        }
        pq.pattern[i] = t.to_char();
        pq.mm++;
        if (i >= primer.len - parent.crit_term_nt)                   pq.crit++;
     }
    pq.Q = pq.mm + pq.crit * 4;
    if constexpr (debugging >= debugging_TRACE)  
        seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             //<< oligo_target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
}

void PCRSplitter::target_pattern(target_q & tq, const oligo_seq_t &full_target, long ampl_beg)// build target match pattern on ref at ampl_beg
{
    // already aligned to the reference sequence: just create the target pattern
    tq.target_pattern = std::string(ref_len, '.');

    for ( int i = 0; i < ref_len; ++i)
    {
        auto s =    full_target[ ampl_beg + i];
        auto r = parent.ref_seq[ ref_beg + i - 1];
        if (mismatch.score(s, r))   // 0 if equal
           tq.target_pattern[i] = s.to_char();   
    }
}

// sequences pre-aligned into the MSA to the reference sequence
target_count& PCRSplitter::check_msa_rec(auto& record)
{
    msa_seq_t& full_target = record.sequence();
    count++;

    target_count & target_c = msa_grouped[{full_target.begin()+msa_beg, 
                                           full_target.begin()+msa_end}]; 
    if (!target_c.count)  // new target sequence
        evaluate_msa_target(target_c.target, full_target);
    
    target_c.count++;  // assume always found
    return target_c;   
}

void build_id_field_patterns(target_count & target_c, std::string & patterns)
{
        patterns.clear();              

        for (auto& pq : target_c.target.patterns)
        {
            patterns += std::format("|{}_Q_{}_mm_{}_N_{}_crit_{}_pat_{}", 
                            pq.primer.name, pq.Q, pq.mm, pq.N, pq.crit, pq.pattern);
        }
        patterns += "|" + target_c.target.target_pattern + "|" ;  // at the end we will add the count    
}

void PCRSplitter::write_grouped ()
{
    using types    = seqan3::type_list<oligo_seq_t, std::string>;
    using fields   = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using record_t = seqan3::sequence_record<types, fields>;
    using sgr_t    = decltype(grouped)::value_type;  // std::unordered_map<oligo_seq_t, target_count>;
    using sgr_p    = sgr_t*;
    using fasta_out= seqan3::sequence_file_output<fields, seqan3::type_list<seqan3::format_fasta> >;

       auto name = parent.fasta_name + "." + pcr_name;

       std::filesystem::path gr  = parent.dir / (name + ".grouped.fasta" );

       std::filesystem::path grm    = parent.dir / (name + ".monthly.fasta" );
       //std::filesystem::path grmc   = parent.dir / (name + ".monthly-continent.fasta" );
       //std::filesystem::path grmcc  = parent.dir / (name + ".monthly-continent-country.fasta" );
       //std::filesystem::path grmccc = parent.dir / (name + ".monthly-continent-country-clade.fasta" );

       std::filesystem::path grd    = parent.dir / (name + ".daily.fasta"   );
       std::filesystem::path grdc   = parent.dir / (name + ".daily-continent.fasta" );
       std::filesystem::path grdcc  = parent.dir / (name + ".daily-continent-country.fasta" );
       std::filesystem::path grdccc = parent.dir / (name + ".daily-continent-country-clade.fasta" );

    seqan3::debug_stream << "Going to write: " << grdccc << "\n" ;
    fasta_out  file_e_gr  {gr  },  
               file_e_grm {grm }; //, file_e_grmc {grmc}, file_e_grmcc{grmcc}, file_e_grmccc{grmccc};
    fasta_out  file_e_grd {grd }, file_e_grdc {grdc}, file_e_grdcc{grdcc}, file_e_grdccc{grdccc};
    
    std::vector<sgr_p> gr_v;
    gr_v.reserve(grouped.size());
    for (sgr_t& sgr : grouped) gr_v.push_back(&sgr);  // vector of pointers to msa_grouped target_count sequences

    std::sort(gr_v.begin(), gr_v.end(), []( sgr_p &a,  sgr_p &b)    {return a->second.count > b->second.count;});

    std::string patterns, id;
    patterns.reserve(4000);
    id.reserve(4000);
    // count and logg ignored sequences due to too many Ns and to empty patterns
    long many_N = 0, empty_patterns = 0;

    for (sgr_p           sg :  gr_v            )   // pointers to msa_grouped target_count   seq: target_count
    {
        target_count & target_c = sg->second;

        // discard sequences with too many N.
        if (std::all_of(target_c.target.patterns.begin(), target_c.target.patterns.end(), 
                        [&](pattern_q &pq) {return pq.N > parent.crit_N;}))        { many_N++; continue;}

        build_id_field_patterns(target_c, patterns);
        if (patterns.empty())                                    {empty_patterns++; continue;}

        for (auto& [year,      yc]:  sg->second.years) 
        for (auto& [month,     mc]:  yc.months       ) {
        for (auto& [day,       dc]:  mc.days         ) {
        for (auto& [continent, ct]:  dc.continents   ) {
        for (auto& [country,   cc]:  ct.countries    )   
        {    
            for (auto& [clade, cl]:  cc.clades       ) 
            {
                if (!cl.id) continue;  // ERROR !!!!!
                auto & pid = *cl.id;

                id = std::format(
                    "|{} |{:04d}-{:02d}-{:02d}|{}|{}|{}|{}|{}|{}|{}",
                    pid.EPI_ISL, 
                    year, month, day,
                    pid.isolate, 
                    pid.continent,   // or combine them 
                    pid.country,
                    pid.region,
                    pid.clade,
                    pid.pango,
                    pid.pango_version    ) + patterns;

                // output to ".daily-continent-country-clade" file_e_grdccc
                file_e_grdccc.push_back<record_t>( record_t{sg->first, std::to_string(cl.count) + id} );
            }
            // output to ".daily-continent-country" file_e_grdcc
            file_e_grdcc.push_back<record_t>( record_t{sg->first, std::to_string(cc.count) + id} );
        }
        // output to ".daily-continent" file_e_grdc
        file_e_grdc.push_back<record_t>( record_t{sg->first, std::to_string(ct.count) + id} );
        }
        // output to ".daily" file_e_grd
        file_e_grd.push_back<record_t>( record_t{sg->first, std::to_string(dc.count) + id} );
        }
        // output to ".monthly" file_e_grm
        file_e_grm.push_back<record_t>( record_t{sg->first, std::to_string(mc.count) + id} );
        }
        // output to ".grouped" file_e_gr
        file_e_gr.push_back<record_t>( record_t{sg->first, std::to_string(sg->second.count) + id} );
    }
    // logg ignored sequences due to too many Ns and to empty patterns
    if constexpr (debugging >= debugging_INFO) 
        seqan3::debug_stream << name << ": Ignored sequences due to too many Ns: " << many_N << " and to empty patterns: " << empty_patterns << '\n';
}

void PCRSplitter::write_msa_grouped ()
{
    using types    = seqan3::type_list<msa_seq_t, std::string>;
    using fields   = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using record_t = seqan3::sequence_record<types, fields>;
    using sgr_t    = decltype(msa_grouped)::value_type;  // std::unordered_map<msa_seq_t, target_count>;
    using sgr_p    = sgr_t*;
    using fasta_out= seqan3::sequence_file_output<fields, seqan3::type_list<seqan3::format_fasta> >;

       auto name = parent.fasta_name + "." + pcr_name;

       std::filesystem::path gr  = parent.dir / (name + ".grouped.fasta" );

       std::filesystem::path grm    = parent.dir / (name + ".monthly.fasta" );
       //std::filesystem::path grmc   = parent.dir / (name + ".monthly-continent.fasta" );
       //std::filesystem::path grmcc  = parent.dir / (name + ".monthly-continent-country.fasta" );
       //std::filesystem::path grmccc = parent.dir / (name + ".monthly-continent-country-clade.fasta" );

       std::filesystem::path grd    = parent.dir / (name + ".daily.fasta"   );
       std::filesystem::path grdc   = parent.dir / (name + ".daily-continent.fasta" );
       std::filesystem::path grdcc  = parent.dir / (name + ".daily-continent-country.fasta" );
       std::filesystem::path grdccc = parent.dir / (name + ".daily-continent-country-clade.fasta" );

    seqan3::debug_stream << "Going to write: " << grdccc << "\n" ;
    fasta_out  file_e_gr  {gr  },  
               file_e_grm {grm }; //, file_e_grmc {grmc}, file_e_grmcc{grmcc}, file_e_grmccc{grmccc};
    fasta_out  file_e_grd {grd }, file_e_grdc {grdc}, file_e_grdcc{grdcc}, file_e_grdccc{grdccc};
    
    std::vector<sgr_p> gr_v;
    gr_v.reserve(msa_grouped.size());
    for (sgr_t& sgr : msa_grouped) gr_v.push_back(&sgr);  // vector of pointers to msa_grouped target_count sequences

    std::sort(gr_v.begin(), gr_v.end(), []( sgr_p &a,  sgr_p &b)    {return a->second.count > b->second.count;});

    std::string patterns, id;
    patterns.reserve(4000);
    id.reserve(4000);

    for (sgr_p           sg :  gr_v            )   // pointers to msa_grouped target_count   seq: target_count
    {
        target_count & target_c = sg->second;
        
        // discard sequences with too many N.
        if (std::all_of(target_c.target.patterns.begin(), target_c.target.patterns.end(), 
                        [&](pattern_q &pq) {return pq.N > parent.crit_N;})) 
            return;

        build_id_field_patterns(target_c, patterns);
        if (patterns.empty()) continue;

        for (auto& [year,      yc]:  sg->second.years) 
        for (auto& [month,     mc]:  yc.months       ) {
        for (auto& [day,       dc]:  mc.days         ) {
        for (auto& [continent, ct]:  dc.continents   ) {
        for (auto& [country,   cc]:  ct.countries    )   
        {    
            for (auto& [clade, cl]:  cc.clades       ) 
            {
                if (!cl.id) continue;  // ERROR !!!!!
                auto & pid = *cl.id;

                id = std::format(
                    "|{} |{:04d}-{:02d}-{:02d}|{}|{}|{}|{}|{}|{}|{}",
                    pid.EPI_ISL, 
                    year, month, day,
                    pid.isolate, 
                    pid.continent,   // or combine them 
                    pid.country,
                    pid.region,
                    pid.clade,
                    pid.pango,
                    pid.pango_version    ) + patterns;

                // output to ".daily-continent-country-clade" file_e_grdccc
                file_e_grdccc.push_back<record_t>( record_t{sg->first, std::to_string(cl.count) + id} );
            }
            // output to ".daily-continent-country" file_e_grdcc
            file_e_grdcc.push_back<record_t>( record_t{sg->first, std::to_string(cc.count) + id} );
        }
        // output to ".daily-continent" file_e_grdc
        file_e_grdc.push_back<record_t>( record_t{sg->first, std::to_string(ct.count) + id} );
        }
        // output to ".daily" file_e_grd
        file_e_grd.push_back<record_t>( record_t{sg->first, std::to_string(dc.count) + id} );
        }
        // output to ".monthly" file_e_grm
        file_e_grm.push_back<record_t>( record_t{sg->first, std::to_string(mc.count) + id} );
        }
        // output to ".grouped" file_e_gr
        file_e_gr.push_back<record_t>( record_t{sg->first, std::to_string(sg->second.count) + id} );
    }
}

// class SplitCoVfasta

bool SplitCoVfasta::extract_msa_seq(const msa_seq_t &full_msa_seq, 
                                          msa_seq_t &msa_fragment, 
                                        oligo_seq_t &reconstructed_seq,
                        long msa_beg, long msa_end, int tent_len )  // = 0
{
    seqan3::debug_stream << "\nReconstructing sequence from MSA pos " << msa_beg << " to " << msa_end << '\n';
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

void SplitCoVfasta::set_ref_seq( )  // set reference sequence , todo: set/check positions of primers?
{
    seqan3::sequence_file_input<OLIGO> file_in{dir / "MN908947.3.fasta"};

    for (auto&& ref_rec: file_in)
   { 
        ref_seq = std::move(ref_rec.sequence());
        seqan3::debug_stream << "\n\nReference: " << ref_rec.id() << ", reference lenth = " << ref_seq.size() << '\n';
        return;
   }
}

void SplitCoVfasta::set_msa_ref_pos( )
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
    
    // for (auto & pcr : pcrs) pcr.set_msa_ref_pos(); Set/Check correct positions of primers on the reference sequence
    for (auto & pcr : pcrs) 
    {
        // set PCR target positions, first check ref_beg and ref_end were alrready set (in read_oligos)
        if (pcr.ref_beg == pcr.ref_end) // run time error if not set
            throw std::runtime_error{"PCR " + pcr.pcr_name + " has no reference positions set"};
        pcr.msa_beg = msa_pos[pcr.ref_beg-1];
        pcr.msa_end = msa_pos[pcr.ref_end-1];
        pcr.msa_len = pcr.msa_end - pcr.msa_beg + 1;
        
        // check all primers
        for (auto & primer : pcr.all_oligos)
        {
            if (primer.ref_beg == primer.ref_end) // run time error if not set
                throw std::runtime_error{"Primer " + primer.name + " has no reference positions set"};
            primer.msa_beg = msa_pos[primer.ref_beg-1];
            primer.msa_end = msa_pos[primer.ref_end-1];
            primer.msa_len = primer.msa_end - primer.msa_beg + 1;
            extract_msa_seq(msa_ref, primer.msa_ref, primer.ref_seq, primer.msa_beg, primer.msa_end, primer.seq.size());
        }

    }

}

void parsed_id::parse_id(const std::string& id) 
{
    // parse the id of the fasta record, like:
    // hCoV-19/USA/CO-CDPHE-2010040598/2020|2020-10-02|2021-09-04 
  try
  {
    std::size_t country_beg = 8;
    std::size_t country_end = id.find('/',   country_beg) - 1;          
    std::size_t country_len = country_end - country_beg   + 1;
    std::size_t isolate_beg = country_end                 + 2;
    std::size_t isolate_end = id.find('/', isolate_beg)   - 1;
    std::size_t isolate_len = isolate_end - isolate_beg   + 1;

    country = id.substr(country_beg, country_len);
    isolate = id.substr(isolate_beg, isolate_len);

    // parse the two dates and todo take the early as the collection date in pid.year, pid.month, pid.day
    std::size_t date_beg = id.find('|', isolate_end) + 1;
    std::size_t date_end = id.find('|', date_beg) - 1;
    std::size_t date_len = date_end - date_beg + 1;
    std::string date1 = id.substr(date_beg, date_len);
    
    // save the collection date in pid.year, pid.month, pid.day
    year  = std::stoi(date1.substr(0, 4));
    month = std::stoi(date1.substr(5, 2));
    day   = std::stoi(date1.substr(8, 2));/* code */
  }
    catch(const std::exception& e)
    {
        if constexpr (debugging >= debugging_ERROR) 
            seqan3::debug_stream << "Error parsing missing in metadata id: " << id 
                                << ", becouse of: " << e.what() << '\n';
    }
 

}

const std::string_view SplitCoVfasta::isolate(const std::string_view virus_name)
{ 
    // extract isolate from virus_name like BTC-4694 from hCoV-19/United Arab Emirates/BTC-4694/2021
    // between the second and third '/'. Return empty string_view if not found.

    std::size_t beg = virus_name.find('/', 8);
    if (beg == std::string_view::npos) return virus_name;  

    std::size_t end = virus_name.find('/', ++beg);
    return virus_name.substr(beg, end - beg); 

    // if the second '/' is found in the last position, beg will = size() and the result is an empty string_view
}

std::ifstream SplitCoVfasta::open_metadata() const
{
    // extract path of metadata file from the fasta file and add name metadata.tsv. Check if it exists.
    std::filesystem::path metadata_fn = dir / "metadata.tsv";
    if (!std::filesystem::exists(metadata_fn))
    {
        seqan3::debug_stream << "Metadata file " << metadata_fn << " does not exist\n";
        return std::ifstream{};
    }
    seqan3::debug_stream << "Reading metadata file " << metadata_fn << '\n';

    // open the metadata file and return the ifstream
    std::ifstream metadata_file{metadata_fn};
    if (!metadata_file.is_open())
    {
        seqan3::debug_stream << "Could not open metadata file " << metadata_fn << '\n';
    }
    return metadata_file;
}

std::unordered_map<std::string, size_t> SplitCoVfasta::parse_metadata_header(std::ifstream& metadata_file)
{   
    // ----------------------------------------------------------------
    // 1. Parse header line to identify column positions
    // ----------------------------------------------------------------
    std::unordered_map<std::string, size_t> col_index; 
    std::string header_line;
    if (!std::getline(metadata_file, header_line))
    {
        seqan3::debug_stream << "The metadata file is empty.\n";
        return col_index;
    }
    if constexpr (debugging)
        seqan3::debug_stream << "Header: " << header_line << '\n';

    // Tokenize the header by tab. We'll store the index in a map keyed by column name.
    // e.g., col might be "Virus name" or "Collection date", etc.

    size_t current_idx = 0;
    std::istringstream iss{header_line};
    std::string col;
    
    while (std::getline(iss, col, '\t'))    col_index[col] = current_idx++;

    return col_index;
}

// read metadata
// GISAID seems to be changing the publication format of the main source of sequences.
// Aparentely, there are no big alignment files, only fasta files and metadata files.
// We will need to parse both and keep the info in sync.
void SplitCoVfasta::read_metadata()
{
    // Virus name	Last vaccinated	Passage details/history	Type	Accession ID	Collection date	Location	Additional location information	Sequence length	Host	Patient age	Gender	Clade	Pango lineage	Pango version	Variant	AA Substitutions	Submission date	Is reference?	Is complete?	Is high coverage?	Is low coverage?	N-Content	GC-Content
    // hCoV-19/United Arab Emirates/BTC-4694/2021		Original	betacoronavirus	EPI_ISL_5142592	2021-01-31	Asia / United Arab Emirates / Abu Dhabi		29884	Human	51	Male	GRY	B.1.1.7	PANGO-v1.23	Former VOC Alpha GRY (B.1.1.7+Q.*) first detected in the UK	(Spike_H69del,NS3_L15F,NS8_Q27stop,NSP3_T183I,Spike_T716I,NSP6_S106del,N_R203K,Spike_A570D,NSP13_K460R,NSP4_F17L,Spike_N501Y,NSP3_I1412T,NS8_R52I,Spike_P681H,Spike_Y144del,NSP6_G107del,NSP3_A890D,Spike_D1118H,NSP6_F108del,NS8_Y73C,N_G204R,Spike_V70del,NSP12_P323L,Spike_D614G,N_D3L,Spike_S982A,N_S235F)	2021-10-14		True	True			0.379400348012
    // hCoV-19/England/QEUH-2C8C50B/2021		Original	betacoronavirus	EPI_ISL_7299806	2021-11-30	Europe / United Kingdom / England		29764	Human	unknown	unknown	GK	AY.4	consensus call	Former VOC Delta GK (B.1.617.2+AY.*) first detected in India	(Spike_T29S,N_G215C,NSP3_A1711V,Spike_T95I,N_D63G,N_R203M,NSP12_G671S,Spike_G142D,NS3_S26L,Spike_P681R,Spike_R158del,NSP3_P1228L,Spike_Y170H,NS7a_V82A,NSP3_A488S,Spike_F157del,NSP4_T492I,NSP14_A394V,Spike_T19R,NS7a_T120I,N_D377Y,M_I82T,Spike_D950N,NSP15_T114A,NS7b_T40I,NSP2_S369F,NSP13_P77L,Spike_E156G,NSP3_P1469S,NSP5_R76K,Spike_T478K,NSP6_T77A,NSP12_P323L,Spike_D614G,Spike_L452R,NSP4_V167L)	2021-12-07		True	True			0.379430836945
    // hCoV-19/USA/MI-CDC-STM-000048060/2021		Original	betacoronavirus	EPI_ISL_1691108	2021-03-31	North America / USA / Michigan	

    auto metadata_file = open_metadata();
    if (!metadata_file.is_open()) return;

    auto col_index = parse_metadata_header(metadata_file);
    if (col_index.empty()) return;

    // For convenience, find column indexes (if they exist) once:
    // We'll do lookups like col_index.find("Virus name"). 
    // (some might not exist in older metadata versions)
    auto idx_virusname_it  = col_index.find("Virus name");       // hCoV-19/United Arab Emirates/BTC-4694/2021	
    auto idx_EPI_ISL_it    = col_index.find("Accession ID");     // EPI_ISL_5142592
    auto idx_date_it       = col_index.find("Collection date");  // 2021-01-31
    auto idx_location_it   = col_index.find("Location");         // Asia / United Arab Emirates / Abu Dhabi
    auto idx_clade_it      = col_index.find("Clade");            // GRY
    auto idx_pango_it      = col_index.find("Pango lineage");    // B.1.1.7
    auto idx_pangover_it   = col_index.find("Pango version");    // PANGO-v1.23

    // Retrieve needed columns, if present: If not found, skip that field.

    // ----------------------------------------------------------------
    // 2. Prepare to parse lines in a memory-friendly manner
    // ----------------------------------------------------------------
    // We will parse them in-place, storing string_views that point into `line`.
    std::string line;
    line.reserve(4096);                     // a guess; adjust if lines are large
    std::vector<std::string_view> cols;
    cols.reserve(64);                       // typical # of columns (over-reserve to avoid repeated allocations)

    // For 17M lines, we might want to disable sync with stdio for speed: (some do this at the start of main)
    // std::ios::sync_with_stdio(false);  std::cin.tie(nullptr);

    long date_errors = 0;
    size_t lines_parsed = 0;                // todo print progress every 100k? lines or so
    size_t duplicates = 0;                  // count duplicates
    size_t empty_virusname = 0;             // count empty virus names
    
    long m{(1L<<19)-1};                            // for progress printing
    seqan3::debug_stream << "\nchunk - m= " << m << "\n" ; 

    while (std::getline(metadata_file, line))
    {
        if (line.empty()) continue;
        if (!(++lines_parsed & m))                         // progress indication every 2^19 (m) sequences 
        {
            if constexpr (debugging >= debugging_INFO) seqan3::debug_stream << "\tlines_parsed= " << lines_parsed  << "\n" ;
        }
        //if constexpr (debugging >= break_by) if (lines_parsed >= break_point) break;     // set for debugging only  !!!!!!!!!!!! 
        
        // 2.1. Tokenize this line by tab into `cols` as string_views
        cols.clear();
        {
            // We create views into `line` in-place
            char * str_data = line.data();
            size_t len      = line.size();
            size_t start    = 0;

            for (size_t i = 0; i < len; ++i)
            {
                if (line[i] == '\t')
                {
                    cols.emplace_back(str_data + start, i - start);
                    start = i + 1;
                }
            }
            if (start < len)                // last column after final tab (or if no tabs at all)
                cols.emplace_back(str_data + start, len - start);
        }
        // We now have a vector of string_view columns. Let's extract the relevant fields.

                               // skip if no column or no data - we need at least virus name, expected first field.
        if (idx_virusname_it->second >= cols.size())  { empty_virusname++; continue; }  // skip if no virus name
        std::string_view virus_name_sv = cols[idx_virusname_it->second];
        if (virus_name_sv.empty())                    { empty_virusname++; continue; }  // skip if no virus name

        parsed_id& pid = metadata[std::string{virus_name_sv}];          // reference to the newly inserted parsed_id
        // Create the parsed_id in-place in the map, keyed by virus_name
        // We can do an emplace -> returns pair<iterator,bool>
        // We'll fill the struct with brace initialization.
        // auto [it, inserted] = metadata.try_emplace(std::string{virus_name_sv}, parsed_id{});
        // if already in map, skip or overwrite? 

        if (!pid.isolate.empty()) {duplicates++; continue;}             // Let’s assume we skip re-parsing if it’s already inserted:

        pid.isolate = isolate(virus_name_sv);                           // parse the isolate from the "virus name" 

        if (idx_date_it->second < cols.size())                          // 2.2. Parse Collection date -> pid.year, pid.month, pid.day. Working OK
        {
            std::string_view date_sv = cols[idx_date_it->second];
           
            if (date_sv.size() == 10)                                   // e.g. "2021-03-31"
            {
                auto d = date_sv.data();
              
                auto fc1 = std::from_chars(d    , d +  4, pid.year );  // parse year
                auto fc2 = std::from_chars(d + 5, d +  7, pid.month);  // parse month
                auto fc3 = std::from_chars(d + 8, d + 10, pid.day  );  // parse day

                          // optional: check fc1.ec, fc2.ec, fc3.ec for parse errors
            }
            else          // log debug, and try using chrono to parse the date
            if (date_sv.size() == 7)                                   // e.g. "2021-03"
            {
                auto d = date_sv.data();
              
                auto fc1 = std::from_chars(d    , d +  4, pid.year );   
                auto fc2 = std::from_chars(d + 5, d +  7, pid.month);   
             }
            else         
            if (date_sv.size() == 4)                                   // e.g. "2021"
            {
                auto d = date_sv.data();
              
                auto fc1 = std::from_chars(d    , d +  4, pid.year );   
            }
            else  
            if (date_sv.size() == 6)                                   // e.g. "2021-3"
            {
                auto d = date_sv.data();
              
                auto fc1 = std::from_chars(d    , d +  4, pid.year ); 
                auto fc2 = std::from_chars(d + 5, d +  6, pid.month);   
            }
            else
            if (date_sv.size() == 9)                              // e.g. "2023-12-4" or "2023-9-25"
            {
                auto d = date_sv.data();
              
                auto fc1 = std::from_chars(d    , d +  4, pid.year );   
                if (d[7] == '-')  // 2023-9-25
                {
                    auto fc2 = std::from_chars(d + 5, d +  6, pid.month);   
                    auto fc3 = std::from_chars(d + 8, d +  9, pid.day  );  
                }
                else  // 2023-12-4
                {
                    auto fc2 = std::from_chars(d + 5, d +  7, pid.month);  
                    auto fc3 = std::from_chars(d + 8, d +  9, pid.day  );  
                }
            }
            // 2024-4-2
            else if (date_sv.size() == 8)                              // e.g. "2024-4-2"
            {
                auto d = date_sv.data();
              
                auto fc1 = std::from_chars(d    , d +  4, pid.year );   
                auto fc2 = std::from_chars(d + 5, d +  6, pid.month);   
                auto fc3 = std::from_chars(d + 7, d +  8, pid.day  );  
            }
            else                              // log debug, and try using chrono to parse the date?
            {
                if (date_errors++ < 50)
                    seqan3::debug_stream << "ERROR: [read_metadata] " << date_errors << "- Unexpected date format: " << date_sv << " in " << line << '\n';

                // std::chrono::year_month_day ymd;
                // // try to parse the date using chrono    from_stream 
                // try
                // {
                //     std::istringstream ss{std::string{date_sv}};
                //     std::chrono::from_stream(ss, "%F", ymd);  // no from_stream available in gcc C++20 ??
                //     pid.year  = static_cast<int     >(ymd.year ());
                //     pid.month = static_cast<unsigned>(ymd.month());
                //     pid.day   = static_cast<unsigned>(ymd.day  ());
                //     seqan3::debug_stream << "[read_metadata] Parsed date using chrono: " << pid.year << '-' << pid.month << '-' << pid.day << '\n';
                // }
                // catch (const std::exception& e)
                // {
                //     seqan3::debug_stream << "[read_metadata] Failed to parse date using chrono: " << e.what() << '\n';
                // }
            }
        }

        if (idx_location_it->second < cols.size())        // 2.3. Parse Location -> pid.continent, pid.country, pid.region
        {
            std::string_view loc_sv = cols[idx_location_it->second];

            //     e.g. "North America / USA / Michigan",  We'll do a split by '/' to get up to 3 parts

            auto continent_end = loc_sv.find('/');
            if (continent_end == std::string_view::npos)   
                pid.continent = loc_sv;                           // log an error? just set continent to the whole string
            else
            {
                pid.continent = loc_sv.substr(0, continent_end - 1);
                auto country_start = continent_end + 2;
                if (country_start < loc_sv.size())                    // there is a country
                {  
                    auto country_end   = loc_sv.find('/', country_start);
                    if (country_end == std::string_view::npos)        // log an error? just set country to the whole string
                        pid.country = loc_sv.substr(country_start);
                    else
                    {
                        pid.country = loc_sv.substr(country_start, country_end - 1 - country_start);
                        auto region_start = country_end + 2;
                        if (region_start < loc_sv.size())  // there is a region
                        {  
                            pid.region  = loc_sv.substr(region_start);
                        }
                    }
                }
            }
        }

        if ( idx_clade_it->second < cols.size())                        // 2.4. Clade
             pid.clade = cols[idx_clade_it->second];  

        if (idx_pango_it->second < cols.size())                         // 2.5. Pango lineage
            pid.pango = cols[idx_pango_it->second];

        if (idx_pangover_it->second < cols.size())                      // 2.6. Pango version
            pid.pango_version  = cols[idx_pangover_it->second];

        if (idx_EPI_ISL_it->second < cols.size())                       // 2.7. EPI_ISL
            pid.EPI_ISL = cols[idx_EPI_ISL_it->second];
    }

    if constexpr (debugging >= debugging_INFO)
        seqan3::debug_stream << "[read_metadata] Parsed "  << lines_parsed 
                             << " lines. metadata.size()=" << metadata.size() 
                             << ", duplicates=" << duplicates
                             << ", empty_virusname=" << empty_virusname
                             << ", date_errors=" << date_errors << '\n';
}

void parsed_id::parse_id_allnuc(const std::string& id) 
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
    std::size_t country_end = id.find('/',   country_beg) - 1;          
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

    isolate = id.substr(isolate_beg, isolate_len);                 
    EPI_ISL = id.substr(EPI_ISL_beg, EPI_ISL_len); 
    country = id.substr(region_beg) + " - " + id.substr(country_beg, country_len) ;  

    auto s = id.data() ;
    auto y = s + year_beg ;
    auto m = s + month_beg;
    auto d = s + day_beg  ;    
    
    std::from_chars(y, m    , year );   
    std::from_chars(m, d    , month); 
    std::from_chars(d, d + 3, day  );

     seqan3::debug_stream << "\n" << isolate << " - EPI: " << EPI_ISL << " - " 
                          << year << " - "   << month << " - " << day 
                         << " - " << country << '\n' ; 
}

// check format by reading and checking the first record. todo: too simple, make more checks??
GISAID_format SplitCoVfasta::check_format()
{
    // Initialise a file input object with the fasta file.
    seqan3::sequence_file_input<MSA> file_in{fasta};

    for (auto&& ref_rec : file_in)                               // read the first record ONLY
    {    
        if (ref_rec.id().find("EPI_ISL_") != std::string::npos)  // we assume it is msa format
            format = GISAID_format::msa;
        else
            format = GISAID_format::fasta;

        // if debuging print the format found
        if constexpr (debugging) 
            seqan3::debug_stream << "\nFormat: " << (format == GISAID_format::msa ? "msa" : "fasta") << '\n';   

        return format; 
    }
    throw std::runtime_error{"ERROR: Empty fasta file"};    
}

void SplitCoVfasta::split( )
{
    // launched inmediatelly after the SplitCoVfasta constructor, and teh whole SplitCoVfasta is destructed after this call.
    // No data race posible
    // use chrono to logg the time reading metadata and splitting
    auto start = std::chrono::high_resolution_clock::now();
    try
    {
        check_format();
        read_metadata();      // fill metadata map keyed by Virus name
    }
    catch(const std::exception& e)
    {
        throw std::runtime_error{ std::string("ERROR setting format or reading metadata") + e.what()};
    }
    catch(...)
    {
        throw std::runtime_error{ "ERROR setting format or reading metadata"};
    }
    
    auto stop_metadata = std::chrono::high_resolution_clock::now(); 

    try
    {
        if (format == GISAID_format::fasta)    split_fasta();
        else                                     split_msa();        
    }
    catch(const std::exception& e)
    {
        throw std::runtime_error{ std::string("ERROR splitting the fasta file - ") + e.what()};
    }
    catch(...)
    {
        throw std::runtime_error{ "ERROR splitting the fasta file"};
    }

    // take time and logg in hours, minutes, seconds
    auto stop_split = std::chrono::high_resolution_clock::now();
    
    // logg in hours, minutes, seconds: hh:mm:ss format
    auto duration_metadata = std::chrono::duration_cast<std::chrono::seconds>(stop_metadata - start);
    auto duration_split    = std::chrono::duration_cast<std::chrono::seconds>(stop_split    - stop_metadata);
    if constexpr (debugging >= debugging_INFO)
        seqan3::debug_stream << "Reading metadata: " << duration_metadata.count() / 3600 << ':' 
                                                     << (duration_metadata.count() % 3600) / 60 << ':' 
                                                     << duration_metadata.count() % 60 << '\n';
    if constexpr (debugging >= debugging_INFO)
        seqan3::debug_stream << "Split all sequences: " << duration_split.count() / 3600 << ':' 
                                                        << (duration_split.count() % 3600) / 60 << ':' 
                                                         << duration_split.count() % 60 << '\n';
}

void skip_bad(auto& it, const auto& end)
{
    bool bad = true;
    while (bad && it != end)
    {
        try
        {
            ++it;
            bad = false;
        }
        catch(const std::exception& e)
        {
            if constexpr (debugging >= debugging_TRACE) seqan3::debug_stream << "ERROR reading sequence: " << e.what() << '\n';

            try {it.seek_to(it.file_position()+1);} 
            catch(...) {}  
        }
    }
}

void SplitCoVfasta::split_fasta( )
{
    set_ref_seq();                                        // load the reference sequence
    
    seqan3::sequence_file_input<OLIGO> file_in{fasta};    // Initialise a file input object with a FASTA file.

    long t{0L}, m{(1L<<16)-1};                            // for progress printing
    seqan3::debug_stream << "\nchunk - m= " << m << "\n" ; 

    std::string id ;

    decltype(file_in.begin()) it;
    try { it = file_in.begin(); } catch(...) {skip_bad(it, file_in.end());}            // some seq. are proteins not nucleotides

    while(it != file_in.end())  // read each sequence in the file sequencially !! One by one record!
    {    
            auto && record = *it;
            id = record.id();    

            if constexpr (debugging >= debugging_TRACE)  seqan3::debug_stream << "\n" << record.id() << '\n' ;
            
            std::string_view virus_name = record.id();
            virus_name = virus_name.substr(0, virus_name.find('|'));          // if not found, keep the whole virus_name

            if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "Virus name: " << virus_name << '\n';

            parsed_id& pid = metadata[std::string{virus_name}]; // reference to the (posibly newly inserted) parsed_id, always stable living in metadata map

            if (pid.isolate.empty())                                  // if isolate was not set from metadata file, parse the id
            {
                if constexpr (debugging >= debugging_WARNING) seqan3::debug_stream << "Not found in metadata. Going Parsing id: " << record.id() << '\n';
                pid.parse_id(record.id());
            }
            if constexpr (debugging >= debugging_TRACE+3)        seqan3::debug_stream << "Parsed id: " << pid.isolate << '\n';

            //for (auto & pcr : pcrs)  
            std::for_each(std::execution::par_unseq, pcrs.begin(), pcrs.end(), [&](PCRSplitter& pcr) 
            {
                // each in paralell PCR check recive a const reference to the record, 
                // and modify only that PCR specific data. No data race posible.
                // Every PCR keep a link to the parent SplitCoVfasta object, but it is const, read only. 

                if constexpr (debugging >= debugging_TRACE+3) seqan3::debug_stream << "Trace before check_rec: " ;

                auto tcop = pcr.check_rec(record);  // check in parallel for each PCR/target - that tc belong to current PCR - no data race
                if (tcop) 
                {
                    target_count& tc = *tcop;
                    if constexpr (debugging >= debugging_TRACE)   seqan3::debug_stream << "Trace after check_rec returning target_pattern: " << tc.target.target_pattern <<"\n" ;
                
                    tc.update_counts_from(pid);
                }
                else 
                {if constexpr (debugging >= debugging_TRACE)   seqan3::debug_stream << "Trace after check_rec FAILED\n" ;;}
            }
            );  

            try
            {
                ++it;
            }
            catch(const std::exception& e)
            {
                if constexpr (debugging >= debugging_ERROR) seqan3::debug_stream << "ERROR reading the sequence after: " << t << ": " << id << " becouse of: " << e.what() << '\n';
                skip_bad(it, file_in.end());
            }

        if (!(++t & m))                                              // progress indication every 2^18 (m) sequences 
        {
            seqan3::debug_stream << "\tT= " << t  << "\n" ;
            for (auto & pcr : pcrs)
                seqan3::debug_stream << pcr.pcr_name <<"= " << pcr.count 
                                     << ". Grouped: "    << pcr.grouped.size() << "\n" ; 
        }
        if constexpr (debugging >= break_by) if (t>=break_point) break;     // set for debugging only  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    
    seqan3::debug_stream << "\nTotal= " << t  << "\n" ; 

    for (auto & pcr : pcrs)                // write grouped sequences for each PCR/target
    {
        seqan3::debug_stream << pcr.pcr_name <<"= " << pcr.count 
                             << ". Grouped: "    << pcr.grouped.size() << ". " ; 
        pcr.write_grouped();
        // by debugging level print the amplicon_pos_beg (map of positiion: count) in descending order of counts (values)
        if constexpr (debugging >= debugging_NOT_USED)
        {
            //std::multimap<int, int, std::greater<>> inverted{pcr.amplicon_pos_beg.begin(), pcr.amplicon_pos_beg.end()};
            // print time durations in seconds
            seqan3::debug_stream << "Amplicon_reference pos beg: " << pcr.ref_beg << '\n'; 
            for (auto& [pos, count] : pcr.amplicon_pos_beg)
                seqan3::debug_stream << "pos: "  << pos << ", count: " << count << '\n';
            seqan3::debug_stream << "Detected with Hint1: " << pcr.count_hint1 << " used " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint1_time_used).count() << " s, wasted " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint1_time_wasted).count() << " s\n";
            seqan3::debug_stream << "Detected with Hint2: " << pcr.count_hint2 << " used " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint2_time_used).count() << " s, wasted " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint2_time_wasted).count() << " s\n";
            seqan3::debug_stream << "Detected with Hint3: " << pcr.count_hint3 << " used " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint3_time_used).count() << " s, wasted " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint3_time_wasted).count() << " s\n";
            seqan3::debug_stream << "Detected with Hint4: " << pcr.count_hint4 << " used " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint4_time_used).count() << " s, wasted " << std::chrono::duration_cast<std::chrono::seconds>(pcr.hint4_time_wasted).count() << " s\n";
            seqan3::debug_stream << "Detected with full seq: "<<pcr.count_full << " used " << std::chrono::duration_cast<std::chrono::seconds>(pcr.full_seq_time_used).count() << " s, wasted " << std::chrono::duration_cast<std::chrono::seconds>(pcr.full_seq_time_wasted).count() << " s\n";
            seqan3::debug_stream << "Not detected: "         << pcr.count_not_found << ", with lenths: \n";
            // print the not_found_lenghts
            for (auto& [len, count] : pcr.not_found_lenghts)
                seqan3::debug_stream << "len: " << len << ", count: " << count << '\n';
            }
    }

}

void SplitCoVfasta::split_msa( )
{
    set_msa_ref_pos();    // load the reference, set MSA positions, etc.
 
    // Initialise a file input object with a FASTA file.
    seqan3::sequence_file_input<MSA> file_in{fasta};

    long t{0L}, m{(1L<<18)-1};  // for progress printing
    seqan3::debug_stream << "\nchunk - m= " << m << "\n" ; 


    for (auto && record : file_in)             // read each sequence in the file
    {
        if constexpr (debugging)
            seqan3::debug_stream << "\n" << record.id() << '\n' ;

        std::string_view virus_name = record.id();
        virus_name = virus_name.substr(0, virus_name.find('|'));

        parsed_id& pid = metadata[std::string{virus_name}]; // reference to the newly inserted parsed_id
        if (pid.isolate.empty())  // if isolate is not set, parse the id
        {
            if constexpr (debugging)
                seqan3::debug_stream << "Not found in metadata. Going Parsing id: " << record.id() << '\n';
            pid.parse_id_allnuc(record.id());
        }

        //for (auto & pcr : pcrs)   
        std::for_each(std::execution::par_unseq, pcrs.begin(), pcrs.end(), [&](auto& pcr) 
        {
            target_count& tc = pcr.check_msa_rec(record);
            tc.update_counts_from(pid);
        });

        if (!(++t & m))                      // print a dot every 2^18 sequences for progress indication
        {
            // seqan3::debug_stream << '.' ;
        
            seqan3::debug_stream << "\tT= " << t  << "\n" ;
            for (auto & pcr : pcrs)
                seqan3::debug_stream << pcr.pcr_name <<"= " << pcr.count 
                                     << ". Grouped: "    << pcr.msa_grouped.size() << "\n" ; 
        }
        if (t>17000) break;
    }
    seqan3::debug_stream << "\nTotal= " << t  << "\n" ; 

    for (auto & pcr : pcrs)                // write grouped sequences for each PCR/target
    {
        seqan3::debug_stream << pcr.pcr_name <<"= " << pcr.count 
                             << ". Grouped: "    << pcr.msa_grouped.size() << ". " ; 
        pcr.write_msa_grouped();
        // by debugging level print the amplicon_pos_beg (map of positiion: count) in descending order of counts (values)
        if constexpr (debugging >= debugging_INFO)
        {
            std::multimap<int, int, std::greater<>> inverted{pcr.amplicon_pos_beg.begin(), pcr.amplicon_pos_beg.end()};
            for (auto& [pos, count] : inverted)
                seqan3::debug_stream << "pos: " << pos << ", count: " << count << '\n';
        }


    }

}

void target_count::update_counts_from(parsed_id& pid)
{
    if (!count) return;

    year_count& yc = years[pid.year];
    yc.count++;
    
    month_count& mc = yc.months[pid.month];
    mc.count++;
    
    day_count& dc = mc.days[pid.day];
    dc.count++;
                
    continent_count& ct = dc.continents[pid.continent];
    ct.count++;

    country_count& cc = ct.countries[pid.country];
    cc.count++;

    clade_count& cd = cc.clades[pid.clade];
    cd.count++;            
    if (!cd.id) cd.id = &pid;
}

}  // namespace cov

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
}
*/
