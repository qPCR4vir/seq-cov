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

#include "seq-cov.hpp"

/*
struct dna_deg : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna15; // instead of dna5
    template <typename alph>
    using sequence_container = seqan3::dna15_vector;
};
*/


// struct SplitGene
 
SplitCoVfasta::SplitGene::SplitGene(SplitCoVfasta const &parent, std::string gene, bool split, bool group)
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

bool SplitCoVfasta::SplitGene::read_oligos(const std::filesystem::path& path_oligos)
{
        if (path_oligos.empty()) return false;   // todo: more checks?

        seqan3::debug_stream << "\nGoing to load: " << path_oligos.string();
        seqan3::sequence_file_input file_in{path_oligos};
        std::string fw, rv;
        beg = end = 0;
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
            pr.beg = std::stoi(seq_pos.substr(0, dash_pos));
            pr.end = std::stoi(seq_pos.substr(dash_pos + 1));
            // todo check if beg, end are valid
            seqan3::debug_stream << " with beg: " << pr.beg << " and end: " << pr.end;
            pr.seq  = primer.sequence();
            pr.match = std::round(parent.match * double(pr.seq.size()) / 100.0);
            seqan3::debug_stream << " with minimum matches:" << pr.match;

            if (pr.beg < pr.end)  // one forward primer/prbe
            {
                if ( !beg || beg > pr.beg)  // extern forward primer
                {
                    beg = pr.beg;
                    forw_idx = f_primers.size();
                }   
                f_primers.push_back(pr);                     
                continue;
            }
            else  // one reverse primer/prbe
            {
                if ( !end || end < pr.end)  // extern reverse primer
                {
                    end = pr.end;
                    rev_idx = r_primers.size();
                }
                r_primers.push_back(pr);
                continue;
            }
        }
        return true;
}

bool check(auto& record)  /// record identified and ...?
{
    // seqan3::debug_stream << '\n' << record.id();
    if (ignore) return false;
    if (!full_msa)
    {
        if (!record.id().starts_with(start)) return false;

        if (split) file_fasta_split.push_back(record);  // todo deprecate
    }

    if (group) check_rec(record);
    
    return true;
}

bool SplitCoVfasta::SplitGene::check_rec(auto& record)
{
    sequence_type &sq = record.sequence();
    count++;
    int flank = parent.flank;
    // >Gene name|Isolate name|YYYY-MM-DD|Isolate ID|Passage details/history|Type^^
    //           1            212345678910
    std::string id = std::move(record.id());  // each record is checked only once

    std::size_t date_beg = id.find('|', 1 + id.find('|')) + 1;           // pos next to 2nd |
    std::string date     = id.substr(date_beg, id.find('|', date_beg) - date_beg);  // 9 or 10
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

    SeqGr& sg1 = grouped[cur_seq];                // first try
    if (sg1.count++)                              // count not 0 ==> duplicate seq. More than 99% of cases.
    {
        SeqGr& sgd = daily   [std::move(date)][cur_seq];            
        if (!sgd.count++)                         // but it is new in daily, count 0 ==> new seq for this day, blank grd 
        {
            sgd.id = id;
            sgd.beg = sg1.beg ;
            sgd.end = sg1.end ;
        }
        SeqGr& sgm = monthly[std::move(month)][std::move(cur_seq)];  // we are about to return - avoid copy          
        if (!sgm.count++)                                            // but it is new in monthly
        {
            sgm.id = std::move(id);
            sgm.beg = sg1.beg ;
            sgm.end = sg1.end ;
        }
        return true;
    }

    // new, unknown seq if limited as the first seq. sg1 is a new, blank group that is already in grouped

    SeqGr sg=set_seq_pos(sq);   // first try to find a true align to correct position - find limits
    
    if (sg.beg == notfound || sg.end == notfound)
    {
        // failed!!! NNN??, is new in grouped, but we report it for completeneds of staistic ??!
        sg1.id = std::move(id);
        sg1.beg = sg.beg ;
        sg1.end = sg.end ;                               
        
        // sg1.count already set to 1 becouse ir is new in global grouped ==> new in daily and monthly too

        daily   [std::move(date)][          cur_seq ] = sg1;            
        monthly[std::move(month)][std::move(cur_seq)] = sg1;              

        //seqan3::debug_stream  << start <<  seqan3::dna5_vector{bg, en} << " -- Failed! " 
        //                      << " (" << lbeg << ", " << lend  <<")\n" ;
        return true; // todo ???????????????
    }
    //seqan3::debug_stream << start <<  seqan3::dna5_vector{bg, en} 
    //                     << " (" << lbeg << ", " << lend  <<")\n" ;

    int nlend= std::min<int>(sg.end + flank, sq.size() );
    int nlbeg= std::max<int>(sg.beg - flank, 0);
    
    if (lbeg == nlbeg && lend == nlbeg)  // new seq in the same pos, reuse the positions  ????????
    {
        sg1.id = id;
        sg1.beg = sg.beg - nlbeg;
        sg1.end = sg.end - nlbeg;

        daily   [std::move(date)][          cur_seq ] = sg1;            
        monthly[std::move(month)][std::move(cur_seq)] = sg1;   

        return true;
    }
        
    // new seq in a new pos. Lets make a position correction and try againg. 

    grouped.erase(cur_seq);
    //seqan3::debug_stream << " New try ----- (" << lbeg << ", " << lend  <<")\n" ;

    auto new_seq = sequence_type{sq.begin()+nlbeg, sq.begin()+nlend};
    SeqGr&   sgr = grouped[new_seq];

    if (sgr.count++)                                // duplicate seq. More than 99% of cases.
    {
        SeqGr& sgd = daily   [std::move(date)][new_seq];            
        if (!sgd.count++)                           // but it is new in daily
        {
            sgd.id = id;
            sgd.beg = sgr.beg ;
            sgd.end = sgr.end ;
        }
        SeqGr& sgm = monthly[std::move(month)][std::move(new_seq)];            
        if (!sgm.count++)  // but it is new in monthly
        {
            sgm.id = std::move(id);
            sgm.beg = sgr.beg ;
            sgm.end = sgr.end ;
        }
        return true;
    }
    
    //seqan3::debug_stream << " It was new !!!!!!!\n" ;
    sgr.beg = sg.beg - nlbeg;
    sgr.end = sg.end - nlbeg;
    sgr.id = std::move(id);

    daily   [std::move(date)][          new_seq ] = sgr;            
    monthly[std::move(month)][std::move(new_seq)] = sgr;   
    
    return true;
}

SeqGr SplitCoVfasta::SplitGene::set_seq_pos(const sequence_type& s)
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
    // todo  align with ref_seq

    return sg;
}

void SplitCoVfasta::SplitGene::write_grouped ()
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
            auto id = "m_" + date + 
                    +"-15_x_"+ std::to_string(gr.count) 
                    + "_" + std::to_string(gr.beg) 
                    + "_" + std::to_string(gr.end)
                    + "_" + gr.id ;
            file_e_grm.push_back(record_t{std::move(seq), std::move(id)});
        }

}
        
// class SplitCoVfasta
 
void SplitCoVfasta::split_fasta( )
{
    
    // Initialise a file input object with a FASTA file.
    seqan3::sequence_file_input file_in{fasta};

    long t{0L}, m{(1L<<18)-1};
    seqan3::debug_stream << "\nchunk - m= " << m << "\n" ; 

    for (auto & record : file_in)             // read each sequence in the file
    {
        for (auto & gene : genes)             // for each sequence, check each gene/target
            if (gene.check(record)) 
                break;

        if (!(++t & m))                      // print a dot every 2^18 sequences for progress indication
            seqan3::debug_stream << '.' ;
        {
            seqan3::debug_stream << "\tT= " << t  << "\n" ;
            for (auto & gene : genes)
                seqan3::debug_stream << gene.gene <<"= " << gene.count 
                                     << ". Grouped: "    << gene.grouped.size() << "\n" ; 
        }
        //if (t>10000000) break;
    }
    seqan3::debug_stream << "\nTotal= " << t  << "\n" ; 

    for (auto & gene : genes)                // write grouped sequences for each gene/target
    {
        seqan3::debug_stream << gene.gene <<"= " << gene.count 
                             << ". Grouped: "    << gene.grouped.size() << ". " ; 
        gene.write_grouped();
    }

}
