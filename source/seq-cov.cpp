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
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/sequence_file/all.hpp>      // for sequence_file_input and sequence_file_output
#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>

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
 
SplitGene::SplitGene(SplitCoVfasta const &parent, std::string gene)
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
        beg = end = 0;
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
                all_oligos.push_back(pr);                  
                continue;
            }
            else  // one reverse primer/prbe
            {
                if ( !end || end < pr.beg)  // extern reverse primer
                {
                    end = pr.beg;
                    rev_idx = r_primers.size();
                }
                r_primers.push_back(pr);
                all_oligos.push_back(pr);    
                continue;
            }
        }
        /*for (auto & primer : f_primers) all_oligos.push_back(&primer);
        for (auto & primer : r_primers) all_oligos.push_back(&primer);
        for (auto & primer : probes_s) all_oligos.push_back(&primer);
        for (auto & primer : probes_a) all_oligos.push_back(&primer);*/

        return true;
}

bool SplitGene::reconstruct_seq(const msa_seq_t& s, oligo_seq_t& seq, 
                                                int& beg, int& end, std::vector<int>& msa_pos, int tent_len)
{
    seq.clear();
    seq.reserve(tent_len);
    msa_pos.clear();
    msa_pos.reserve(tent_len);
    bool reverse = (beg > end);
    // go through the sequence and eliminate gaps to reconstruct the original sequence
    if (reverse) 
    {
        for ( int i = parent.msa_pos[beg-1]; i >= parent.msa_pos[end-1]; --i)
        {
            if (s[i].holds_alternative<seqan3::gap>() ) continue;
            seq.push_back(s[i].convert_unsafely_to<oligo_seq_alph>().complement());
            msa_pos.push_back(i);
        }
        return true;
    }
    for ( int i = parent.msa_pos[beg-1]; i <= parent.msa_pos[end-1]; ++i)
    {
        if (s[i].holds_alternative<seqan3::gap>() ) continue;
        seq.push_back(s[i].convert_unsafely_to<oligo_seq_alph>());
        msa_pos.push_back(i);
    }
    return true;
}

void SplitGene::align(pattern_q & pq, msa_seq_t& target)
{
    static long count{0L};
    count++;
    seqan3::debug_stream << "\n " << count << ": To do - align Primer: " << pq.primer.name << " \n" ;
                        //<< target <<'\n' ;
                        //<< pattern << '\n'  << target << '\n' ;  // << "Misatches: " << mm << ", Ns: " << N << ", crit: " << crit << '\n';
}
void SplitGene::target_pattern(target_q & tq, msa_seq_t& sq, long msa_beg, long msa_end)
{
    tq.target_pattern.clear();
    tq.target_pattern.reserve(end-beg);

    for ( int i = msa_beg; i <= msa_end; ++i)
    {
        if (parent.msa_ref[i].holds_alternative<seqan3::gap>() &&
                        sq[i].holds_alternative<seqan3::gap>()    )           continue;

        if (parent.msa_ref[i].holds_alternative<seqan3::gap>() ||
                        sq[i].holds_alternative<seqan3::gap>()    )          
        {    
            tq.target_pattern.push_back(sq[i].to_char());     
            continue;  
        }      
        auto s =             sq[i].convert_unsafely_to<oligo_seq_alph>();
        auto r = parent.msa_ref[i].convert_unsafely_to<oligo_seq_alph>();
        if (!mismatch.score(s, r))   
        { 
            tq.target_pattern.push_back('.');  
            continue;
        }
        tq.target_pattern.push_back(s.to_char());   
    }
}

void SplitGene::evaluate_target(target_q & tq, msa_seq_t& sq, long msa_beg, long msa_end)
{
    for (auto & primer : all_oligos)  // todo: parallelize (std::execution::par)
    {
        tq.patterns.emplace_back(primer);
    }
    std::for_each(std::execution::par, tq.patterns.begin(), tq.patterns.end(), [&](pattern_q &pq)
    {
        evaluate_target_primer(pq, sq);
    });
    target_pattern(tq, sq, msa_beg, msa_end);

    /*for (auto & primer : f_primers)
    {
        evaluate_target_primer(tq, primer, sq);
    }
    for (auto & primer : r_primers)
    {
        evaluate_target_primer(tq, primer, sq);
    }
    for (auto & probe : probes_s)
    {
        evaluate_target_primer(tq, probe, sq);
    }
    for (auto & probe : probes_a)
    {
        evaluate_target_primer(tq, probe, sq);
    }
    /* for (auto & pq : tq.patterns)
        seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             //<< target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 */}

void SplitGene::evaluate_target_primer(pattern_q &pq, cov::msa_seq_t &sq)
{
    // pattern_q &pq = tq.patterns.emplace_back(primer);
    cov::oligo &primer = pq.primer;
    oligo_seq_t target;
    reconstruct_seq(sq, target, primer.beg, primer.end, msa_target_pos, primer.seq.size());
    if (target.size() != primer.seq.size())
        return align(pq, sq); // todo: try to align the primer with the target sequence
 
    pq.pattern = std::string(primer.seq.size(), '.');
    int len = primer.seq.size();
    for (int i = 0; i < len; ++i)
    {
        if (!mismatch.score(primer.seq[i], target[i]))   continue;

        if  (target[i] == 'N'_dna15)  
        {
            pq.N++;
            pq.pattern[i] = 'N';
            continue;
        }
        pq.pattern[i] = target[i].to_char();
        pq.mm++;
        if (len - i <= parent.crit_term_nt)                   pq.crit++;
     }
    pq.Q = pq.mm + pq.crit * 4;
/*     seqan3::debug_stream << "\nPrimer: " << pq.primer.name << ":\n" 
                             << pq.primer.seq <<'\n' 
                             << pq.pattern << '\n'
                             << target << '\n' 
                             << "Q = " << pq.Q << ", Missatches: " << pq.mm << ", Ns: " << pq.N << ", crit: " << pq.crit << '\n';
 */
}

target_count& SplitGene::check_rec(auto& record)
{
    
    msa_seq_t& sq = record.sequence();
    count++;
    long msa_beg{0}, msa_end{0};
    if ( beg < end) { msa_beg = parent.msa_pos[beg-1]; msa_end = parent.msa_pos[end-1]; }
    else            { msa_beg = parent.msa_pos[end-1]; msa_end = parent.msa_pos[beg-1]; }

    target_count & target_c = grouped[{sq.begin()+msa_beg, sq.begin()+msa_end}];  // todo: check if it is new
    if (!target_c.count)  // new target sequence
        evaluate_target(target_c.target, sq, msa_beg, msa_end);
    
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
    
    // for (auto & gene : genes) gene.set_ref_pos(); \todo: check correct positions of primers on the reference sequence

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


        std::for_each(std::execution::par, genes.begin(), genes.end(), [&](auto& gene) 
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
        });

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