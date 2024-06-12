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

#include "gui.hpp"
#include "seq-cov.hpp"
using namespace cov;
GUI::GUI() : nana::form{nana::api::make_center(1000, 350)}
{
    caption("Split-CoV-fasta. v3.02.01");

    input_file.tip_string("Original fasta file:").multi_lines(false);
    flank.range(0, 100, 1);
    flank.value("5");
    match.range(50.0, 100.0, 1.0);
    match.value("70.0");
    period.range(0, 12, 1);
    period.value("3");

    run_split.events().click([&]()
    {
        std::filesystem::path fasta{input_file.text()};
        if (!std::filesystem::is_regular_file(fasta)) return;  // todo msg

        SplitCoVfasta sp{fasta, flank.to_int(), match.to_double(), from.text(), to.text()};

        auto Add_Gene = [&sp](auto &g)
        {
        if (!g.forw.text().empty() && !g.rev.text().empty())
            sp.add_gene(g.input_file.text(),
                        g.gene.text(),     
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
                <height=25 input arrange=[variable, 35, 50,  75,50,  60, 60] gap=7> 
                <height=10  >
                <height=30 file >
                <height=10  >
                <height=133 vertical genes gap=10> 
                <height=10  >
                <height=30 no_dates arrange=[ 35,80,  35,80,  45,40,180] gap=10> 
                >   
        )");
    p["input"] << input_file_label  << "Flank:"  << flank 
                                   << "Min match %" << match     << set  << run_split ;
    p["file"]  << input_file ;
    p["genes"] << E << N << S << R;
    p["dates"] << "From date: " << from << " To date: " << to
               << "Separate by: " << period << " months \n(0 - all time in 1 period).";
    p.collocate();
};

GeneGUI::GeneGUI(nana::window parent, std::string gene_name, std::string fw, std::string rv)
    : nana::panel<false>(parent)
    {
        gene.multi_lines(false).reset(gene_name);
        forw.tip_string("Seq. forward primer").multi_lines(false).reset(fw);
        rev .tip_string("Seq. reverse primer").multi_lines(false).reset(rv);
        input_file.tip_string("fasta file with primers").multi_lines(false);
        set.events().click([&]()
        {
            nana::filebox fb{*this, true};
            fb.title("Select a fasta file with primer sequences for gene " + this->gene.text());
            fb.add_filter("fasta file", "*.fasta");
            const auto&files = fb.show();
            if (files.empty()) return;   
            
            this->input_file.reset(files[0].string());
            seqan3::debug_stream << "\nGoing to load: " << files[0].string();
            seqan3::sequence_file_input file_in{files[0]};
            std::string fw, rv;
            int fw_beg{0}, fw_end{0}, rv_beg{0}, rv_end{0};
            for (auto & primer : file_in)
            {
                std::string id = std::move(primer.id());
                seqan3::debug_stream << "\nGoing to check: " << id << "\n" << primer.sequence();

                // parse beg, end from id = >SARS_NF+A -13900 MN908947.3: Seq pos: 28775-28794
                std::string seq_pos = id.substr(id.find("Seq pos: ") + 9);
                size_t dash_pos = seq_pos.find("-");
                int beg = std::stoi(seq_pos.substr(0, dash_pos));
                int end = std::stoi(seq_pos.substr(dash_pos + 1));
                // todo check if beg, end are valid
                seqan3::debug_stream << " with beg: " << beg << " and end: " << end;

                if (beg < end)  // one forward primer/prbe
                {
                    if (fw_beg == 0 || beg < fw_beg)
                    {
                        fw_beg = beg;
                        fw_end = end;
                        auto e = primer.sequence() | seqan3::views::to_char;
                        fw = std::string{e.begin(), e.end()}; 
                        this->forw.reset(fw);
                    }                        
                    continue;
                }
                else  // one reverse primer/prbe
                {
                    if (rv_beg == 0 || beg > rv_beg)
                    {
                        rv_beg = beg;
                        rv_end = end;
                        auto e = primer.sequence() | seqan3::views::to_char;
                        rv = std::string{e.begin(), e.end()}; 
                        this->rev.reset(rv);
                    }
                    continue;
                }
            }
        });

        plc.div(R"( <min=500 all arrange=
                       [27,     50,      200,    200,   30,     variable] gap=7> )");
        plc["all"] << "Gene" << gene << forw << rev << set << input_file;
        plc.collocate();
    }