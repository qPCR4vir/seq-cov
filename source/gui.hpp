
#pragma once

#include <nana/gui.hpp>
#include <nana/gui/widgets/button.hpp>
#include <nana/gui/widgets/group.hpp>
#include <nana/gui/widgets/label.hpp>
#include <nana/gui/widgets/progress.hpp>
#include <nana/gui/widgets/textbox.hpp>
#include <nana/gui/widgets/combox.hpp>
#include <nana/gui/widgets/spinbox.hpp>
#include <nana/gui/place.hpp>

class GeneGUI: public nana::panel<false>
{
  public:
    GeneGUI(nana::window parent, std::string gene_name="", std::string fw="", std::string rv="");

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
    GUI() ;
};

