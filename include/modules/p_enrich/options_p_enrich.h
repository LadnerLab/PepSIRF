#ifndef OPTIONS_P_ENRICH_HH_INCLUDED
#define OPTIONS_P_ENRICH_HH_INCLUDED
#include "options.h"

#include <string>

class options_p_enrich : public options
{
 public:

    options_p_enrich()
        : DEFAULT_OUTPUT_DIRNAME{ "paired" },
          DEFAULT_OUT_FNAME_JOIN{ "~" }
          {}

    const std::string DEFAULT_OUTPUT_DIRNAME;
    const std::string DEFAULT_OUT_FNAME_JOIN;
};

#endif // OPTIONS_P_ENRICH_HH_INCLUDED
