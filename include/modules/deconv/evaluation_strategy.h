#ifndef EVALUATION_STRATEGY_HH_INCLUDED
#define EVALUATION_STRATEGY_HH_INCLUDED

/**
 * Defines certain evaluation strategies for
 * handling different situations. 
 * If we need to score or filter species 
 * we have different ways we can do that. 
 **/
namespace evaluation_strategy
{
    /**
     * Strategy for scoring species. 
     **/
    enum score_strategy
    {
        /**
         * For integer scoring, the score of 
         * each species is the number of peptides 
         * that species shares a kmer with.
         **/
        INTEGER_SCORING = 0,

        /**
         * For fractional scoring,
         * each species is given a score based 
         * on the ratio of species its shared peptides
         * share a kmer with. So a peptide that shares a kmer
         * with 3 species is given a score of 1/3, and a peptide
         * who only shares a kmer with a given species is given a 
         * score of 1/1. In this way peptides that are more unique to
         * a species are given a higher score.
         **/
        FRACTIONAL_SCORING,

        /**
         * Summation scoring is used in conjuction 
         * with peptides that have count information.
         * So if species1 shares 7 kmers with a peptide,
         * then that species is given a score of 7. 
         * Note that this is done for each peptide a species 
         * shares kmers with and the species is given the score
         * of the sum of these scores.
         **/
        SUMMATION_SCORING
    };

    /**
     * Strategies for filtering peptides, 
     * i.e. if either a species' score or 
     * count falls below a threshold that 
     * species is removed from consideration.
     **/
    enum filter_strategy
    {
        /**
         * Species whose scores fall below
         * a certain threshold will be removed.
         **/
        SCORE_FILTER = 0,

        /**
         * Species whose count falls 
         * below a certain threshold will be 
         * removed.
         **/
        COUNT_FILTER
    };

    /**
     * Strategies specifying how to break ties between species. 
     * Generally one of two outcomes will result from a tie:
     * -Two species are tied and have enough overlapped peptides to be reported 
     *  together. 
     * -Two or more species are tied and the species that shares 
     * the smallest overlap with others is reported.
     **/
    enum tie_eval_strategy
    {
        /**
         * Species that share a certain percentage
         * of their kmers are reported together, or 
         * if they do not meet a given threshold they are reported
         * separately.
         **/
        PERCENT_TIE_EVAL = 0,

        /**
         * Species that share a certain number of 
         * peptides are considered tied.
         **/
        INTEGER_TIE_EVAL,

        /**
         * Species must have a number of 
         * unique shared peptides
         **/
        SUMMATION_SCORING_TIE_EVAL
    };

    enum tie_distance_strategy
    {
        INTEGER_DISTANCE = 0,
        RATIO_DISTANCE
    };

}; //namespace evaluation_strategy

#endif // EVALUATION_STRATEGY_HH_INCLUDED
