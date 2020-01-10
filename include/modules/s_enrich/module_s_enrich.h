#ifndef MODULE_S_ENRICH_HH_INCLUDED
#define MODULE_S_ENRICH_HH_INCLUDED
#include "module.h"
#include "options_s_enrich.h"
#include "peptide_score.h"
#include "peptide_scoring.h"
#include <vector>
#include <fstream>

/**
 * Base case for value_constrained_by
 **/
template<typename Val, typename UnaryPred>
bool value_constrained_by( const Val& arg, const UnaryPred funct )
{
    return funct( arg );
}

/**
 * Template to determine whether a value is constrained by 
 * each in a list of unary predicates. 
 * @tparam Val the type of the value to test
 * @tparam UnaryPred the type of the current unary predicate to test, must 
 *         be a function compatible with  the signature 
 *         []( const Val val ) -> bool
 * @tparam preds the rest of the predicates to test
 * @returns true if every predicate returns true given the value arg
 **/
template<typename Val, typename UnaryPred, typename... Predicates>
bool value_constrained_by( const Val& arg,
                           const UnaryPred funct,
                           Predicates... preds
                         )
{
    constexpr const bool is_same = std::is_same<decltype( funct( arg ) ), bool>::value;

    static_assert( is_same, "funct must be a functor that takes value of type Val "
                            "and returns bool."
                 );

    return funct( arg ) && value_constrained_by( arg, preds... );
}


/**
 * Capture the values in the range [begin, end) for which 
 * pred returns true, passing the return value of get to pred.
 * @tparam InputIterator The type of the iterator used as input
 * @tparam OutputIterator the type of the iterator used as output
 * @tparam Pred the predicate used to perform the test for each value
 * @tparam Get A function that taks a value of type InputIterator::value_type
 *         and returns a member of a type for which pred is defined.
 * @pre Pred has a signature 
 *      equivalent to pred( InputIterator::value_type )->bool
 * @param begin The first item in the range [begin,end)
 * @param end The last item in the range [begin,end)
 * @param dest The location to store the values for which pred is true,
 * @note Behavoir undefined if dest is not large enoughto hold all values for which 
 *       pred is true.
 * @returns OutputIterator pointing after the last item that was evaluated 
 *          as true 
 **/
template<typename InputIterator,
         typename OutputIterator,
         typename Pred,
         typename Get
         >
OutputIterator valid_for( const InputIterator begin,
                          const InputIterator end,
                          OutputIterator dest,
                          const Pred pred,
                          Get get
                        )
{
    using iter_val = typename InputIterator::value_type;
    using get_ret_type = decltype( get( std::declval<iter_val>() ) );
    using pred_ret_type = decltype( pred( std::declval<get_ret_type>() ) );

    static_assert( std::is_same<pred_ret_type,
                                bool>::value,
                   "Pred must be a unary predicate that returns bool."
                 );

    for( auto curr = begin; curr != end; ++curr )
        {
            if( pred( get( *curr ) ) )
                {
                    *dest = *curr;
                    ++dest;
                }
            
        }
    return dest;
}

/**
 * Capture the values in the range [begin, end) for which 
 * pred returns true.
 * @tparam InputIterator The type of the iterator used as input
 * @tparam OutputIterator the type of the iterator used as output
 * @tparam Pred the predicate used to perform the test for each value
 * @pre Pred has a signature 
 *      equivalent to pred( InputIterator::value_type )->bool
 * @param begin The first item in the range [begin,end)
 * @param end The last item in the range [begin,end)
 * @param dest The location to store the values for which pred is true,
 * @note Behavoir undefined if dest is not large enoughto hold all values for which 
 *       pred is true.
 * @returns OutputIterator pointing after the last item that was evaluated 
 *          as true 
 **/
template<typename InputIterator,
         typename OutputIterator,
         typename Pred
         >
OutputIterator valid_for( const InputIterator begin,
                          const InputIterator end,
                          OutputIterator dest,
                          const Pred pred
                        )
{
    using iter_val = typename InputIterator::value_type;
    auto get = []( const iter_val& val ) -> const iter_val& { return val; };
    return valid_for( begin, end, dest, pred, get );
}


class module_s_enrich : public module
{
public:

    void run( options *opts );

    /**
     * Get the score data about a peptide, returning each of its scores in 
     * a vector of peptide_score types. 
     * @param zscore_data contains data for peptide zscores.
     * @param norm_score_data Data containing normalized scores for each peptide
     * @param raw_score_data Data containing raw scores for each peptide.
     * @pre zscore_data, norm_score_data, and (if included) raw_score_data 
     *      all contain data for the same peptides and samples.
     * @pre zscore_data and norm_score_data are not null
     * @note if raw_score_data is null, the return peptide_scores will 
     *       have their raw_score member set to 0.
     **/
    std::vector<peptide_score<std::string>>
    get_enrichment_candidates( const peptide_score_data_sample_major *zscore_data,
                               const peptide_score_data_sample_major *norm_score_data,
                               const peptide_score_data_sample_major *raw_score_data,
                               const std::string sample_id
                             );

    /**
     * Write the names of probes to an output stream, one per line..
     * @param stream The stream to write probe names to
     * @param probes A vector of peptide_scores to whose names
     *        should be written to the stream.
     **/
    void write_probe_names( std::ostream &stream,
                            const std::vector<peptide_score<std::string>>& probes
                          );


};

#endif // MODULE_S_ENRICH_HH_INCLUDED 
