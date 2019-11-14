#ifndef SPECIES_ID_HH_INCLUDED
#define SPECIES_ID_HH_INCLUDED

/**
 * A class to represent the ID for a species.
 **/
template <class IDType>
class species_id
{
 public:
    /**
     * Default constructor.
     **/
    species_id() = default;

    /**
     * Argument constructor.
     * @param init_id The id to initialize this member with.
     **/
    species_id( const IDType& init_id )
        : id( init_id ) {}

    /**
     * Return a mutable reference to this member's id.
     **/
    IDType& get_id()
        {
            return id;
        }

    /**
     * Return a constant reference to this member's id.
     **/
    const IDType& get_id() const
    {
        return id;
    }

    /**
     * Set this class member's id to new_id.
     * @param new_id The new id for this class member.
     **/
    void set_id( const IDType& new_id )
    {
        id = new_id;
    }

 private:
    IDType id;


};

#endif // SPECIES_ID_HH_INCLUDED
