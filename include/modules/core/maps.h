#include <parallel_hashmap/phmap.h>
#include <mutex>

#define NMSP phmap
#define MAPNAME phmap::parallel_flat_hash_map
    using phmap::parallel_flat_hash_map;
#define EXTRAARGS_PMAP , NMSP::container_internal::hash_default_hash<K>, \
        NMSP::container_internal::hash_default_eq<K>, \
        std::allocator<std::pair<const K, V>>, 8, std::mutex


#define EXTRAARGS_SMAP , NMSP::container_internal::hash_default_hash<K>, \
        NMSP::container_internal::hash_default_eq<K>, \
        std::allocator<std::pair<const K, V>>, 1, NMSP::NullMutex

#define EXTRAARGS_PSET , NMSP::container_internal::hash_default_hash<K>, \
        NMSP::container_internal::hash_default_eq<K>, \
        std::allocator<K>, 8, std::mutex

#define EXTRAARGS_SSET , NMSP::container_internal::hash_default_hash<K>, \
        NMSP::container_internal::hash_default_eq<K>, \
        std::allocator<K>, 1, NMSP::NullMutex

#define EXTRAARGS_PFSET , NMSP::container_internal::hash_default_hash<K>, \
        NMSP::container_internal::hash_default_eq<K>, \
        std::allocator<K>, 8, std::mutex

template <class K, class V>
using parallel_map = MAPNAME< K, V EXTRAARGS_PMAP>;

template <class K, class V>
using sequential_map = MAPNAME< K, V EXTRAARGS_SMAP>;

template <class K>
using parallel_node_set = phmap::parallel_node_hash_set< K EXTRAARGS_PFSET>;

template <class K>
using parallel_set = phmap::parallel_flat_hash_set< K EXTRAARGS_PSET>;

template <class K>
using sequential_set = phmap::parallel_flat_hash_set< K EXTRAARGS_SSET>;
