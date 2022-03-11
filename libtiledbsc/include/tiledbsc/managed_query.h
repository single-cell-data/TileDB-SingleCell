#ifndef TILEDBSC_MANAGED_QUERY_H
#define TILEDBSC_MANAGED_QUERY_H

#include <memory>
#include <vector>

#include <tiledb/tiledb>
#include <tiledbsc/buffer_set.h>

#include "tiledbsc_export.h"


namespace tiledbsc {

// TODO use allocator which does not zero initialize
// TODO config - 4 MB
constexpr size_t TILEDBSC_DEFAULT_ALLOC = 524288;

struct BufferGroup {
    std::map<std::string, std::shared_ptr<BufferSet>> buffers;
};

class TILEDBSC_EXPORT QueryResult {
public:
    QueryResult(std::shared_ptr<BufferGroup> buffers) : buffers_(buffers) {};
    std::shared_ptr<BufferGroup> buffers();
    void to_arrow(void* schema, void* array);

    std::optional<std::shared_ptr<BufferSet>> get(std::string name);

    size_t nbuffers();

    std::vector<std::string> names();

private:
    std::shared_ptr<BufferGroup> buffers_;

};

// Forward declaration
class MQAux;

/**
 *
 * Class encapsulating a TileDB Query with all buffers and batching
 * handled internally.
 *
 */
class TILEDBSC_EXPORT ManagedQuery {

public:
    ManagedQuery(
        const std::shared_ptr<tiledb::Array> array,
        size_t initial_alloc = TILEDBSC_DEFAULT_ALLOC);

    ~ManagedQuery();

    std::unique_ptr<QueryResult> execute();

    /** TODO subset on attributes **/
    void select_attributes(std::vector<std::string> attrs);

    /** TODO subset on dimensions **/
    void select_dimensions(std::vector<std::string> attrs);

    /**
     * Select a set of points to query on a given dimension
     */
    template <typename T>
    void select_points(uint32_t dim_idx, std::vector<T> points) {
        for (auto r = points.begin(); r != points.end(); r++) {
            query_->add_range<T>(dim_idx, *r, *r);
        }
    }

    /**
     * Select a set of (start,end) ranges to query on a given dimension
     */
    template <typename T>
    void select_ranges(uint32_t dim_idx, std::vector<std::array<T, 2>> ranges) {
        for (auto r = ranges.begin(); r != ranges.end(); r++) {
            query_->add_range<T>(dim_idx, *r[0], *r[1]);
        }
    }

    std::unique_ptr<tiledb::Query> query_;

private:
    /* methods */
    static std::shared_ptr<tiledb::Array> check_array(std::shared_ptr<tiledb::Array>);

    void allocate_buffers();
    void set_buffers();
    void resize_result_buffers();
    void validate_query();
    void complete_query();

    /* fields */
    std::shared_ptr<BufferGroup> buffers_;

    std::shared_ptr<tiledb::Array> array_;

    size_t initial_alloc_;

    /***/
    // if set, use only attribute or dim with these names
    std::optional<std::set<std::string>> use_attrs_dims_ = std::nullopt;

    /* private implementation for helpers */
    friend class MQAux;
    std::unique_ptr<MQAux> impl_;
};

};

#endif