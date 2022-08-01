#include "tiledbsc/soma_collection_query.h"
#include "tiledbsc/logger_public.h"
#include "tiledbsc/soma_collection.h"

#include "thread_pool/thread_pool.h"

namespace tiledbsc {
using namespace tiledb;

SOMACollectionQuery::SOMACollectionQuery(SOMACollection* soco) {
    std::vector<ThreadPool::Task> tasks;
    ThreadPool pool{threads_};

    for (auto& [name, soma] : soco->get_somas()) {
        LOG_DEBUG(fmt::format("Get SOMA query: {}", name));
        soma_queries_[name] = soma->query(name);
    }
}

std::optional<SOCOBuffers> SOMACollectionQuery::next_results() {
    submitted_ = true;

    LOG_DEBUG(fmt::format("[SOMACollectionQuery] Start queries."));

    std::vector<ThreadPool::Task> tasks;
    ThreadPool pool{threads_};

    for (auto& [name, sq] : soma_queries_) {
        if (!sq->is_complete()) {
            LOG_DEBUG(
                fmt::format("[SOMACollectionQuery] Queue query for {}", name));
            tasks.emplace_back(pool.execute([&]() {
                sq->next_results();
                return Status::Ok();
            }));
        }
    }

    // Block until all tasks complete
    pool.wait_all(tasks).ok();

    LOG_DEBUG(fmt::format("[SOMACollectionQuery] Queries done."));

    SOCOBuffers results;

    for (auto& [name, sq] : soma_queries_) {
        if (sq->results().has_value()) {
            LOG_DEBUG(fmt::format(
                "[SOMACollectionQuery] SOMA {} has {} results.",
                name,
                sq->results()->begin()->second.begin()->second->size()));
            results[name] = *sq->results();
        }
    }

    if (results.empty()) {
        return std::nullopt;
    }
    return results;
}

/*
std::optional<SOCOBuffers> SOMACollectionQuery::next_results() {
    submitted_ = true;

    SOCOBuffers results;

    for (auto& [name, sq] : soma_queries_) {
        LOG_DEBUG(fmt::format("[SOMACollectionQuery] Run query for {}", name));
        auto soma_results = sq->next_results();

        if (soma_results.has_value()) {
            results[name] = *sq->next_results();
        }
    }

    LOG_DEBUG(fmt::format(
        "[SOMACollectionQuery] Queries done. SOMA result count = {}",
        results.size()));

    if (results.empty()) {
        return std::nullopt;
    }

    return results;
}
*/

}  // namespace tiledbsc
