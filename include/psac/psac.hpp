#pragma once

#include <mpi.h>
#include <vector>
#include <alphabet.hpp>
#include <suffix_array.hpp>
#include <suffix_tree.hpp>
#include <mxx/utils.hpp>

namespace psac {
template <typename return_t = unsigned long>
struct sa_result {
    std::vector<return_t>                 sa;
    std::optional<std::vector<return_t>>  lcp; // present only if build_sa_lcp was called
};

template <typename return_t = unsigned long, typename Container>
[[nodiscard]] sa_result<return_t>
build_sa(MPI_Comm mpi_comm,
         const Container& input,
         std::size_t resolution_threshold)
{
    mxx::comm comm{mpi_comm};
    mxx::print_node_distribution(comm);

    suffix_array<typename Container::value_type, return_t> sa{comm};
    sa.construct(input.begin(), input.end(), resolution_threshold);

    return { std::move(sa.local_SA), std::nullopt };
}

template <typename return_t = unsigned long, typename Container>
[[nodiscard]] sa_result<return_t>
build_sa_lcp(MPI_Comm mpi_comm,
             const Container& input,
             std::size_t resolution_threshold)
{
    mxx::comm comm{mpi_comm};
    mxx::print_node_distribution(comm);

    suffix_array<typename Container::value_type, return_t, /*_CONSTRUCT_LCP=*/true> sa{comm};
    sa.construct(input.begin(), input.end(), resolution_threshold);

    return { std::move(sa.local_SA), std::move(sa.local_LCP) };
}

} // namespace psac
