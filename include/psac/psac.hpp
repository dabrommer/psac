#pragma once

#include <mpi.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <alphabet.hpp>
#include <check_suffix_array.hpp>
#include <check_suffix_tree.hpp>
#include <suffix_array.hpp>
#include <suffix_tree.hpp>

#include <mxx/comm.hpp>
#include <mxx/env.hpp>
#include <mxx/file.hpp>
#include <mxx/timer.hpp>
#include <mxx/utils.hpp>

struct sa_builder {
public:
    sa_builder() = default;

    void construct_sa_lcp(const MPI_Comm& mpi_com, std::string& local_str, std::size_t resolval_threshold) {
        mxx::comm comm = mxx::comm(mpi_com);
        mxx::print_node_distribution(comm);

        mxx::timer t;
        double start = t.elapsed();

        // construct SA+LCP
        using index_t = uint64_t;
        suffix_array<char, index_t, true> sa(comm);
        sa.construct(local_str.begin(), local_str.end(), resolval_threshold);
        double end = t.elapsed() - start;
        if (comm.rank() == 0) {
            std::cerr << "PSAC time: " << end << " ms" << std::endl;
        }

        sa_vec = sa.local_SA;
        lcp_vec = sa.local_LCP;
    }

    std::vector<unsigned long>& get_sa() {
        return sa_vec;
    }

    std::vector<unsigned long>& get_lcp() {
        return lcp_vec;
    }

private:
    std::vector<unsigned long> sa_vec;
    std::vector<unsigned long> lcp_vec;
};
