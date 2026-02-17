// include MPI
//#include <mpi.h>

// C++ includes
#include <fstream>
#include <iostream>
#include <string>

// distributed suffix array construction
#include <suffix_array.hpp>
#include <check_suffix_array.hpp>
#include <alphabet.hpp>

// suffix tree construction
#include <suffix_tree.hpp>
#include <check_suffix_tree.hpp>

// parallel file block decompose
#include <mxx/env.hpp>
#include <mxx/comm.hpp>
#include <mxx/file.hpp>
#include <mxx/utils.hpp>
// Timer
#include <mxx/timer.hpp>

typedef uint64_t index_t;

struct sa_builder {
public:
    sa_builder() = default;

    void construct_sa_lcp(const MPI_Comm &mpi_com, std::string &local_str) {

        mxx::comm comm = mxx::comm(mpi_com);
        mxx::print_node_distribution(comm);

        mxx::timer t;
        double start = t.elapsed();

        // construct SA+LCP
        suffix_array<char, index_t, true> sa(comm);
        sa.construct(local_str.begin(), local_str.end(), true);
        double end = t.elapsed() - start;
        if (comm.rank() == 0)
            std::cerr << "PSAC time: " << end << " ms" << std::endl;

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
