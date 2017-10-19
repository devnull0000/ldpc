#pragma once
#include <memory>
#include <vector>
#include <assert.h>

namespace ldpc {

class tanner_graph_t
{
public:
        struct check_node_t;

        struct bit_node_t
        {
                std::vector<check_node_t *>     check_nodes;
                std::vector<int>                my_idx_in_cn;           //index of this bit_node in check_nodes[i].bit_nodes
                std::vector<float>              in_estimates;
                std::vector<float>              out_estimates;
                float                           apriori_prob = 0.f;
        };
        typedef std::shared_ptr<bit_node_t>     bit_node_pt;

        struct check_node_t
        {
                std::vector<bit_node_t *>       bit_nodes;
                std::vector<int>                my_idx_in_bn;
        };
        typedef std::shared_ptr<check_node_t>   check_node_pt;

        template <int c_check_bits, int c_max_hc>
        tanner_graph_t(int const (&sparsed_h_matrix)[c_check_bits][c_max_hc]);

        void            reset(float const * const probs);
        bit_node_t *    bit_node(int const column_idx);
        check_node_t *  check_node(int const row_idx);

private:
        std::vector<check_node_pt>      m_rows;
        std::vector<bit_node_pt>        m_columns;
};

} //ldpc


#include "ldpc_tanner_graph.impl.hpp"
