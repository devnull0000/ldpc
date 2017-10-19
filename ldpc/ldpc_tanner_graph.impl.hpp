//do not include directly

namespace ldpc {

template<int c_check_bits, int c_max_hc>
tanner_graph_t::tanner_graph_t(int const (& sparsed_h_matrix)[c_check_bits][c_max_hc])
{
        int num_cols = 0;

        m_rows.resize(c_check_bits);
        for( int ri = 0; ri < c_check_bits; ++ri ) {
                m_rows[ri] = std::make_shared<check_node_t>();
                for( int mi = 0; mi < c_max_hc; ++mi ) {
                        int const col_idx = sparsed_h_matrix[ri][mi];
                        if( col_idx == -1 ) break;
                        num_cols = std::max(num_cols, col_idx);
                }
        }
        num_cols += 1;

        m_columns.resize(num_cols);
        for( int ci = 0; ci < num_cols; ++ci )
                m_columns[ci] = std::make_shared<bit_node_t>();

        for( int ri = 0; ri < c_check_bits; ++ri ) {
                check_node_t * const cn = m_rows[ri].get();
                for( int mi = 0; mi < c_max_hc; ++mi ) {
                        int const col_idx = sparsed_h_matrix[ri][mi];
                        if( col_idx == -1 ) break;

                        bit_node_t * const bn = m_columns[col_idx].get();
                        cn->my_idx_in_bn.push_back((int)bn->check_nodes.size());
                        bn->my_idx_in_cn.push_back((int)cn->bit_nodes.size());
                        cn->bit_nodes.push_back(bn);
                        bn->check_nodes.push_back(cn);
                        bn->out_estimates.push_back(0.f);
                        bn->in_estimates.push_back(0.f);
                }
        }
}

inline float log_ll(float p)
{
        assert(p<=1 && p >= 0);
        if (p == 0) p = 0.000001f;
        else if (p == 1) p = 0.999999f;
        
        float r = logf(p/(1-p));
        if (std::isnan(r)) {
                int a = 10; ++a;
                assert(false);
        }
        return r;
}
        
inline
void tanner_graph_t::reset(float const * const probs) {
        int const num_cols = (int)this->m_columns.size();
        for( int ci = 0; ci < num_cols; ++ci ) {
                bit_node_t * const bn = m_columns[ci].get();
                bn->apriori_prob = log_ll(probs[ci]);

                int const num_estimates = (int)bn->out_estimates.size();
                for( int cni = 0; cni < num_estimates; ++cni )
                        bn->out_estimates[cni] = log_ll(probs[ci]);
        }
}

inline
tanner_graph_t::bit_node_t * tanner_graph_t::bit_node(int const column_idx) {
        return this->m_columns[column_idx].get();
}

inline
tanner_graph_t::check_node_t * tanner_graph_t::check_node(int const row_idx) {
        return this->m_rows[row_idx].get();
}


} //ldpc
