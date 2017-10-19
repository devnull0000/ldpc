//do not include directly!

namespace ldpc {
        
template<int c_n, int c_k, int c_max_hc>
bool
decode_t<c_n, c_k, c_max_hc>::decode_bsc_hard(bool (&msg)[c_k], bool const (&codeword_orig)[c_n], int const num_iterations)
{
        bool codeword[c_n];
        memcpy(codeword, codeword_orig, sizeof(codeword));

        //be optimistic and check initially, probably everything is just fine right now!
        if (this->check_codeword(codeword)) {
                bit_ops::column_reorder_dec(codeword, m_column_idx);
                memcpy(msg, &codeword[c_n-c_k], sizeof(msg));
                return true;
        }

        for (int iter = 0; iter < num_iterations; ++iter) {
                uint16_t vote_for_zero[c_n] = {0};
                uint16_t vote_for_one[c_n] = {0};
                for (int ri = 0; ri < c_check_bits; ++ri) {     //ri == h matrix row index
                        for (int ci = 0; ci < c_max_hc; ++ci) {
                                int const cw_idx = m_h[ri][ci];
                                if (cw_idx == -1)
                                        break;

                                bool val = false;
                                for (int ci2 = 0; ci2 < c_max_hc; ++ci2) {
                                        int const cw_idx2 = m_h[ri][ci2];
                                        if (cw_idx2 == -1)
                                                break;
                                        if (ci != ci2)
                                                val = val != codeword[cw_idx2];
                                }
                                if (val)
                                        ++vote_for_one[cw_idx];
                                else
                                        ++vote_for_zero[cw_idx];
                        }
                }
                for (int ci = 0; ci < c_n; ++ci)
                        if (vote_for_zero[ci] > vote_for_one[ci]) {
                                if( codeword[ci] ) {
                                        codeword[ci] = false;
                                        //break;
                                }
                        }
                        else {
                                if( !codeword[ci] ) {
                                        codeword[ci] = true;
                                        //break;
                                }
                        }

                if (this->check_codeword(codeword)) {
                        bit_ops::column_reorder_dec(codeword, m_column_idx);
                        memcpy(msg, &codeword[c_n-c_k], sizeof(msg));
                        return true;
                }
        }

        //return best of what we had (and it's an original value)
        bit_ops::column_reorder_dec(codeword, m_column_idx);
        memcpy(msg, &codeword[c_n-c_k], sizeof(msg));
        return false;
}

template<int c_n, int c_k, int c_max_hc>
bool
decode_t<c_n, c_k, c_max_hc>::decode_bsc_hard(uint8_t (&msg)[c_message_bytes], uint8_t (&codeword)[c_codeword_bytes], int const num_iterations)
{
        bool msg_bits[c_k];
        bool cw_bits[c_n];
        bit_ops::convert_from_bytes(cw_bits, codeword);
        bool const r = this->decode_bsc_hard(msg_bits, cw_bits, num_iterations);
        bit_ops::convert_from_bit_array(msg, msg_bits);
        return r;
}

template<int c_n, int c_k, int c_max_hc>
bool
decode_t<c_n, c_k, c_max_hc>::check_codeword(bool const (&codeword)[c_n])
{
        for (int i = 0; i < c_check_bits; ++i) {
                bool const val = bit_ops::mad_rows_sparse(m_h[i], codeword);
                if (val)
                        return false;
        }
        return true;
}

//*********************         decode_soft_t           ****************************************
        
template<int c_n, int c_k, int c_max_hc>
void
decode_soft_t<c_n, c_k, c_max_hc>::compute_out_estimates(check_node_t * const cn)
{
        int const num_bn = (int) cn->bit_nodes.size();
        for( int lnode_idx = 0; lnode_idx < num_bn; ++lnode_idx ) {       //left_node index (in equation p = ...) for which we compute now the probability
                bit_node_t * const lnode = cn->bit_nodes[lnode_idx];
                float llr = 0.f;
                float group = 1.f;
                for (int j = 0; j < num_bn; ++j) {
                        if (j != lnode_idx) {
                                bit_node_t * const bj = cn->bit_nodes[j];
                                float const p = bj->out_estimates[cn->my_idx_in_bn[j]];
                                float t = tanhf(p/2.f);
                                group *= t;
                        }
                }
                
                if (group == 1.f)
                        group = 0.999999f;
                if (group == -1.f)
                        group = -0.999999f;
                
                llr = log((1.f + group) / (1.f - group));
                if (std::isnan(llr) || std::isinf(llr)) {
                        assert(false);
                        int a = 10; ++a; (void)a;
                }
                
                lnode->in_estimates[cn->my_idx_in_bn[lnode_idx]] = llr;
        }
}


template<int c_n, int c_k, int c_max_hc>
void
decode_soft_t<c_n, c_k, c_max_hc>::compute_in_estimates(bit_node_t * const bn)
{
        int const num_cn = (int)bn->check_nodes.size();
        for( int cio = 0; cio < num_cn; ++cio ) {               //cio = Check node Index Output
                float prob = 0.f;
                for (int cii = 0; cii < num_cn; ++cii)
                        if (cio != cii)
                                prob += bn->in_estimates[cii];
                bn->out_estimates[cio] = bn->apriori_prob + prob;
        }
        
}

template<int c_n, int c_k, int c_max_hc>
bool
decode_soft_t<c_n, c_k, c_max_hc>::check_estimates_llr(bool (&msg)[c_k])
{
        float estimates[c_n] = {0};
        for( int ci = 0; ci < c_n; ++ci ) {
                bit_node_t * const bn = m_graph.bit_node(ci);
                float prob = 0.f;
                int const num_cn = (int) bn->check_nodes.size();
                for( int cio = 0; cio < num_cn; ++cio )               //cio = Check node Index Output
                        prob += bn->in_estimates[cio];
                
                estimates[ci] = bn->apriori_prob + prob;
        }
        bool const ret = this->check_codeword_and_write_llr(msg, estimates);
        return ret;
}
        
template<int c_n, int c_k, int c_max_hc>
bool
decode_soft_t<c_n, c_k, c_max_hc>::decode_bsc_soft(bool (&msg)[c_k], float const (&probs)[c_n], int const num_iterations)
{
        if( this->check_codeword_and_write(msg, probs) )
                return true;
        
        m_graph.reset(probs);
        
        for( int iter = 0; iter < num_iterations; ++iter ) {
                for( int ri = 0; ri < c_check_bits; ++ri ) {
                        check_node_t * const cn = m_graph.check_node(ri);
                        this->compute_out_estimates(cn);
                }

                if (this->check_estimates_llr(msg))
                        return true;

                for( int ci = 0; ci < c_n; ++ci ) {
                        bit_node_t * const bn = m_graph.bit_node(ci);
                        this->compute_in_estimates(bn);
                }
        }

        {
                bool codeword[c_n];
                bit_ops::rasterize_probs(codeword, probs);
                bit_ops::column_reorder_dec(codeword, this->m_column_idx);
                memcpy(msg, &codeword[c_n - c_k], sizeof(msg));
        }
        return false;
}
        
template<int c_n, int c_k, int c_max_hc>
bool
decode_soft_t<c_n, c_k, c_max_hc>::decode_bsc_soft(uint8_t (&msg)[c_message_bytes], float const (&probs)[c_n], int const num_iterations)
{
        bool msg_bits[c_k];
        bool const r = this->decode_bsc_soft(msg_bits, probs, num_iterations);
        bit_ops::convert_from_bit_array(msg, msg_bits);
        return r;
}

template<int c_n, int c_k, int c_max_hc>
bool
decode_soft_t<c_n, c_k, c_max_hc>::check_codeword_and_write_llr(bool (&msg)[c_k], float const (&llrs)[c_n])
{
        bool codeword[c_n];
        bit_ops::rasterize_llr(codeword, llrs);
        if (this->check_codeword(codeword)) {
                bit_ops::column_reorder_dec(codeword, this->m_column_idx);
                memcpy(msg, &codeword[c_n-c_k], sizeof(msg));
                return true;
        } else
                return false;
}

template<int c_n, int c_k, int c_max_hc>
bool
decode_soft_t<c_n, c_k, c_max_hc>::check_codeword_and_write(bool (&msg)[c_k], float const (&probs)[c_n])
{
        bool codeword[c_n];
        bit_ops::rasterize_probs(codeword, probs);
        if (this->check_codeword(codeword)) {
                bit_ops::column_reorder_dec(codeword, this->m_column_idx);
                memcpy(msg, &codeword[c_n-c_k], sizeof(msg));
                return true;
        } else
                return false;
}
        

}
