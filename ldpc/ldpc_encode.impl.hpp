//! do not include this file directly

namespace ldpc {
    
template<int c_n, int c_k>
void
encode_t<c_n, c_k>::encode(bool (&cw)[c_n], bool const (&msg)[c_k])
{
        for(int i = 0; i < (c_n - c_k); ++i)
                cw[i] = bit_ops::mad_rows(m_g[i], msg);
        memcpy(&cw[c_n - c_k], msg, sizeof(msg));
        bit_ops::column_reorder_enc(cw, m_column_idx);
}

template<int c_n, int c_k>
void
encode_t<c_n, c_k>::encode(uint8_t (&r)[c_codeword_bytes], uint8_t const (&msg)[c_message_bytes])
{
        bool msg_bits[c_k];
        bool cw_bits[c_n];
        ldpc::bit_ops::convert_from_bytes(msg_bits, msg);
        this->encode(cw_bits, msg_bits);
        ldpc::bit_ops::convert_from_bit_array(r, cw_bits);
}
        
template<int c_n, int c_k>
std::vector<uint8_t>
encode_t<c_n, c_k>::encode(std::vector<uint8_t> const & data)
{
        int const seg_length = bit_ops::align_division_up(c_k,8);
        int const seg_num = (int)data.size() / seg_length;
        int const seg_num_au = bit_ops::align_division_up((int)data.size(), seg_length);
        
        int const out_seg_length = bit_ops::align_division_up(c_n, 8);
        
        std::vector<uint8_t> ret; ret.resize(seg_num_au * out_seg_length);
        for (int si = 0; si < seg_num; ++si) {
                this->encode(*(uint8_t (*)[bit_ops::align_division_up(c_n,8)])&ret[si*out_seg_length], *(uint8_t const (*)[bit_ops::align_division_up(c_k,8)])&data[si*seg_length]);
        }
        
        if (seg_num != seg_num_au) {
                //encode last message with zero-padding
                uint8_t msg[bit_ops::align_division_up(c_k,8)] = {0};
                memcpy(msg, data.data() + seg_num * seg_length, (int)data.size() - seg_num * seg_length);
                this->encode(*(uint8_t (*)[bit_ops::align_division_up(c_n,8)])&ret[seg_num*out_seg_length], msg);
        }
        return ret;
}
        
}
