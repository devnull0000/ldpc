#pragma once

namespace ldpc { namespace bit_ops{
        
template<int c_n>
bool mad_rows(bool const (&a)[c_n], bool const (&b)[c_n])
{
        bool res = false;
        for (int i = 0; i < c_n; ++i)
                res = res != (a[i] && b[i]);
        return res;
}
        
template<int c_m, int c_n>
bool mad_row_col(bool const (&r)[c_n], bool const (&c)[c_m][c_n], int const column_idx)
{
        bool res = false;
        for (int i = 0; i < c_n; ++i) {
                res = res != r[i] && c[i][column_idx];
        }
        return res;
}

template<int c_max_num_cols, int c_n>
bool mad_rows_sparse(int const (&a)[c_max_num_cols], bool const (&b)[c_n])
{
        bool res = false;
        for (int si = 0; si < c_max_num_cols; ++si) {
                int const i = a[si];
                if (i == -1) break;
                bool const bit = b[i];
                res = res != bit;
        }
        return res;
}

template<int c_n>
void column_reorder_enc(bool (&a)[c_n], int const (&column_idx)[c_n])           //or front
{
        bool r[c_n];
        for (int i = 0; i < c_n; ++i)
                r[column_idx[i]] = a[i];
        memcpy(a, r, sizeof(r));
}

template<int c_n>
void column_reorder_dec(bool (&a)[c_n], int const (&column_idx)[c_n])           //or front
{
        bool r[c_n];
        for (int i = 0; i < c_n; ++i)
                r[i] = a[column_idx[i]];
        memcpy(a, r, sizeof(r));
}


constexpr int align_division_up(int const a, int const b) {
        int r = a/b;
        if (r * b == a) return r; else return r+1;
}


//converts any data to boolean array;
//if data size in bits is longer than c_n, the least significant bits are moved to r
template<int c_n>
void convert_from_bytes(bool (&r)[c_n], uint8_t const (&data)[align_division_up(c_n,8)])
{
        for (int byte_idx = 0; byte_idx < align_division_up(c_n, 8); ++byte_idx)
                for (int bit_idx = 0; bit_idx < std::min(8, c_n - byte_idx*8); ++bit_idx)
                        r[byte_idx * 8 + bit_idx] = !!((data[byte_idx] >> bit_idx) & 0x1);
}

//converts boolean array to data back. if there is bit-space in r, most significant beats are filled with zero
template<int c_n>
void convert_from_bit_array(uint8_t (&r)[align_division_up(c_n,8)], bool const (&ba)[c_n])
{
        memset(r, 0, sizeof(r));
        for (int byte_idx = 0; byte_idx < align_division_up(c_n, 8); ++byte_idx)
                for (int bit_idx = 0; bit_idx < std::min(8, c_n - byte_idx*8); ++bit_idx)
                        r[byte_idx] |= int(ba[byte_idx * 8 + bit_idx]) << bit_idx;
}

template <int c_n>
void rasterize_llr(bool (&codeword)[c_n], float const (&llrs)[c_n])
{
        for (int i = 0; i < c_n; ++i)
                codeword[i] = ((llrs[i] <= 0) ? 0.f : 1.f);
}

template <int c_n>
void rasterize_probs(bool (&codeword)[c_n], float const (&probs)[c_n])
{
        for (int i = 0; i < c_n; ++i)
                codeword[i] = ((probs[i] >= 0.5) ? 1.f : 0.f);
}
        

}} //namespaces
