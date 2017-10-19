#pragma once
#include <limits.h>
#include <functional>

//random matrix generator with specified column weight
namespace ldpc
{
namespace v02
{

static constexpr int c_n = 144;         //number of bits in codeword
static constexpr int c_k = 64;          //number of bits in message (useful, which we encode)
static constexpr int c_m = c_n - c_k;       //number of check bits

static constexpr int c_wc = 5;          //weight of column (number of 1s in each column)
//static constexpr int c_wr = 5;          //approximate number of 1s in each row. can be changed, but won't be less than c_min_wr
constexpr int c_min_wr = 2;             //minimum weight of a row

//removes all loops up to loop_length (inclusive). i.e. girth of the graph will be not less then loop_length+1
//there is a possibility that the task is impossible...
constexpr int c_loop_length = 8;

constexpr int64_t c_max_loop_iters = 55; //INT_MAX;   //maximum number of iterations to try to remove loops


inline void gen_H_matrix(bool (&mx)[c_m][c_n])
{
        memset(mx, 0, sizeof(mx));
        int row_idx = 0;
        for (int col_idx = 0; col_idx < c_n; ++col_idx) {
                for (int sub_row_idx = 0; sub_row_idx < c_wc; ++sub_row_idx)
                        mx[(row_idx + sub_row_idx)%c_m][col_idx] = true;
                row_idx += c_wc;
        }
}

typedef int16_t ms_t;       //measure type; should be > c_n; keep as small as possible
struct point_t
{
        ms_t r, c;      //row, column
};
inline bool operator== (point_t const l, point_t const r) { return l.r == r.r && l.c == r.c; }
inline bool operator!= (point_t const l, point_t const r) { return l.r != r.r || l.c != r.c; }

/*
typedef struct loop_t
{
        point_t p[c_loop_length];  //maximum which can be detected
        int     length = 0;        //actual
};*/

struct loops_t
{
        uint16_t num[c_loop_length / 2 - 1];         //only even numbers; skip 0 & 1 indices
};

void operator+= (loops_t & l, loops_t const & r) {
        for (int i = 0; i < (c_loop_length / 2 - 1); ++i)
                l.num[i] += r.num[i];
}


//returns length of the found loop, or INT_MAX if no loop were found; maximum depth examined c_loop_length
inline int find_loop_bit(bool const (&mx)[c_m][c_n], loops_t & loops, point_t const start,
                           ms_t const from_row, ms_t const to_column, int const len);


inline int find_loop_check(bool const (&mx)[c_m][c_n], loops_t & loops, point_t const start,
                                ms_t const from_column, ms_t const to_row, int const len)
{
        if (to_row== start.r) {
                loops.num[len/2 - 2]++;
                return len;
        }
        if (len >= c_loop_length)
                return INT_MAX;

        int min_loop = INT_MAX;
        for (ms_t j = 0; j < c_n; ++j) {
                if (mx[to_row][j] && j != from_column) {
                        int const l = find_loop_bit(mx ,loops, start, to_row, j, len + 1);
                        min_loop = std::min(min_loop, l);
                }
        }
        return min_loop;
}

int find_loop_bit(bool const (&mx)[c_m][c_n], loops_t & loops, point_t const start,
                  ms_t const from_row, ms_t const to_column, int const len)
{
        int min_loop = INT_MAX;
        for (ms_t i = 0; i < c_m; ++i) {
                if (mx[i][to_column] && i != from_row) {
                        int const l = find_loop_check(mx, loops, start, to_column, i, len + 1);
                        min_loop = std::min(min_loop, l);
                }
        }
        return min_loop;
}

inline int find_loop(bool const (&mx)[c_m][c_n], loops_t & loops, point_t const start)
{
        //we are in bit node, go to check ones
        memset(&loops, 0, sizeof(loops));
        int const ret = find_loop_bit(mx, loops, start, start.r, start.c, 1);
        return ret;
}


inline bool can_move_1_from_here(bool (&mx)[c_m][c_n], point_t const pt)
{
        int wr = 0;
        for (int j = 0; j < c_n; ++j)
                wr += mx[pt.r][j];
        return (wr - 1) >= c_min_wr;
}

#if 1
inline bool is_better(loops_t const & ol, loops_t const & nl, bool & can_shake)
{
        for (int i = 0; i < sizeof(ol.num)/sizeof(ol.num[0]); ++i) {
                if (nl.num[i] < ol.num[i])
                        return true;
                else if (nl.num[i] > ol.num[i])
                        return false;
        }
        return false;
}
#else
inline bool is_better(loops_t const & ol, loops_t const & nl, bool & can_shake)
{
        if (can_shake) { // && nl.num[sizeof(ol.num)/sizeof(ol.num[0]) - 2] == 0) {
                can_shake = false;
                return true;    //try to move from local non optimum
        }

        for (int i = 0; i < sizeof(ol.num)/sizeof(ol.num[0]); ++i) {
                int dist = ol.num[i] - nl.num[i];
                if ( dist > 0 )
                        return true;
                else if (dist < 0) {
                        return false;
                }
        }
        return false;
}
#endif

inline bool handle_loop(bool (&mx)[c_m][c_n], point_t const loop_pt, loops_t const & loops, bool & can_shake)
{
        int start = rand();
        bool was_enhancment = false;
        for (int ci = 0; ci < c_m; ++ci) {
                int i = (start + ci) % c_m;
                if (!mx[i][loop_pt.c] && i != loop_pt.r && can_move_1_from_here(mx, loop_pt)) {
                        //try to unleash the loop
                        mx[loop_pt.r][loop_pt.c] = false;
                        mx[i][loop_pt.c] = true;

                        loops_t new_loops;
                        find_loop(mx, new_loops, point_t{(ms_t)i, loop_pt.c});
                        if (!is_better(loops, new_loops, can_shake)) { //we didn't enhance the situation, so rollback
                                mx[loop_pt.r][loop_pt.c] = true;
                                mx[(ms_t)i][loop_pt.c] = false;
                        }
                        else {
                                was_enhancment = true;
                                break;
                        }
                }
        }
        return was_enhancment;
}

typedef std::function<void (loops_t const &)> loop_stat_fn_t;
inline void remove_loops(bool (&mx)[c_m][c_n], loop_stat_fn_t const & stat_fn)
{
        int no_enhancement_counter = 0;
        bool can_shake = false;

        for(int64_t num_iter = 0; num_iter < c_max_loop_iters; ++num_iter) {
                bool was_enhancement = false;
                int num_loops = 0;
                loops_t all_loops;
                memset(&all_loops, 0, sizeof(all_loops));

                for (int i = 0; i < c_m; ++i) {
                        int rand_j = rand();
                        for (int wj = 0; wj < c_n; ++wj) {
                                int const j = (rand_j + wj) % c_n;
                                if (mx[i][j]) {
                                        loops_t loops;
                                        int const loop = find_loop(mx, loops, point_t{(ms_t)i, (ms_t)j});
                                        if (loop != INT_MAX) {  //we found a loop!!!
                                                //try to remove it, we can do this by moving 1 higher or below of where we are
                                                //but only in case the place where we move is free and if we don't create another loop
                                                //of the same length or shorter
                                                was_enhancement |= handle_loop(mx, point_t{(ms_t)i, (ms_t)j}, loops, can_shake);
                                                all_loops.num[loop/2-2]++;
                                                num_loops++;
                                        }
                                }
                        }
                        //printf("iteration: %03d row: %03d\n", (int)num_iter, i);
                }

                if (!was_enhancement) { ++no_enhancement_counter; can_shake = true; }
                else no_enhancement_counter = 0;

                //print_matrix(mx);
                printf("iteration: %d was_enhancement: %d total loops: %d\n", (int)num_iter, (int)was_enhancement, num_loops);
                stat_fn(all_loops);


                if (!num_loops || no_enhancement_counter > 3) //mission completed!
                        break;
        }

}





}
}
