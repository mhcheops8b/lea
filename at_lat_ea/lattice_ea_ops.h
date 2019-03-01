#pragma once

#include "lattice_ea.h"
#include <map>

struct Lattice_EA_ops {

    Lattice_EA_ops(const Lattice_EA& lea)
        : lea(lea)
    {
        init();
    }

    block_based_elem_type idx_to_elem(int idx);
    int elem_to_idx(const block_based_elem_type& elem);

    int get_size();

    const Lattice_EA& get_lea() { return lea; }

    // order relation, operations
    int order(int e1, int el2);
    int oplus(int el1, int el2);
    int orthosupplement(int el);
    int inf(int el1, int el2);
    int sup(int el1, int el2);
    int impl(int el1, int el2);

    // must be called before accessing any operation
    // for large algebras might be slow
    void generate_ops();

private:
    void init();
    int elem_arr_idx(int row, int col);

    std::map<int, block_based_elem_type> map_idx_elem;
    std::map<block_based_elem_type, int> map_elem_idx;
    // size * size arrays, matrices
    std::vector<int> arr_order, arr_oplus, arr_inf, arr_sup;
    // size  array, vector
    std::vector<int> vec_orthosupplement;
    int size; // number of elements in lattice_ea

    const Lattice_EA lea;
};
