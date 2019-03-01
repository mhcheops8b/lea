#pragma once

#include "mv_block.h"
#include <string>
#include <vector>

// vector of pairs (atom_id, atom_mult)
using lattice_ea_elem_type = std::vector<std::pair<int, int>>;

struct block_based_elem_type {
    int block_id; // index of the block in the blocks
    int elem_id; // id of element in the current block

    std::string toString();

    friend bool operator<(const block_based_elem_type& bel1, const block_based_elem_type& bel2);
};

struct Lattice_EA {
    const std::vector<MV_Block> blocks;

    MV_Block get_MV_block(const block_based_elem_type& bel) const;
    int get_elem_id(const block_based_elem_type& bel) const;

    std::vector<int> get_block_ids(int atom_id) const;

    std::vector<int> get_block(const std::vector<int>& atom_ids) const;

    bool is_elem(const lattice_ea_elem_type& elem) const;

    int common_atom(int block1_id, int block2_id) const;

    bool are_same(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const;
    bool are_leq_proto(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const;
    bool are_leq(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const;
    int get_id_total(const block_based_elem_type& bel) const;

    block_based_elem_type get_cannonical(const block_based_elem_type& bel) const;
    bool is_cannonical(const block_based_elem_type& bel) const;

    block_based_elem_type oplus(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const;

    block_based_elem_type orthosupplement(const block_based_elem_type& bel) const;

    block_based_elem_type inf(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const;
    block_based_elem_type sup(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const;

    block_based_elem_type impl(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const;

    block_based_elem_type get_block_representation(const block_based_elem_type& bel, int block_id) const;
};
