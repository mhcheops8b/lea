#include "lattice_ea.h"
#include <algorithm>
#include <iostream>

/* std::ostream& operator<<(std::ostream& os, const block_based_elem_type& bel)
{
    os << '{' << bel.block_id << ", " << bel.elem_id << '}';
    return os;
} */

std::string block_based_elem_type::toString()
{
    int log_length = 0;
    int max = std::max(block_id, elem_id);
    while (max > 0) {
        max /= 10;
        log_length++;
    }

    char aa[10] = "% d";
    aa[1] = (char)(48 + log_length);
    //std::cout << std::string(aa) << '\n';

    char buf1[128], buf2[128];
    std::sprintf(buf1, aa, block_id);
    std::sprintf(buf2, aa, elem_id);

    return "{ " + std::string(buf1) + ", " + std::string(buf2) + " }";
}

bool operator<(const block_based_elem_type& bel1, const block_based_elem_type& bel2)
{
    return (bel1.block_id < bel2.block_id || (bel1.block_id == bel2.block_id && bel1.elem_id < bel2.elem_id));
}

std::vector<int> Lattice_EA::get_block_ids(int atom_id) const
{
    std::vector<int> v;

    int id = 0;
    for (const auto& block : blocks) {
        bool found = false;
        for (const auto& elem_id : block.ids)
            if (elem_id == atom_id) {
                found = true;
                break;
            }

        if (found)
            v.push_back(id);

        id++;
    }

    return v;
}

std::vector<int> Lattice_EA::get_block(const std::vector<int>& atom_ids) const
{
    std::vector<int> v;
    if (atom_ids.size() < 1)
        return v;

    v = get_block_ids(atom_ids[0]);

    for (size_t i = 1; i < atom_ids.size(); i++) {
        std::vector<int> v2 = get_block_ids(atom_ids[i]);
        std::vector<int> v_intersection(blocks.size());

        std::vector<int>::iterator it = std::set_intersection(v.begin(), v.end(), v2.begin(), v2.end(), v_intersection.begin());
        v_intersection.resize(it - v_intersection.begin());

        v = v_intersection;
    }

    return v;
}

bool Lattice_EA::is_elem(const std::vector<std::pair<int, int>>& elem) const
{
    if (elem.size() < 1)
        return false;

    std::vector<int> ids(elem.size());
    for (size_t i = 0; i < elem.size(); i++)
        ids[i] = elem[i].first;

    std::vector<int> block_ids = get_block(ids);

    int size = block_ids.size();

    if (size == 0) {
        std::cerr << "No block containing all atoms in elem" << '\n';
        return false;
    }

    if (size >= 2) {
        if (elem.size() != 1) {
            std::cerr << "More than one block shares the more atoms" << '\n';
            return false;
        } else {
            MV_Block b{ blocks[block_ids[0]] };

            int elem_id = b.get_atom_order_by_id(elem[0].first);

            if (elem[elem_id].second < 0 || elem[elem_id].second > b.orders[elem_id])
                return false;

            return true;
        }
    }

    if (size == 1) {
        MV_Block b{ blocks[block_ids[0]] };

        for (size_t i = 0; i < elem.size(); i++) {
            int elem_id = b.get_atom_order_by_id(elem[i].first);

            if (elem[elem_id].second < 0 || elem[elem_id].second > b.orders[elem_id])
                return false;
        }
        return true;
    }

    return true;
}

int Lattice_EA::common_atom(int block1_id, int block2_id) const
{
    int max_size = std::max(blocks[block1_id].size, blocks[block2_id].size);
    std::vector<int> result(max_size);

    std::vector<int>::iterator it = std::set_intersection(
        blocks[block1_id].ids.begin(), blocks[block1_id].ids.end(),
        blocks[block2_id].ids.begin(), blocks[block2_id].ids.end(),
        result.begin());

    result.resize(it - result.begin());

    if (result.size() == 1)
        return result[0];
    else
        return 0;
}

bool Lattice_EA::are_same(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
{
    // is equal to 0
    if (bel1.elem_id == 0 && bel2.elem_id == 0)
        return true;
    //
    // is equal to 1
    if ((bel1.elem_id == blocks[bel1.block_id].get_max_id() - 1) && (bel2.elem_id == blocks[bel2.block_id].get_max_id() - 1))
        return true;

    // are from the same block
    if (bel1.block_id == bel2.block_id && bel1.elem_id == bel2.elem_id)
        return true;

    if (bel1.block_id != bel2.block_id) {
        int common_atom_id = common_atom(bel1.block_id, bel2.block_id);
        // are from different blocks without common atom
        if (common_atom_id == 0)
            return false;

        // are from different blocks with common atom
        int atom_order = blocks[bel1.block_id].get_atom_order_by_id(common_atom_id);
        for (int i = 1; i <= atom_order; i++) {
            // are equal to i.atom
            MV_Block::annotated_type el{ { common_atom_id, i } };
            if (bel1.elem_id == blocks[bel1.block_id].get_id(el) && bel2.elem_id == blocks[bel2.block_id].get_id(el))
                return true;
            // are equal to (i.atom)'
            MV_Block::annotated_type el_c1 = blocks[bel1.block_id].get_annotated_elem_from_id(blocks[bel1.block_id].get_max_id() - 1),
                                     el_c2 = blocks[bel2.block_id].get_annotated_elem_from_id(blocks[bel2.block_id].get_max_id() - 1);

            el_c1[blocks[bel1.block_id].get_atom_index_by_id(common_atom_id)].second = atom_order - i;
            el_c2[blocks[bel2.block_id].get_atom_index_by_id(common_atom_id)].second = atom_order - i;

            if (bel1.elem_id == blocks[bel1.block_id].get_id(el_c1) && bel2.elem_id == blocks[bel2.block_id].get_id(el_c2))
                return true;
        }
    }

    return false;
}

/*
 *
bool Lattice_EA::are_leq(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
 *
 * For each representant r1 of bel1
 *
 * 	For each representant r2 of bel2
 *
 * 		if r1 and r2 are in the same block -> answer
 * 		if r1 and r2 are int the blocks sharing an atom -> answer
 */

bool Lattice_EA::are_leq(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
{
    if (bel1.block_id == bel2.block_id || common_atom(bel1.block_id, bel2.block_id) != 0)
        return are_leq_proto(bel1, bel2);

    for (int b1 = 0; b1 < (int)blocks.size(); ++b1)
        for (int e1 = 0; e1 < blocks[b1].get_max_id(); ++e1)
            if (are_same({ b1, e1 }, bel1))
                for (int b2 = 0; b2 < (int)blocks.size(); ++b2)
                    for (int e2 = 0; e2 < blocks[b2].get_max_id(); ++e2)
                        if (are_same({ b2, e2 }, bel2)) {
                            if (b1 == b2 || common_atom(b1, b2) != 0)
                                return are_leq_proto({ b1, e1 }, { b2, e2 });
                        }

    return false;
}

bool Lattice_EA::are_leq_proto(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
{
    if (bel1.block_id == bel2.block_id) {
        return blocks[bel1.block_id].are_leq(bel1.elem_id, bel2.elem_id);
    } else {
        // first equals 0
        if (bel1.elem_id == 0)
            return true;
        // second equalt 1
        if (bel2.elem_id == blocks[bel2.block_id].get_max_id() - 1)
            return true;

        int common_atom_id = common_atom(bel1.block_id, bel2.block_id);
        if (common_atom_id == 0) {
            // no common atom
            return false;
        } else {
            // common atom
            // second represents a, 2a, .. na, (na)', ..., (2a)', (a)'
            // and first is less than or equal to its block1 representation
            int atom_order = blocks[bel1.block_id].get_atom_order_by_id(common_atom_id);
            for (int i = 1; i <= atom_order; i++) {
                // are equal to i.atom
                MV_Block::annotated_type el{ { common_atom_id, i } };

                if ((bel2.elem_id == blocks[bel2.block_id].get_id(el) && blocks[bel1.block_id].are_leq(bel1.elem_id, blocks[bel1.block_id].get_id(el))) || (bel1.elem_id == blocks[bel1.block_id].get_id(el) && blocks[bel2.block_id].are_leq(blocks[bel2.block_id].get_id(el), bel2.elem_id))) {
                    return true;
                }
                // are equal to (i.atom)'
                MV_Block::annotated_type el_c1 = blocks[bel1.block_id].get_annotated_elem_from_id(blocks[bel1.block_id].get_max_id() - 1),
                                         el_c2 = blocks[bel2.block_id].get_annotated_elem_from_id(blocks[bel2.block_id].get_max_id() - 1);

                el_c1[blocks[bel1.block_id].get_atom_index_by_id(common_atom_id)].second = atom_order - i;
                el_c2[blocks[bel2.block_id].get_atom_index_by_id(common_atom_id)].second = atom_order - i;

                if ((bel2.elem_id == blocks[bel2.block_id].get_id(el_c2) && blocks[bel1.block_id].are_leq(bel1.elem_id, blocks[bel1.block_id].get_id(el_c1))) || (bel1.elem_id == blocks[bel1.block_id].get_id(el_c1) && blocks[bel2.block_id].are_leq(blocks[bel2.block_id].get_id(el_c2), bel2.elem_id))) {
                    return true;
                }
            }
        }
    }
    return false;
}

int Lattice_EA::get_id_total(const block_based_elem_type& bel) const
{
    int total_id = 0;
    for (int i = 0; i < bel.block_id; i++)
        total_id += blocks[bel.block_id].get_max_id();

    total_id += bel.elem_id;
    return total_id;
}

block_based_elem_type Lattice_EA::get_cannonical(const block_based_elem_type& bel) const
{
    // is equal 0
    if (bel.elem_id == 0)
        return { 0, 0 };

    // is equal 1
    if (bel.elem_id == blocks[bel.block_id].get_max_id() - 1)
        return { (int)blocks.size() - 1, blocks.back().get_max_id() - 1 };

    for (int i = 0; i < bel.block_id; i++)
        for (int j = 0; j < blocks[i].get_max_id(); j++)
            if (are_same({ i, j }, bel))
                return { i, j };

    return bel;
}

block_based_elem_type Lattice_EA::oplus(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
{
    if (bel1.elem_id == 0)
        return get_cannonical(bel2);

    if (bel2.elem_id == 0)
        return get_cannonical(bel1);

    if (bel1.block_id == bel2.block_id) {
        MV_Block b = blocks[bel1.block_id];
        MV_Block::elem_type oplus = b.oplus(b.get_elem_from_id(bel1.elem_id), b.get_elem_from_id(bel2.elem_id));
        if (oplus[0] == -1) {
            // undefined
            return { -1, -1 };
        } else {
            return get_cannonical({ bel1.block_id, b.get_id(oplus) });
        }

    } else {

        block_based_elem_type b1 = get_cannonical(bel1),
                              b2 = get_cannonical(bel2);

        if (b1.block_id == b2.block_id) {
            return oplus(b1, b2);
        } else {
            //block_based_elem_type b2_c = {b2.block_id, blocks[b2.block_id].orthosupplement(b2.elem_id)};
            auto b2_c = orthosupplement(b2);
            /*
			std::cout << "\n  b1 " << b1;
			std::cout << "\n  b2_c " << b2_c;
			std::cout << '\n';
*/

            if (are_leq(b1, b2_c)) {
                block_based_elem_type b2_in_b1_block = get_block_representation(b2, b1.block_id);

                block_based_elem_type b1_in_b2_block = get_block_representation(b1, b2.block_id);

                /*
				std::cout << "b2_in_b1_block: " << b2_in_b1_block << '\n';
				std::cout << "b1_in_b2_block: " << b1_in_b2_block << '\n';
				*/

                if (b2_in_b1_block.block_id != -1) {
                    MV_Block b{ blocks[b1.block_id] };
                    return get_cannonical({ b1.block_id,
                        b.get_id(b.oplus(b.get_elem_from_id(b1.elem_id), b.get_elem_from_id(b2_in_b1_block.elem_id))) });
                }

                if (b1_in_b2_block.block_id != -1) {
                    MV_Block b{ blocks[b2.block_id] };
                    return get_cannonical({ b2.block_id,
                        b.get_id(b.oplus(b.get_elem_from_id(b1_in_b2_block.elem_id), b.get_elem_from_id(b2.elem_id))) });
                }
            }
        }
    }

    return { -1, -1 };
}

block_based_elem_type Lattice_EA::get_block_representation(const block_based_elem_type& bel, int block_id) const
{
    // from the same block?
    if (block_id == bel.block_id)
        return bel;

    // is 0
    if (bel.elem_id == 0)
        return { block_id, 0 };

    // is 1
    if (bel.elem_id == blocks[bel.block_id].get_max_id() - 1)
        return { block_id, blocks[block_id].get_max_id() - 1 };

    int common_atom_id = common_atom(bel.block_id, block_id);
    if (common_atom_id == 0) {
        // no common atom
        return { -1, -1 };
    } else {
        MV_Block b{ blocks[block_id] };
        int atom_idx = b.get_atom_index_by_id(common_atom_id);
        std::vector<int> v(b.size);

        std::vector<int> v_c(b.size);
        for (int i = 0; i < b.size; i++)
            v_c[i] = b.orders[i];

        for (int i = 1; i <= b.orders[atom_idx]; i++) {
            // atoms
            v[atom_idx] = i;
            int elem_idx = b.get_id(v);
            if (are_same(bel, { block_id, elem_idx }))
                return { block_id, elem_idx };

            // coatoms
            v_c[atom_idx] = b.orders[atom_idx] - i;
            elem_idx = b.get_id(v_c);
            if (are_same(bel, { block_id, elem_idx }))
                return { block_id, elem_idx };
        }
    }
    return { -1, -1 };
}

int Lattice_EA::get_elem_id(const block_based_elem_type& bel) const
{
    return bel.elem_id;
}

MV_Block Lattice_EA::get_MV_block(const block_based_elem_type& bel) const
{
    return blocks[bel.block_id];
}

block_based_elem_type Lattice_EA::orthosupplement(const block_based_elem_type& bel) const
{
    return get_cannonical({ bel.block_id, get_MV_block(bel).orthosupplement(get_elem_id(bel)) });
}

block_based_elem_type Lattice_EA::inf(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
{
    std::vector<block_based_elem_type> common_lower_bounds;

    for (int i = 0; i < (int)blocks.size(); i++)
        for (int e = 0; e < blocks[i].get_max_id(); e++) {
            // for all cannonical elements
            auto p = get_cannonical({ i, e });
            if (p.block_id == i && p.elem_id == e) {
                if (are_leq(p, bel1) && are_leq(p, bel2))
                    common_lower_bounds.push_back(p);
            }
        }

    bool is_first = true;
    block_based_elem_type sup{ -1, -1 };
    for (const auto& p : common_lower_bounds) {
        if (is_first) {
            sup = p;
            is_first = false;
        } else {
            if (are_leq(sup, p))
                sup = p;
        }
    }

    // test whether it it is greatest element
    bool ok = true;
    for (const auto& p : common_lower_bounds)
        if (!are_leq(p, sup)) {
            ok = false;
            break;
        }

    if (ok)
        return sup;
    else
        return { -1, -1 };
}

block_based_elem_type Lattice_EA::sup(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
{
    std::vector<block_based_elem_type> common_upper_bounds;

    for (int i = 0; i < (int)blocks.size(); i++)
        for (int e = 0; e < blocks[i].get_max_id(); e++) {
            // for all cannonical elements
            auto p = get_cannonical({ i, e });
            if (p.block_id == i && p.elem_id == e) {
                if (are_leq(bel1, p) && are_leq(bel2, p))
                    common_upper_bounds.push_back(p);
            }
        }

    bool is_first = true;
    block_based_elem_type inf{ -1, -1 };
    for (const auto& p : common_upper_bounds) {
        if (is_first) {
            inf = p;
            is_first = false;
        } else {
            if (are_leq(p, inf))
                inf = p;
        }
    }

    // test whether it it is smallest element
    bool ok = true;
    for (const auto& p : common_upper_bounds)
        if (!are_leq(inf, p)) {
            ok = false;
            break;
        }

    if (ok)
        return inf;
    else
        return { -1, -1 };
}

bool Lattice_EA::is_cannonical(const block_based_elem_type& bel) const
{
    auto p = get_cannonical(bel);
    if (p.block_id == bel.block_id && p.elem_id == bel.elem_id)
        return true;
    else
        return false;
}

block_based_elem_type Lattice_EA::impl(const block_based_elem_type& bel1, const block_based_elem_type& bel2) const
{
    return oplus(orthosupplement(bel1), inf(bel1, bel2));
}
