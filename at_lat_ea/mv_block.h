#pragma once

#include <iomanip>
#include <iostream>
#include <vector>

struct MV_Block {
    static bool MV_Block_sort_alg(const std::pair<int, int>& p1, const std::pair<int, int>& p2);
    using elem_type = std::vector<int>;
    using annotated_type = std::vector<std::pair<int, int>>;

    int size;
    int max_id;

    std::vector<int> ids;
    std::vector<int> orders;

    MV_Block(int size, const std::vector<int>& ids, const std::vector<int>& orders);

    int get_atom_order_by_id(int id) const;
    int get_atom_index_by_id(int id) const;

    int get_max_id() const noexcept;
    bool is_elem(const MV_Block::elem_type& mv_ea_elem) const;

    bool is_less_or_equal(const MV_Block::elem_type& mv_ea_elem1, const MV_Block::elem_type& mv_ea_elem2) const;
    bool are_leq(int id1, int id2) const;

    int get_id(const MV_Block::elem_type& mv_ea_elem) const;
    // possibly unsorted, or shorter
    int get_id(const MV_Block::annotated_type& mv_ea_elem) const;

    MV_Block::elem_type get_elem_from_id(int id) const;
    MV_Block::annotated_type get_annotated_elem_from_id(int id) const;

    void disp_elem(const MV_Block::elem_type& mv_ea_elem, std::ostream& os = std::cout) const;
    void disp_elem_with_ids(const MV_Block::elem_type& mv_ea_elem) const;
    void disp_elem_annotated(const MV_Block::annotated_type& mv_ea_elem) const;

    MV_Block::elem_type orthosupplement(const MV_Block::elem_type& mv_ea_elem) const;
    int orthosupplement(int id) const;

    MV_Block::elem_type oplus(const MV_Block::elem_type& mv_ea_elem1, const MV_Block::elem_type& mv_ea_elem2) const;
	private:

    int _get_max_id() const;
};
