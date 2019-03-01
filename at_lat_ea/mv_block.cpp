#include "algorithm"
#include "utility"

#include "mv_block.h"

bool MV_Block::MV_Block_sort_alg(const std::pair<int, int>& p1, const std::pair<int, int>& p2)
{
    return p1.first < p2.first;
}

MV_Block::MV_Block(int size, const std::vector<int>& ids, const std::vector<int>& orders)
{
    this->size = size;
    std::vector<std::pair<int, int>> tmp_(size);

    for (int i = 0; i < size; i++)
        tmp_[i] = std::pair<int, int>(ids[i], orders[i]);

    // sort items by atom_id
    std::sort(tmp_.begin(), tmp_.end(), MV_Block_sort_alg);

    // allocate memory 
    this->ids.resize(size);
    this->orders.resize(size);

    for (int i = 0; i < size; i++) {
        this->ids[i] = tmp_[i].first;
        this->orders[i] = tmp_[i].second;
    }

    // initialize max_id
    max_id = _get_max_id();
}

int MV_Block::get_atom_order_by_id(int id) const
{
    for (int i = 0; i < size; i++)
        if (ids[i] == id)
            return orders[i];

    return 0;
}

int MV_Block::get_atom_index_by_id(int id) const
{
    for (int i = 0; i < size; i++)
        if (ids[i] == id)
            return i;

    return -1;
}

int MV_Block::get_max_id() const noexcept 
{
	return max_id;
}

int MV_Block::_get_max_id() const
{
    int max_value_id = 1;

    for (const auto& v : orders)
        max_value_id *= (v + 1);

    return max_value_id;
}

/* Asserts:
 * 	size == mv_ea_elem.size()
 */
bool MV_Block::is_elem(const MV_Block::elem_type& mv_ea_elem) const
{
    /*
	if (size != (int)mv_ea_elem.size())
	{
		std::cerr << "Size of the element differs from " << size << '\n';
		return false;
	}
	*/

    for (int i = 0; i < size; i++)
        if (mv_ea_elem[i] < 0 || mv_ea_elem[i] > orders[i])
            return false;

    return true;
}

/* Asserts:
 * 	is_elem(mv_ea_elem1)
 * 	is_elem(mv_ea_elem2)
 */
bool MV_Block::is_less_or_equal(const MV_Block::elem_type& mv_ea_elem1, const MV_Block::elem_type& mv_ea_elem2) const
{
    for (int i = 0; i < size; i++)
        if (mv_ea_elem1[i] > mv_ea_elem2[i])
            return false;

    return true;
}

bool MV_Block::are_leq(int id1, int id2) const
{
    return is_less_or_equal(get_elem_from_id(id1), get_elem_from_id(id2));
}

int MV_Block::get_id(const MV_Block::elem_type& mv_ea_elem) const
{
    if (!is_elem(mv_ea_elem))
        return -1; // undefined

    int value = mv_ea_elem[size - 1];
    for (int i = size - 2; i >= 0; --i) {
        value = value * (orders[i] + 1) + mv_ea_elem[i];
    }

    return value;
}

// possibly unsorted, or shorter
int MV_Block::get_id(const MV_Block::annotated_type& mv_ea_elem) const
{
    MV_Block::elem_type v(size);

    for (const auto& p : mv_ea_elem) {
        for (int i = 0; i < size; i++)
            if (ids[i] == p.first)
                v[i] = p.second;
    }

    return get_id(v);
}

MV_Block::elem_type MV_Block::get_elem_from_id(int id) const
{

    int max_value_id = get_max_id();

    std::vector<int> result(size);

    if (id < 0 || id >= max_value_id) {
        for (int i = 0; i < size; i++)
            result[i] = -1; // undefined

        return result;
    }

    for (int i = 0; i < size; i++) {
        result[i] = id % (orders[i] + 1);

        id = (id - result[i]) / (orders[i] + 1);
    }

    return result;
}

MV_Block::annotated_type MV_Block::get_annotated_elem_from_id(int id) const
{

    int max_value_id = get_max_id();

    MV_Block::annotated_type result(size);

    if (id < 0 || id >= max_value_id) {
        for (int i = 0; i < size; i++)
            result[i] = std::pair<int, int>(ids[i], -1); // undefined

    } else {
        for (int i = 0; i < size; i++) {
            int val = id % (orders[i] + 1);
            result[i] = std::pair<int, int>(ids[i], val);

            id = (id - val) / (orders[i] + 1);
        }
    }
    // sort by ids
    std::sort(result.begin(), result.end(), MV_Block_sort_alg);

    return result;
}

void MV_Block::disp_elem(const MV_Block::elem_type& mv_ea_elem, std::ostream& os) const
{
    if (!is_elem(mv_ea_elem)) {
        os << "-1" << '\n';
    } else {
        os << '{';
        for (int i = 0; i < size; i++) {
            if (i)
                os << ',';
            os << mv_ea_elem[i];
        }
        os << '}';
    }
}

void MV_Block::disp_elem_with_ids(const MV_Block::elem_type& mv_ea_elem) const
{
    if (!is_elem(mv_ea_elem)) {
        std::cout << '{';
        for (int i = 0; i < size; i++) {
            if (i)
                std::cout << ',';
            std::cout << '{' << ids[i] << ',' << -1 << '}';
        }
        std::cout << '}';
    } else {
        std::cout << '{';
        for (int i = 0; i < size; i++) {
            if (i)
                std::cout << ',';
            std::cout << '{' << ids[i] << ',' << mv_ea_elem[i] << '}';
        }
        std::cout << '}';
    }
}

void MV_Block::disp_elem_annotated(const MV_Block::annotated_type& mv_ea_elem) const
{
    bool is_first = true;
    for (const auto& p : mv_ea_elem) {
        if (is_first)
            is_first = false;
        else
            std::cout << ',';
        std::cout << '[' << p.first << ", " << p.second << ']';
    }
    std::cout << '\n';
}

MV_Block::elem_type MV_Block::orthosupplement(const MV_Block::elem_type& mv_ea_elem) const
{
    MV_Block::elem_type e{ mv_ea_elem };

    for (int i = 0; i < size; i++)
        e[i] = orders[i] - e[i];

    return e;
}

int MV_Block::orthosupplement(int id) const
{
    return get_id(orthosupplement(get_elem_from_id(id)));
}

MV_Block::elem_type MV_Block::oplus(const MV_Block::elem_type& mv_ea_elem1, const MV_Block::elem_type& mv_ea_elem2) const
{
    MV_Block::elem_type result(size);

    bool ok = true;
    for (int i = 0; i < size; i++) {
        int val = mv_ea_elem1[i] + mv_ea_elem2[i];
        if (val > orders[i]) {
            ok = false;
            break;
        } else {
            result[i] = val;
        }
    }

    if (ok)
        return result;
    else {
	    std::fill(result.begin(), result.end(), -1);
//        for (int i = 0; i < size; i++)
//            result[i] = -1;

        return result;
    }
}
