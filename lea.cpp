#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

template <size_t N>
struct atomic_MV_effect_def {
    std::array<int, N> atom_ids;
    std::array<int, N> atom_orders;
};

template <size_t N>
using MV_EA_elem = std::array<int, N>;

template <size_t N>
struct atomic_MV_effect_algebra {
    const atomic_MV_effect_def<N> atoms;

    //	atomic_MV_effect_algebra(const std::array<std::pair<int, int>, N> &pairs);

    constexpr int get_max_id() const;

    bool is_elem(const MV_EA_elem<N>& mv_ea_elem) const;

    bool is_less_or_equal(const MV_EA_elem<N>& mv_ea_elem1, const MV_EA_elem<N>& mv_ea_elem2) const;

    bool are_summable(const MV_EA_elem<N>& mv_ea_elem1, const MV_EA_elem<N>& mv_ea_elem2) const;

    bool are_summable(const int id1, const int id2) const;

    MV_EA_elem<N> oplus(const MV_EA_elem<N>& el1, const MV_EA_elem<N>& el2) const;

    int oplus(int id1, int id2) const;

    int get_id(const MV_EA_elem<N>& mv_ea_elem) const;

    MV_EA_elem<N> get_elem_from_id(int id) const;

    void get_elem_from_id(int id, MV_EA_elem<N>& mv_ea_elem);

    void disp_elem(const MV_EA_elem<N>& mv_ea_elem) const;
    void disp_elem_with_ids(const MV_EA_elem<N>& mv_ea_elem) const;

    //std::array<std::array<int, this.get_max_id()>, this.get_max_id()> generate_oplus();
    void generate_oplus();
};

template <size_t N>
void atomic_MV_effect_algebra<N>::generate_oplus()
{
    size_t size = get_max_id();

    std::cout << '{';
    for (size_t i = 0; i < size; i++) {
        if (i)
            std::cout << ",\n";
        else
            std::cout << "\n";

        std::cout << "\t{";
        for (size_t j = 0; j < size; j++) {
            if (j)
                std::cout << ", ";
            std::cout << this->oplus(i, j);
        }
        std::cout << '}';
        //std::cout << '\n';
    }
    std::cout << "\n}\n";
}

#if 0
template <size_t N>
atomic_MV_effect_algebra<N>::atomic_MV_effect_algebra(const std::array<std::pair<int, int>, N> &pairs)
{
	int id = 0;
	for (const auto &p: pairs)
	{
		atoms.atom_ids[id] = p.first;
		atoms.atom_orders[id] = p.second;
		id++;
	}
}
#endif

template <size_t N>
constexpr int atomic_MV_effect_algebra<N>::get_max_id() const
{
    int max_value_id = 1;

    for (size_t i = 0; i < N; i++)
        max_value_id *= (atoms.atom_orders[i] + 1);

    return max_value_id;
}

template <size_t N>
bool atomic_MV_effect_algebra<N>::is_elem(const MV_EA_elem<N>& mv_ea_elem) const
{
    for (size_t i = 0; i < N; i++)
        if (mv_ea_elem[i] < 0 || mv_ea_elem[i] > atoms.atom_orders[i])
            return false;

    return true;
}

template <size_t N>
bool atomic_MV_effect_algebra<N>::is_less_or_equal(const MV_EA_elem<N>& mv_ea_elem1, const MV_EA_elem<N>& mv_ea_elem2) const
{
    for (size_t i = 0; i < N; i++)
        if (mv_ea_elem1[i] > mv_ea_elem2[i])
            return false;

    return true;
}

template <size_t N>
bool atomic_MV_effect_algebra<N>::are_summable(const MV_EA_elem<N>& mv_ea_elem1, const MV_EA_elem<N>& mv_ea_elem2) const
{
    if (!is_elem(mv_ea_elem1))
        return false;

    if (!is_elem(mv_ea_elem2))
        return false;

    for (size_t i = 0; i < N; i++)
        if (mv_ea_elem1[i] + mv_ea_elem2[i] > atoms.atom_orders[i])
            return false;

    return true;
}

template <size_t N>
bool atomic_MV_effect_algebra<N>::are_summable(int id1, int id2) const
{
    const int max_id = get_max_id();

    if (id1 < 0 || id1 > max_id)
        return false;

    if (id2 < 0 || id2 > max_id)
        return false;

    MV_EA_elem<N> el1 = get_elem_from_id(id1);
    MV_EA_elem<N> el2 = get_elem_from_id(id2);

    return are_summable(el1, el2);
}

template <size_t N>
MV_EA_elem<N> atomic_MV_effect_algebra<N>::oplus(const MV_EA_elem<N>& mv_ea_elem1, const MV_EA_elem<N>& mv_ea_elem2) const
{
    if (!is_elem(mv_ea_elem1))
        return { -1 };

    if (!is_elem(mv_ea_elem2))
        return { -1 };

    if (!are_summable(mv_ea_elem1, mv_ea_elem2))
        return { -1 };

    MV_EA_elem<N> result;
    for (size_t i = 0; i < N; i++)
        result[i] = mv_ea_elem1[i] + mv_ea_elem2[i];

    return result;
}

template <size_t N>
int atomic_MV_effect_algebra<N>::oplus(const int id1, const int id2) const
{
    const int max_id = get_max_id();

    if (id1 < 0 || id1 >= max_id)
        return -1;

    if (id2 < 0 || id2 >= max_id)
        return -1;

    if (!are_summable(id1, id2))
        return -1;

    return get_id(oplus(get_elem_from_id(id1), get_elem_from_id(id2)));
}

template <size_t N>
int atomic_MV_effect_algebra<N>::get_id(const MV_EA_elem<N>& mv_ea_elem) const
{
    if (!is_elem(mv_ea_elem))
        return -1; // undefined

    int value = mv_ea_elem[N - 1];
    for (int i = N - 2; i >= 0; --i) {
        value = value * (atoms.atom_orders[i] + 1) + mv_ea_elem[i];
    }

    return value;
}

template <size_t N>
MV_EA_elem<N> atomic_MV_effect_algebra<N>::get_elem_from_id(int id) const
{

    int max_value_id = get_max_id();

    if (id < 0 || id >= max_value_id)
        return { -1 }; // undefined

    MV_EA_elem<N> result;
    for (size_t i = 0; i < N; i++) {
        result[i] = id % (atoms.atom_orders[i] + 1);

        id = (id - result[i]) / (atoms.atom_orders[i] + 1);
    }

    return result;
}

template <size_t N>
void atomic_MV_effect_algebra<N>::disp_elem(const MV_EA_elem<N>& mv_ea_elem) const
{
    if (!is_elem(mv_ea_elem)) {
        std::cout << "-1" << '\n';
    } else {
        std::cout << '{';
        for (size_t i = 0; i < N; i++) {
            if (i)
                std::cout << ',';
            std::cout << mv_ea_elem[i];
        }
        std::cout << '}';
    }
}

template <size_t N>
void atomic_MV_effect_algebra<N>::disp_elem_with_ids(const MV_EA_elem<N>& mv_ea_elem) const
{
    if (!is_elem(mv_ea_elem)) {
        std::cout << '{';
        for (size_t i = 0; i < N; i++) {
            if (i)
                std::cout << ',';
            std::cout << '{' << atoms.atom_ids[i] << ',' << -1 << '}';
        }
        std::cout << '}';
    } else {
        std::cout << '{';
        for (size_t i = 0; i < N; i++) {
            if (i)
                std::cout << ',';
            std::cout << '{' << atoms.atom_ids[i] << ',' << mv_ea_elem[i] << '}';
        }
        std::cout << '}';
    }
}

template <size_t N>
bool is_associative(std::array<std::array<int, N>, N>& op, int op_undef = -1)
{
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            for (size_t k = 0; k < N; k++) {
                if (op[i][j] != op_undef && op[op[i][j]][k] != op_undef) {
                    if (op[j][k] == op_undef) {
                        std::cout << "->: ((j = " << j << ") + (k = " << k << ")) undefined.\n";
                        return false;
                    }

                    if (op[i][op[j][k]] == op_undef) {
                        std::cout << "->: i = " << i << " + (j = " << j << " + k = " << k << ") undefined.\n";

                        return false;
                    }

                    if (op[op[i][j]][k] != op[i][op[j][k]]) {
                        std::cout << "->: (i = " << i << " + j = " << j << ") + k = " << k << " != i + (j + k).\n";

                        return false;
                    }
                }

                if (op[j][k] != op_undef && op[i][op[j][k]] != op_undef) {
                    if (op[i][j] == op_undef) {
                        std::cout << "<-: i = " << i << " + j = " << j << " undefined.\n";
                        return false;
                    }

                    if (op[op[i][j]][k] == op_undef) {
                        std::cout << "<-: (i = " << i << " + j = " << j << ") + k = " << k << " undefined.\n";
                        return false;
                    }

                    if (op[op[i][j]][k] != op[i][op[j][k]]) {
                        std::cout << "->: (i = " << i << " + j = " << j << ") + k = " << k << " != i + (j + k).\n";

                        return false;
                    }
                }
            }

    return true;
}

template <size_t N>
bool get_inf(const std::array<std::array<int, N>, N>& order, const int elem_a, const int elem_b, int& inf)
{
    int candidates[N], inf_cand = 0;
    for (size_t i = 0; i < N; i++)
        candidates[i] = 0;

    for (size_t i = 0; i < N; i++) {
        if (order[i][elem_a] && order[i][elem_b]) {
            candidates[i] = 1;

            if (order[inf_cand][i])
                inf_cand = i;
        }
    }

#if 0
	for (size_t i = 0; i < N; i++) 
		if (candidates[i])
		       std::cout << i << ", ";
	std::cout << '\n';
#endif

    for (size_t i = 0; i < N; i++) {
        if (candidates[i] && !order[i][inf_cand]) {
            return false;
        }
    }

    inf = inf_cand;
    return true;
}

template <size_t N>
bool get_sup(const std::array<std::array<int, N>, N>& order, const int elem_a, const int elem_b, int& sup)
{
    int candidates[N], sup_cand = N - 1;
    for (size_t i = 0; i < N; i++)
        candidates[i] = 0;

    for (size_t i = 0; i < N; i++) {
        if (order[elem_a][i] && order[elem_b][i]) {
            candidates[i] = 1;

            if (order[i][sup_cand])
                sup_cand = i;
        }
    }

#if 0
	for (size_t i = 0; i < N; i++) 
		if (candidates[i])
		       std::cout << i << ", ";
	std::cout << '\n';
#endif

    for (size_t i = 0; i < N; i++) {
        if (candidates[i] && !order[sup_cand][i]) {
            return false;
        }
    }

    sup = sup_cand;
    return true;
}

template <size_t N>
bool is_compatible(const std::array<std::array<int, N>, N>& order, const std::array<std::array<int, N>, N>& oplus, const int elem_a, const int elem_b)
{
    int inf{ -1 }, sup{ -1 };

    get_inf(order, elem_a, elem_b, inf);
    get_sup(order, elem_a, elem_b, sup);

    if (inf == -1 || sup == -1)
        return false;

    int minus{ -1 };
    for (size_t i = 0; i < N; i++)
        if (oplus[i][inf] == elem_b)
            minus = i;

    if (minus == -1)
        return false;

    if (sup != oplus[elem_a][minus])
        return false;

    return true;
}

template <size_t N>
bool not_empty_set(const std::array<int, N>& subset)
{
    for (size_t i = 0; i < N; i++)
        if (subset[i])
            return true;

    return false;
}

template <size_t N>
bool mutually_compatible(const std::array<std::array<int, N>, N>& rel_compatibility, const std::array<int, N>& subset)
{
    for (size_t i = 0; i < N; i++)
        if (subset[i]) {
            for (size_t j = 0; j < N; j++)
                if (subset[j]) {
                    if (!rel_compatibility[i][j])
                        return false;
                }
        }
    return true;
}

template <size_t N>
bool compatible_with_elem(const std::array<std::array<int, N>, N>& rel_compatibility, const std::array<int, N>& subset, const int elem)
{
    for (size_t i = 0; i < N; i++)
        if (subset[i]) {
            if (!rel_compatibility[i][elem])
                return false;
        }

    return true;
}

template <size_t N>
void get_blocks_(int level, const std::array<std::array<int, N>, N>& rel_compatibility, std::array<int, N>& subset)
{
    if (level == N - 1) {
        if (not_empty_set(subset) && mutually_compatible(rel_compatibility, subset)) {
            bool ok = true;
            for (size_t i = 1; i < N - 1; i++) {
                if (!subset[i] && compatible_with_elem(rel_compatibility, subset, i)) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                for (size_t i = 0; i < N; i++)
                    if (subset[i])
                        std::cout << i << ", ";

                std::cout << '\n';
            }
        }
    } else {
        for (int i = 0; i <= 1; i++) {
            subset[level] = i;
            get_blocks_<N>(level + 1, rel_compatibility, subset);
            subset[level] = 0;
        }
    }
}

template <size_t N>
void get_blocks(const std::array<std::array<int, N>, N>& rel_compatibility)
{
    std::array<int, N> subset;
    subset[0] = 1;
    subset[N - 1] = 1;

    get_blocks_(1, rel_compatibility, subset);
}

template <size_t N>
struct lattice_ea {

    std::array<std::array<int, N>, N> order;
    std::array<std::array<int, N>, N> oplus;
    std::array<int, N> orthosupplement;
    std::array<std::array<int, N>, N> inf;
    std::array<std::array<int, N>, N> sup;

    int impl(const int a, const int b) const;

    template <size_t M>
    friend bool is_filter(const lattice_ea<M>& lea, const std::array<int, M>& subset);
    template <size_t M>
    friend bool is_fantastic(const lattice_ea<M>& lea, const std::array<int, M>& subset);
    template <size_t M>
    friend bool is_implicative(const lattice_ea<M>& lea, const std::array<int, M>& subset);
    template <size_t M>
    friend bool is_positive_implicative(const lattice_ea<M>& lea, const std::array<int, M>& subset);
    template <size_t M>
    friend bool is_strong(const lattice_ea<M>& lea, const std::array<int, M>& subset);
};

template <size_t N>
int lattice_ea<N>::impl(const int a, const int b) const
{
    return oplus[orthosupplement[a]][inf[a][b]];
}

template <size_t N>
bool is_filter(const lattice_ea<N>& lea, const std::array<int, N>& subset)
{
    for (size_t i = 0; i < N; i++)
        if (subset[i])
            for (size_t j = 0; j < N; j++) {
                if (lea.impl(i, j) == -1)
                    return false;

                if (subset[lea.impl(i, j)] && !subset[j])
                    return false;
            }

    return true;
}

template <size_t N>
bool is_fantastic(const lattice_ea<N>& lea, const std::array<int, N>& subset)
{
    for (size_t c = 0; c < N; c++)
        if (subset[c]) {
            for (size_t a = 0; a < N; a++)
                for (size_t b = 0; b < N; b++) {
#if 0
					if (!is_compatible(lea.order, lea.oplus, a, b))
					{
					int val = lea.impl(lea.impl(lea.impl(a,b),b),a);
					if ( val != N - 1)
						std::cout << "FF: a = " << a << ", b = " << b << ", ((a->b)->b)->a = " << val << '\n';
					}
#endif
                    int val1 = lea.impl(b, a);
                    if (val1 == -1)
                        return false;

                    if (lea.impl(a, b) == -1 || lea.impl(lea.impl(a, b), b) == -1 || lea.impl(lea.impl(lea.impl(a, b), b), a) == -1)
                        return false;

                    if (subset[lea.impl(c, lea.impl(b, a))] && !subset[lea.impl(lea.impl(lea.impl(a, b), b), a)])
                        return false;
                }
        }

    return true;
}

template <size_t N>
bool is_implicative(const lattice_ea<N>& lea, const std::array<int, N>& subset)
{
    for (size_t a = 0; a < N; a++)
        for (size_t b = 0; b < N; b++) {
            if (lea.impl(a, b) == -1)
                return false;

            if (subset[lea.impl(a, b)]) {

                for (size_t c = 0; c < N; c++) {
                    if (lea.impl(b, c) == -1 || lea.impl(a, lea.impl(b, c)) == -1 || lea.impl(a, c) == -1)
                        return false;

                    if (subset[lea.impl(a, lea.impl(b, c))] && !subset[lea.impl(a, c)])
                        return false;
                }
            }
        }

    return true;
}

template <size_t N>
bool is_positive_implicative(const lattice_ea<N>& lea, const std::array<int, N>& subset)
{
    for (size_t a = 0; a < N; a++)
        if (subset[a]) {
            for (size_t b = 0; b < N; b++)
                for (size_t c = 0; c < N; c++)
                    if (subset[lea.impl(a, lea.impl(lea.impl(b, c), b))] && !subset[b])
                        return false;
        }

    return true;
}

template <size_t N>
bool is_strong(const lattice_ea<N>& lea, const std::array<int, N>& subset)
{
    if (!is_filter(lea, subset))
        return false;

    for (size_t a = 0; a < N; a++)
        if (subset[a]) {
            for (size_t b = 0; b < N; b++)
                if (!subset[lea.impl(b, a)])
                    return false;
        }

    return true;
}

template <size_t N>
void get_filters_(int level, const lattice_ea<N>& lea, std::array<int, N>& subset, bool (&f)(const lattice_ea<N>& lea, const std::array<int, N>& subset))
{
    if (level == N - 1) {
        if (f(lea, subset)) {

            for (size_t i = 0; i < N; i++)
                if (subset[i])
                    std::cout << i << ", ";

            std::cout << '\n';
        }

#if 0
			if (not_empty_set(subset) && mutually_compatible(rel_compatibility, subset))
			{
				bool ok = true;
				for (size_t i = 1; i < N - 1; i++)
				{
					if (!subset[i] && compatible_with_elem(rel_compatibility, subset, i))
					{
						ok = false;
						break;
					}
				}

				if (ok)
				{
					for (size_t i = 0; i < N; i++)
						if (subset[i])
							std::cout << i << ", ";

					std::cout << '\n';
				}


			}
#endif
    } else {
        for (int i = 0; i <= 1; i++) {
            subset[level] = i;
            get_filters_<N>(level + 1, lea, subset, f);
            subset[level] = 0;
        }
    }
}

template <size_t N>
void get_filters(const lattice_ea<N>& lea, bool (&f)(const lattice_ea<N>&, const std::array<int, N>&))
{
    std::array<int, N> subset;
    subset[N - 1] = 1;

    get_filters_(0, lea, subset, f);
}

template <size_t N>
void compute_induced_order(const std::array<std::array<int, N>, N>& oplus, std::array<std::array<int, N>, N>& order)
{
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++) {
            bool found = false;
            if (i != j) {
                for (size_t k = 0; k < N; k++)
                    if (static_cast<int>(j) == oplus[i][k]) {
                        found = true;
                        break;
                    }

            } else
                found = true;

            if (found)
                order[i][j] = 1;
            else
                order[i][j] = 0;
        }
}

template <size_t N>
void compute_orthosupplement(const std::array<std::array<int, N>, N>& oplus, std::array<int, N>& orthosupplement)
{
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            if (oplus[i][j] == N - 1)
                orthosupplement[i] = j;
}

template <size_t N>
void print_binary(const std::array<std::array<int, N>, N>& op)
{
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++)
            std::cout << std::setw(3) << op[i][j] << ", ";

        std::cout << '\n';
    }
}

template <size_t N>
void print_unary(const std::array<int, N>& op)
{
    for (size_t i = 0; i < N; i++) {
        std::cout << std::setw(2) << op[i] << ", ";
    }

    std::cout << '\n';
}

template <size_t N>
void compute_rel_compatibility(const std::array<std::array<int, N>, N>& order, const std::array<std::array<int, N>, N>& oplus,
    std::array<std::array<int, N>, N>& rel_compatibility)
{
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (is_compatible(order, oplus, i, j)) {
                rel_compatibility[i][j] = 1;
            } else
                rel_compatibility[i][j] = 0;
        }
    }
}

template <size_t N>
bool compute_inf(const std::array<std::array<int, N>, N>& order, std::array<std::array<int, N>, N>& inf)
{
    bool ok = true;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            int v{ -1 };
            if (get_inf(order, i, j, v)) {
                inf[i][j] = v;
            } else {
                inf[i][j] = -1;
                ok = false;
            }
        }
    }
    return ok;
}

template <size_t N>
bool compute_sup(const std::array<std::array<int, N>, N>& order, std::array<std::array<int, N>, N>& sup)
{
    bool ok = true;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            int v{ -1 };
            if (get_sup(order, i, j, v)) {
                sup[i][j] = v;
            } else {
                sup[i][j] = -1;
                ok = false;
            }
        }
    }
    return ok;
}

template <size_t N>
void print_impl(const lattice_ea<N>& lea)
{
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            std::cout << std::setw(3) << lea.impl(i, j) << ", ";
        }
        std::cout << '\n';
    }
}

template <size_t N>
void do_lea2(const std::string& msg, const std::array<std::array<int, N>, N>& oplus)
{
    std::array<std::array<int, N>, N> order;
    compute_induced_order(oplus, order);

    std::array<int, N> orthosupplement;
    compute_orthosupplement(oplus, orthosupplement);

    std::array<std::array<int, N>, N> inf, sup;

    if (!compute_inf(order, inf)) {
        return;
    }

    if (!compute_sup(order, sup)) {
        return;
    }

    lattice_ea<N> lea{ order, oplus, orthosupplement, inf, sup };

    std::cout << msg; //"Example 3.12:\n";

    // test for ((a->b)->b)->a vs. b
    for (size_t a = 0; a < N; a++)
        for (size_t b = 0; b < N; b++) {
            if (lea.impl(a, b) == -1 || lea.impl(lea.impl(a, b), b) == -1 || lea.impl(lea.impl(lea.impl(a, b), b), a) == -1)
                break;

            int val{ lea.impl(lea.impl(lea.impl(a, b), b), a) };
            int val2{ lea.impl(b, a) };
            if (!order[val2][val])
                //			if (val != val2)
                std::cout << "b->a !<=! ((a->b)->b)->a: a = " << a << ", b = " << b << ", b->a = " << val2 << ", ((a->b)->b)->a = " << val << '\n';
        }
}

template <size_t N>
void do_lea(const std::string& msg, const std::array<std::array<int, N>, N>& oplus, bool do_all = true)
{
    std::cout << msg; //"Example 3.12:\n";

    std::array<std::array<int, N>, N> order;
    compute_induced_order(oplus, order);
    std::cout << " Order: \n";
    print_binary(order);

    std::array<int, N> orthosupplement;
    compute_orthosupplement(oplus, orthosupplement);
    std::cout << "Orthosupplement: ";
    print_unary(orthosupplement);

    std::array<std::array<int, N>, N> inf, sup;

    if (!compute_inf(order, inf)) {
        std::cout << "No lattice - inf!\n";
        if (!do_all)
            return;
    }

    std::cout << "Inf: \n";
    print_binary(inf);

    std::cout << "--\n";
    if (!compute_sup(order, sup)) {
        std::cout << "No lattice - sup!\n";
        if (!do_all)
            return;
    }

    std::cout << "Blocks: \n";
    std::array<std::array<int, N>, N> compatibility;
    compute_rel_compatibility(order, oplus, compatibility);
    get_blocks(compatibility);

    std::cout << "Sup: \n";
    print_binary(sup);

    lattice_ea<N> lea{ order, oplus, orthosupplement, inf, sup };

    std::cout << "Filters: \n";
    get_filters(lea, is_filter);

    std::cout << "Fantastic filters: \n";
    get_filters(lea, is_fantastic);

    std::cout << "Implicative filters: \n";
    get_filters(lea, is_implicative);

    std::cout << "Positive implicative filters: \n";
    get_filters(lea, is_positive_implicative);

    std::cout << "Strong filters: \n";
    get_filters(lea, is_strong);
    std::cout << "--\n";
    print_impl(lea);

    // test for (a->b)->b vs. b
    for (size_t a = 0; a < N; a++)
        for (size_t b = 0; b < N; b++) {
            if (lea.impl(a, b) == -1 || lea.impl(lea.impl(a, b), b) == -1)
                break;

            int val{ lea.impl(lea.impl(a, b), b) };
            if (!order[b][val])
                std::cout << "b !<=! (a->b)->b: a = " << a << ", b = " << b << ", (a->b)->b = " << val << '\n';
        }

    // test for ((a->b)->b)->a vs. b
    for (size_t a = 0; a < N; a++)
        for (size_t b = 0; b < N; b++) {
            if (lea.impl(a, b) == -1 || lea.impl(lea.impl(a, b), b) == -1 || lea.impl(lea.impl(lea.impl(a, b), b), a) == -1)
                break;

            int val{ lea.impl(lea.impl(lea.impl(a, b), b), a) };
            int val2{ lea.impl(b, a) };
            //			if (!order[val2][val])
            if (val != val2)
                std::cout << "b->a !<=! ((a->b)->b)->a: a = " << a << ", b = " << b << ", b->a = " << val2 << ", ((a->b)->b)->a = " << val << '\n';
        }
}

int main()
{
    const size_t ex_3_12_size = 10;
    std::array<std::array<int, ex_3_12_size>, ex_3_12_size> ex_3_12_oplus{ { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
        { 1, 3, 4, -1, 9, -1, -1, -1, -1, -1 },
        { 2, 4, -1, 9, -1, 6, -1, 8, -1, -1 },
        { 3, -1, 9, -1, -1, -1, -1, -1, -1, -1 },
        { 4, 9, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 5, -1, 6, -1, -1, -1, -1, 3, 9, -1 },
        { 6, -1, -1, -1, -1, -1, -1, 9, -1, -1 },
        { 7, -1, 8, -1, -1, 3, 9, -1, -1, -1 },
        { 8, -1, -1, -1, -1, 9, -1, -1, -1, -1 },
        { 9, -1, -1, -1, -1, -1, -1, -1, -1, -1 } } };

    std::cout << "Example 3.12:\n";

    std::array<std::array<int, ex_3_12_size>, ex_3_12_size> ex_3_12_order;
    compute_induced_order(ex_3_12_oplus, ex_3_12_order);
    std::cout << " Order: \n";
    print_binary(ex_3_12_order);

    std::array<int, ex_3_12_size> ex_3_12_orthosupplement;
    compute_orthosupplement(ex_3_12_oplus, ex_3_12_orthosupplement);
    std::cout << "Orthosupplement: ";
    print_unary(ex_3_12_orthosupplement);

    std::cout << "Blocks: \n";
    std::array<std::array<int, ex_3_12_size>, ex_3_12_size> ex_3_12_compatibility;
    compute_rel_compatibility(ex_3_12_order, ex_3_12_oplus, ex_3_12_compatibility);
    get_blocks(ex_3_12_compatibility);

    std::array<std::array<int, ex_3_12_size>, ex_3_12_size> ex_3_12_inf, ex_3_12_sup;

    if (!compute_inf(ex_3_12_order, ex_3_12_inf))
        std::cout << "No lattice - inf!\n";

    std::cout << "Inf: \n";
    print_binary(ex_3_12_inf);

    std::cout << "--\n";
    if (!compute_sup(ex_3_12_order, ex_3_12_sup))
        std::cout << "No lattice - sup!\n";

    std::cout << "Sup: \n";
    print_binary(ex_3_12_sup);

    lattice_ea<ex_3_12_size> lea_3_12{ ex_3_12_order, ex_3_12_oplus, ex_3_12_orthosupplement, ex_3_12_inf, ex_3_12_sup };

    std::cout << "Filters: \n";
    get_filters(lea_3_12, is_filter);

    std::cout << "Fantastic filters: \n";
    get_filters(lea_3_12, is_fantastic);

    std::cout << "Implicative filters: \n";
    get_filters(lea_3_12, is_implicative);

    std::cout << "Positive implicative filters: \n";
    get_filters(lea_3_12, is_positive_implicative);

    std::cout << "Strong filters: \n";
    get_filters(lea_3_12, is_strong);
    std::cout << "--\n";
    print_impl(lea_3_12);

    std::array<std::array<int, 14>, 14> ex_3_13_oplus{ { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 },
        { 1, 2, -1, 4, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 2, -1, -1, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 3, 4, 13, -1, -1, 6, -1, 8, -1, -1, -1, -1, -1, -1 },
        { 4, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 5, -1, -1, 6, -1, -1, -1, 2, 13, -1, -1, -1, -1, -1 },
        { 6, -1, -1, -1, -1, -1, -1, 13, -1, -1, -1, -1, -1, -1 },
        { 7, -1, -1, 8, -1, 2, 13, -1, -1, 10, -1, 12, -1, -1 },
        { 8, -1, -1, -1, -1, 13, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 9, -1, -1, -1, -1, -1, -1, 10, -1, -1, -1, 6, 13, -1 },
        { 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 13, -1, -1 },
        { 11, -1, -1, -1, -1, -1, -1, 12, -1, 6, 13, -1, -1, -1 },
        { 12, -1, -1, -1, -1, -1, -1, -1, -1, 13, -1, -1, -1, -1 },
        { 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } } };

    std::array<std::array<int, 14>, 14> test_ex_3_13_order{ { { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 } } };

    std::array<std::array<int, 14>, 14> ex_3_13_order;

    compute_induced_order(ex_3_13_oplus, ex_3_13_order);

    const int N = 14;
    /* compare test */
    std::cout << "Comparing orders: \n";
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            if (ex_3_13_order[i][j] != test_ex_3_13_order[i][j])
                std::cout << "i = " << i << ", j = " << j << '\n';

    std::cout << "Associativity 3.13: ";
    if (is_associative(ex_3_13_oplus, -1))
        std::cout << "OK";
    else
        std::cout << "Not ok.";
    std::cout << '\n';

    std::array<int, 14> test_ex_3_13_orthosupplement = {
        { 13, 4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9, 0 }
    };
    std::array<int, 14> ex_3_13_orthosupplement;

    compute_orthosupplement(ex_3_13_oplus, ex_3_13_orthosupplement);

    std::cout << "Comparing orthosupplements: \n";
    for (size_t i = 0; i < N; i++)
        if (ex_3_13_orthosupplement[i] != test_ex_3_13_orthosupplement[i])
            std::cout << i << '\n';

    std::array<std::array<int, 14>, 14> ex_3_13_inf, ex_3_13_sup;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int v;
            if (get_inf(ex_3_13_order, i, j, v)) {
                ex_3_13_inf[i][j] = v;
                std::cout << std::setw(3) << v << ", ";
            } else
                std::cout << std::setw(3) << '?' << ", ";
        }
        std::cout << '\n';
    }

    std::cout << "--\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int v;
            if (get_sup(ex_3_13_order, i, j, v)) {
                ex_3_13_sup[i][j] = v;
                std::cout << std::setw(3) << v << ", ";
            } else
                std::cout << std::setw(3) << '?' << ", ";
        }
        std::cout << '\n';
    }

    std::array<std::array<int, 14>, 14> ex_3_13_compatible;
    std::cout << "--\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (is_compatible(ex_3_13_order, ex_3_13_oplus, i, j)) {
                std::cout << std::setw(2) << 1 << ", ";
                ex_3_13_compatible[i][j] = 1;
            } else {
                std::cout << std::setw(2) << 0 << ", ";
                ex_3_13_compatible[i][j] = 0;
            }
        }
        std::cout << '\n';
    }

    get_blocks(ex_3_13_compatible);

    lattice_ea<14> lea_3_13{ ex_3_13_order, ex_3_13_oplus, ex_3_13_orthosupplement, ex_3_13_inf, ex_3_13_sup };

    std::cout << "Filters: \n";
    get_filters(lea_3_13, is_filter);

    std::cout << "Fantastic filters: \n";
    get_filters(lea_3_13, is_fantastic);

    std::cout << "Implicative filters: \n";
    get_filters(lea_3_13, is_implicative);

    std::cout << "Positive implicative filters: \n";
    get_filters(lea_3_13, is_positive_implicative);

    std::cout << "Strong filters: \n";
    get_filters(lea_3_13, is_strong);

    std::cout << "--\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << std::setw(3) << lea_3_13.impl(i, j) << ", ";
        }
        std::cout << '\n';
    }

    std::cout << "--\n";

    std::array<std::array<int, 18>, 18> ex_3_14_oplus{ {
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 },
        { 1, 2, -1, 4, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 2, -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 3, 4, 17, -1, -1, 6, -1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 4, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 5, -1, -1, 6, -1, -1, -1, 2, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 6, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 7, -1, -1, 8, -1, 2, 17, -1, -1, 10, -1, 12, -1, 15, 16, -1, -1, -1 },
        { 8, -1, -1, -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 9, -1, -1, -1, -1, -1, -1, 10, -1, -1, -1, 6, 17, -1, -1, -1, -1, -1 },
        { 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, -1, -1, -1 },
        { 11, -1, -1, -1, -1, -1, -1, 12, -1, 6, 17, -1, -1, 14, -1, 16, -1, -1 },
        { 12, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, 16, -1, -1, -1, -1 },
        { 13, -1, -1, -1, -1, -1, -1, 15, -1, -1, -1, 14, 16, 9, 6, 10, 17, -1 },
        { 14, -1, -1, -1, -1, -1, -1, 16, -1, -1, -1, -1, -1, 6, -1, 17, -1, -1 },
        { 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 16, -1, 10, 17, -1, -1, -1 },
        { 16, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, -1 },
        { 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    } };

    do_lea(std::string{ "Example 3.14:\n" }, ex_3_14_oplus);

    std::cout << "Associativity: ";
    if (is_associative(ex_3_14_oplus, -1))
        std::cout << "OK";
    else
        std::cout << "Not ok.";
    std::cout << '\n';

    std::array<std::array<int, 18>, 18> ex_3_14_order{ { { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 } } };

    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 18; j++) {
            int v;
            if (get_inf(ex_3_14_order, i, j, v)) {
                std::cout << std::setw(3) << v << ", ";
            } else
                std::cout << std::setw(3) << '?' << ", ";
        }
        std::cout << '\n';
    }

    std::cout << "--\n";
    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 18; j++) {
            int v;
            if (get_sup(ex_3_14_order, i, j, v)) {
                std::cout << std::setw(3) << v << ", ";
            } else
                std::cout << std::setw(3) << '?' << ", ";
        }
        std::cout << '\n';
    }

    std::cout << "--\n";
    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 18; j++) {
            if (is_compatible(ex_3_14_order, ex_3_14_oplus, i, j)) {
                std::cout << std::setw(2) << 1 << ", ";
            } else
                std::cout << std::setw(2) << 0 << ", ";
        }
        std::cout << '\n';
    }

#if 0
	int v;
	get_inf(ex_3_14_order, 6, 15, v);	
	std::cout << v << '\n';

#endif

    std::array<std::array<int, 18>, 18> ex_3_15_oplus{ {
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 },
        { 1, 2, -1, 4, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 2, -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 3, 4, 17, -1, -1, 6, -1, 8, -1, -1, -1, -1, -1, -1, 13, -1, 15, -1 },
        { 4, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 5, -1, -1, 6, -1, -1, -1, 2, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 6, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 7, -1, -1, 8, -1, 2, 17, 14, 13, 10, -1, 12, -1, -1, -1, 6, 5, -1 },
        { 8, -1, -1, -1, -1, 17, -1, 13, -1, -1, -1, -1, -1, -1, -1, -1, 6, -1 },
        { 9, -1, -1, -1, -1, -1, -1, 10, -1, -1, -1, 6, 17, -1, -1, -1, -1, -1 },
        { 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, -1, -1, -1 },
        { 11, -1, -1, -1, -1, -1, -1, 12, -1, 6, 17, -1, -1, -1, -1, -1, -1, -1 },
        { 12, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1 },
        { 14, -1, -1, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 17, 2, -1 },
        { 15, -1, -1, -1, -1, -1, -1, 6, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1 },
        { 16, -1, -1, 15, -1, -1, -1, 5, 6, -1, -1, -1, -1, 17, 2, -1, -1, -1 },
        { 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    } };

    std::cout << "$**$*$\n";
    const std::string msg{ "Example 3.15:\n" };
    do_lea(msg, ex_3_15_oplus);

    std::array<std::array<int, 14>, 14> ex_3_15a_oplus{ {
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 },
        { 1, 2, -1, 4, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 2, -1, -1, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 3, 4, 13, -1, -1, 6, -1, 8, -1, -1, 9, -1, 11, -1 },
        { 4, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 5, -1, -1, 6, -1, -1, -1, 2, 13, -1, -1, -1, -1, -1 },
        { 6, -1, -1, -1, -1, -1, -1, 13, -1, -1, -1, -1, -1, -1 },
        { 7, -1, -1, 8, -1, 2, 13, 10, 9, -1, -1, 6, 5, -1 },
        { 8, -1, -1, -1, -1, 13, -1, 9, -1, -1, -1, -1, 6, -1 },
        { 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 13, -1 },
        { 10, -1, -1, 9, -1, -1, -1, -1, -1, -1, -1, 13, 2, -1 },
        { 11, -1, -1, -1, -1, -1, -1, 6, -1, -1, 13, -1, -1, -1 },
        { 12, -1, -1, 11, -1, -1, -1, 5, 6, 13, 2, -1, -1, -1 },
        { 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    } };

    do_lea(std::string("Example 3.15a:\n"), ex_3_15a_oplus);

    std::array<std::array<int, 8>, 8> ex_chajda_oplus{ {
        { 0, 1, 2, 3, 4, 5, 6, 7 },
        { 1, 5, 4, -1, 7, -1, -1, -1 },
        { 2, 4, -1, 6, -1, 7, -1, -1 },
        { 3, -1, 6, 5, -1, -1, 7, -1 },
        { 4, 7, -1, -1, -1, -1, -1, -1 },
        { 5, -1, 7, -1, -1, -1, -1, -1 },
        { 6, -1, -1, 7, -1, -1, -1, -1 },
        { 7, -1, -1, -1, -1, -1, -1, -1 },

    } };

    do_lea(std::string("Example Chajda:\n"), ex_chajda_oplus);

    std::array<std::array<int, 9>, 9> ex_a2_b2{ { { 0, 1, 2, 3, 4, 5, 6, 7, 8 },
        { 1, 4, 3, 6, -1, 7, -1, 8, -1 },
        { 2, 3, 5, 7, 6, -1, 8, -1, -1 },
        { 3, 6, 7, 8, -1, -1, -1, -1, -1 },
        { 4, -1, 6, -1, -1, 8, -1, -1, -1 },
        { 5, 7, -1, -1, 8, -1, -1, -1, -1 },
        { 6, -1, 8, -1, -1, -1, -1, -1, -1 },
        { 7, 8, -1, -1, -1, -1, -1, -1, -1 },
        { 8, -1, -1, -1, -1, -1, -1, -1, -1 } } };

    do_lea(std::string("Example a2-b2:\n"), ex_a2_b2);

    std::array<std::array<int, 12>, 12> ex_a2_b2_b2_c2{ { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
        { 1, 4, 6, -1, -1, 9, 8, -1, -1, 11, -1, -1 },
        { 2, 6, 5, 7, 8, -1, 9, 10, 11, -1, -1, -1 },
        { 3, -1, 7, 4, -1, 10, -1, 8, -1, -1, 11, -1 },
        { 4, -1, 8, -1, -1, 11, -1, -1, -1, -1, -1, -1 },
        { 5, 9, -1, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
        { 6, 8, 9, -1, -1, -1, 11, -1, -1, -1, -1, -1 },
        { 7, -1, 10, 8, -1, -1, -1, 11, -1, -1, -1, -1 },
        { 8, -1, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 9, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 10, -1, -1, 11, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } } };

    do_lea(std::string("Example a2-b2=b2-c2:\n"), ex_a2_b2_b2_c2);

#if 0
	{
		// generate (a2-b2-(c3)-d4-e5)
		//
		for (int a = 0; a <= 2; a++)
		for (int b = 0; b <= 2; b++)
		for (int c = 0; c <= 3; c++)
		{
			std::array<int, 5> el1{a,b,c};
			for (int d = 0; d <= 3; d++)
			for (int e = 0; e <= 4; e++)
			for (int f = 0; f <= 5; f++)
			{
				std::array<int, 5> el2{0,0,d,e,f};

			}

		}
	}
#endif

    {
        atomic_MV_effect_algebra<3> MV_EA{
            { { 1, 2, 3 }, { 2, 2, 3 } }
        };

        std::cout << MV_EA.atoms.atom_ids[1] << '\n';

        std::cout << MV_EA.is_less_or_equal({ 0, 2, 1 }, { 1, 1, 2 }) << '\n';

        std::cout << MV_EA.are_summable({ 0, 2, 1 }, { 1, 1, 2 }) << '\n';

        int id = MV_EA.get_id({ 0, 0, 0 });
        std::cout << id << '\n';
        for (int i = 0; i < 36; i++) {
            std::cout << "Id = " << i << '\n';
            auto mv_ea_elem = MV_EA.get_elem_from_id(i);
            MV_EA.disp_elem(mv_ea_elem);
            std::cout << '\n';
            MV_EA.disp_elem_with_ids(mv_ea_elem);
            std::cout << '\n';
            std::cout << MV_EA.get_id(mv_ea_elem) << '\n';
        }

        std::cout << MV_EA.get_id({ 1, 0, 0 }) << '\n';
        std::cout << MV_EA.get_id({ 0, 2, 1 }) << '\n';
        std::cout << MV_EA.get_id({ 2, 2, 3 }) << '\n';
        std::cout << MV_EA.get_id({ -1 }) << '\n';

        MV_EA.generate_oplus();

        atomic_MV_effect_algebra<2> MV_EA2{
            { { 1, 2 }, { 2, 2 } }
        };

        MV_EA2.generate_oplus();

        std::array<std::array<int, 9>, 9> oplus333{
            { { 0, 1, 2, 3, 4, 5, 6, 7, 8 },
                { 1, 2, -1, 4, 5, -1, 7, 8, -1 },
                { 2, -1, -1, 5, -1, -1, 8, -1, -1 },
                { 3, 4, 5, 6, 7, 8, -1, -1, -1 },
                { 4, 5, -1, 7, 8, -1, -1, -1, -1 },
                { 5, -1, -1, 8, -1, -1, -1, -1, -1 },
                { 6, 7, 8, -1, -1, -1, -1, -1, -1 },
                { 7, 8, -1, -1, -1, -1, -1, -1, -1 },
                { 8, -1, -1, -1, -1, -1, -1, -1, -1 } }

        };

        do_lea("test", oplus333);
    }

    return 0;

    {
        std::fstream fs;

        fs.open("../ea11.txt", std::fstream::in);

        if (!fs.good()) {
            std::cout << "Error openning file ea11.txt\n";
        } else {
            char buffer[1024];
            int row = 0;
            int id = 1;
            std::array<std::array<int, 11>, 11> oplus;
            while (!fs.eof()) {
                fs.getline(buffer, 1024);
                if (!strncmp("Order", buffer, 5))
                    continue;

                if (!strncmp("Left", buffer, 4)) {
                    fs.getline(buffer, 1024);
                    continue;
                }

                if (!strncmp("Right", buffer, 5)) {
                    fs.getline(buffer, 1024);
                    continue;
                }

                if (strlen(buffer) == 0)
                    continue;

                if (!fs.good())
                    continue;

                if (row == 11) {
                    std::string msg = "Example " + std::to_string(id) + '\n';
                    do_lea2(msg, oplus);
                    row = 0;
                    id++;
                } else {
                    char* pom = buffer;
                    for (int col = 0; col < 11; col++) {
                        while (*pom && (*pom == ' ' || *pom == '\t'))
                            pom++;

                        if (*pom == '.') {
                            oplus[row][col] = -1;
                            pom++;
                            continue;
                        }

                        int num = 0;

                        char* pom_beg = pom;
                        while (*pom && (*pom >= '0' && *pom <= '9')) {
                            num = 10 * num + static_cast<int>(*pom - '0');
                            pom++;
                        }
                        if (pom_beg != pom)
                            oplus[row][col] = num;
                    }
                    row++;
                }
            }
        }

        fs.close();
    }

    return 0;

    std::cout << "Associativity: ";
    if (is_associative(ex_3_15_oplus, -1))
        std::cout << "OK";
    else
        std::cout << "Not ok.";
    std::cout << '\n';

    std::array<std::array<int, 18>, 18> ex_3_15_order{ { { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1 },
        { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 },
        { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 },
        { 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 } } };

    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 18; j++) {
            int v;
            if (get_inf(ex_3_15_order, i, j, v)) {
                std::cout << std::setw(3) << v << ", ";
            } else
                std::cout << std::setw(3) << '?' << ", ";
        }
        std::cout << '\n';
    }

#if 0
	int v;
	get_inf(ex_3_15_order, 10, 6, v);	
	std::cout << v << '\n';
#endif
    return 0;
}
