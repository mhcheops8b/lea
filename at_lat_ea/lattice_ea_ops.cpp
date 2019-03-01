#include "lattice_ea_ops.h"

int Lattice_EA_ops::elem_arr_idx(int row, int col)
{
    return row * size + col;
}

void Lattice_EA_ops::init()
{

    int id = 0;
    // initialize maps
    for (int b1 = 0; b1 < (int)lea.blocks.size(); b1++)
        for (int e1 = 0; e1 < lea.blocks[b1].get_max_id(); e1++)
            if (lea.is_cannonical({ b1, e1 })) {
                map_idx_elem.insert({ id, { b1, e1 } });
                map_elem_idx.insert({ { b1, e1 }, id });
                id++;
            }
    // insert undefined
    map_idx_elem.insert({ -1, { -1, -1 } });
    map_elem_idx.insert({ { -1, -1 }, -1 });

    // initialize size
    size = id;
}

void Lattice_EA_ops::generate_ops()
{

    std::cout << "Initializing..." << '\n';
    // allocated memory
    arr_order.resize(size * size);
    arr_inf.resize(size * size);
    arr_sup.resize(size * size);
    arr_oplus.resize(size * size);

    vec_orthosupplement.resize(size);
    int arr_idx = 0;
    int vec_idx = 0;
    for (int i = 0; i < size; i++) {

        arr_idx += i;
        std::cout << i << " ";
        for (int j = i; j < size; j++) {
            // initialize arr_ord
            if (i == j) {
                arr_order[arr_idx] = 1;
            } else {
                // i ! = j
                if (lea.are_leq(map_idx_elem[i], map_idx_elem[j])) {
                    // i < j
                    arr_order[arr_idx] = 1;
                    arr_order[elem_arr_idx(j, i)] = 0;

                } else {
                    if (lea.are_leq(map_idx_elem[j], map_idx_elem[i])) {
                        // j < i
                        arr_order[elem_arr_idx(j, i)] = 1;
                        arr_order[arr_idx] = 0;
                    } else {
                        arr_order[arr_idx] = 0;
                        arr_order[elem_arr_idx(j, i)] = 0;
                    }
                }
            }
            // initialize arr_inf
            arr_inf[arr_idx] = map_elem_idx[lea.inf(map_idx_elem[i], map_idx_elem[j])];
            // initialize arr_sup
            arr_sup[arr_idx] = map_elem_idx[lea.sup(map_idx_elem[i], map_idx_elem[j])];
            // initialize arr_oplus
            arr_oplus[arr_idx] = map_elem_idx[lea.oplus(map_idx_elem[i], map_idx_elem[j])];
            if (i != j) {
                arr_inf[elem_arr_idx(j, i)] = arr_inf[arr_idx];
                arr_sup[elem_arr_idx(j, i)] = arr_sup[arr_idx];
                arr_oplus[elem_arr_idx(j, i)] = arr_oplus[arr_idx];
            }
            arr_idx++;
        }

        // initialize vec_orthosupplement;
        vec_orthosupplement[vec_idx] = map_elem_idx[lea.orthosupplement(map_idx_elem[i])];
        vec_idx++;
        std::flush(std::cout);
    }
    std::cout << "\nFinished" << '\n';
}

int Lattice_EA_ops::get_size()
{
    return size;
}

int Lattice_EA_ops::order(int el1, int el2)
{
    return arr_order[elem_arr_idx(el1, el2)];
}

int Lattice_EA_ops::oplus(int el1, int el2)
{
    return arr_oplus[elem_arr_idx(el1, el2)];
}

int Lattice_EA_ops::orthosupplement(int el)
{
    return vec_orthosupplement[el];
}

int Lattice_EA_ops::inf(int el1, int el2)
{
    return arr_inf[elem_arr_idx(el1, el2)];
}

int Lattice_EA_ops::sup(int el1, int el2)
{
    return arr_sup[elem_arr_idx(el1, el2)];
}

int Lattice_EA_ops::impl(int el1, int el2)
{
    return oplus(orthosupplement(el1), inf(el1, el2));
}

block_based_elem_type Lattice_EA_ops::idx_to_elem(int idx)
{
    return map_idx_elem[idx];
}

int Lattice_EA_ops::elem_to_idx(const block_based_elem_type& elem)
{
    return map_elem_idx[elem];
}
