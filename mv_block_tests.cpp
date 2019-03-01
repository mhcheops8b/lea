#include "at_lat_ea/mv_block.h"
#include <iostream>

int main()
{

    MV_Block b{ 3, { 1, 2, 3 }, { 2, 2, 3 } };

    std::cout << b.size << b.ids[0] << b.ids[1] << b.ids[2] << b.orders[0] << b.orders[1] << b.orders[2];
    std::cout << b.get_max_id();

    std::cout << '\n'
              << b.orders.size() << '\n';

    for (int i = 0; i < b.get_max_id(); i++) {
        std::cout << "Id = " << i << '\n';
        auto mv_ea_elem = b.get_elem_from_id(i);
        b.disp_elem(mv_ea_elem);
        std::cout << '\n';
        b.disp_elem_with_ids(mv_ea_elem);
        std::cout << '\n';
        std::cout << b.get_id(mv_ea_elem) << '\n';
    }

    MV_Block b2{ 3, { 3, 2, 1 }, { 3, 2, 2 } };

    std::cout << b2.size << '\n';

    for (int i = 0; i < b2.get_max_id(); i++) {
        MV_Block::annotated_type el = b2.get_annotated_elem_from_id(i);

        for (const auto& u : el) {
            std::cout << '{' << u.first << ',' << u.second << "},";
        }
        std::cout << '\n';
    }
}
