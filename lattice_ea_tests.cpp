#include "at_lat_ea/lattice_ea.h"
#include "at_lat_ea/lattice_ea_ops.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>

int main()
{
    std::cout << block_based_elem_type{ 100, 12939 }.toString() << '\n';
    Lattice_EA lea1{
        { { 3, { 1, 2, 3 }, { 2, 2, 3 } }, { 3, { 3, 4, 5 }, { 3, 4, 5 } } }
    };

    std::cout << lea1.get_block({ 1, 2 }).size() << "  "
              << lea1.get_block({ 1, 3 }).size() << "  "
              << lea1.get_block({ 2, 4 }).size() << "  "
              << lea1.get_block({ 4, 3 }).size() << "  " << '\n';

    /*
	Lattice_EA lea2{
		{{2, {1,2}, {2,2}}, {2, {2,3}, {2,2}}}
	};
	*/

    /*
    Lattice_EA lea2{
        { { 3, { 1, 2, 3 }, { 2, 2, 3 } },
            { 3, { 3, 4, 5 }, { 3, 2, 1 } },
            { 4, { 4, 6, 7, 8 }, { 2, 2, 1, 2 } } }
    };
    */

    /*
    Lattice_EA lea2{
        { { 3, { 1, 2, 3 }, { 1, 2, 1 } },
            { 3, { 3, 4, 5 }, { 1, 2, 1 } },
            { 3, { 1, 5, 6 }, { 1, 1, 2 } } }
    };
    */

    Lattice_EA lea2{
        { { 3, { 1, 5, 2 }, { 1, 1, 1 } },
            { 3, { 2, 6, 3 }, { 1, 1, 1 } },
            { 3, { 3, 7, 4 }, { 1, 1, 1 } },
            { 3, { 4, 8, 1 }, { 1, 1, 1 } } }
    };

    Lattice_EA_ops l{ lea2 };

    std::cout << "Size = " << l.get_size() << '\n';

    l.generate_ops();

    // print op: order, inf
    std::cout << "Order:\n";
    for (int i = 0; i < l.get_size(); i++) {
        for (int j = 0; j < l.get_size(); j++) {
            std::cout << std::setw(2) << l.order(i, j) << " ";
            //			std::cout << std::setw(2) << l.inf(i,j) << " ";
            //        std::cout << std::setw(2) << l.impl(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << "\nInf:\n";
    for (int i = 0; i < l.get_size(); i++) {
        for (int j = 0; j < l.get_size(); j++) {
            std::cout << std::setw(2) << l.inf(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << "\nSup:\n";
    for (int i = 0; i < l.get_size(); i++) {
        for (int j = 0; j < l.get_size(); j++) {
            std::cout << std::setw(2) << l.sup(i, j) << " ";
        }
        std::cout << '\n';
    }

    std::cout << "\nOplus:\n";
    for (int i = 0; i < l.get_size(); i++) {
        for (int j = 0; j < l.get_size(); j++) {
            std::cout << std::setw(2) << l.oplus(i, j) << " ";
        }
        std::cout << '\n';
    }

    // test b->a <= (((a->b)->b)->a)
    for (int i = 0; i < l.get_size(); i++)
        for (int j = 0; j < l.get_size(); j++) {

            int r1 = l.impl(j, i);
            int r2 = l.impl(l.impl(l.impl(i, j), j), i);

            if (!l.order(r1, r2)) {
                auto bi = l.idx_to_elem(i),
                     bj = l.idx_to_elem(j);
                auto blocki = l.get_lea().blocks[bi.block_id],
                     blockj = l.get_lea().blocks[bj.block_id];

                std::cout << "a = " << bi.toString() << " = ";
                blocki.disp_elem_with_ids(blocki.get_elem_from_id(bi.elem_id));
                std::cout << '\n';

                std::cout << "b = " << bj.toString() << "  = ";
                blockj.disp_elem_with_ids(blockj.get_elem_from_id(bj.elem_id));
                std::cout << '\n';

                std::cout << "Problem!" << '\n';
            }
        }

    // display atom representation of elements
    for (int i = 0; i < l.get_size(); ++i) {
        auto bi = l.idx_to_elem(i);
        auto blocki = l.get_lea().blocks[bi.block_id];

        std::cout << "Id = " << std::setw(2) << i << ": " << bi.toString() << " = ";
        blocki.disp_elem_with_ids(blocki.get_elem_from_id(bi.elem_id));
        std::cout << '\n';
    }
    std::cout << "---------\n\n";

    // display all elements in blocks
    auto block_id = 0;
    for (const auto& block : l.get_lea().blocks) {
        for (int i = 0; i < block.get_max_id(); ++i) {
            std::cout << "{" << block_id << ", " << i << "} = ";
            block.disp_elem_with_ids(block.get_elem_from_id(i));
            std::cout << '\n';
        }
        std::cout << "--\n";
        ++block_id;
    }
    std::cout << "---------\n\n";

    std::cout << l.get_lea().are_same({ 2, 6 }, { 1, 5 }) << '\n';

#if 0
	std::cout << '\n';
	for (int b1 = 0; b1 < (int)lea2.blocks.size(); b1++)
		for (int e1 = 0; e1 < lea2.blocks[b1].get_max_id(); e1++)
			if (lea2.is_cannonical({b1, e1}))
			{
				for (int b2 = 0; b2 < (int)lea2.blocks.size(); b2++)
					for (int e2 = 0; e2 < lea2.blocks[b2].get_max_id(); e2++)
						if (lea2.is_cannonical({b2, e2}))
						{
							std::cout << "a = " << b1 << "," << e1 << '\n';
							std::cout << "b = " << b2 << "," << e2 << '\n';

							// b->a
							auto p1 = lea2.impl({b2, e2}, {b1, e1});
							/*				std::cout << '{' << p1.block_id << ", " << p1.elem_id << "}, ";*/
							auto p2 = lea2.impl(lea2.impl(lea2.impl({b1, e1}, {b2, e2}), {b2, e2}), {b1, e1});
							/*					std::cout << '{' << p2.block_id << ", " << p2.elem_id << "}, \n"; */

							if (!lea2.are_leq(p1, p2))
								std::cout << "Problem: a = " << '{' << p1.block_id << ", " << p1.elem_id << "}, "
										  << "b = " << '{' << p2.block_id << ", " << p2.elem_id << "}\n";
						}
				/*			std::cout << '\n';*/
			}
	return 0;
#endif
#if 0

	block_based_elem_type e1 = {0,2};
	block_based_elem_type e2 = {1,3};

	std::cout << "a = " << e1 << ", " << "b = " << e2 << '\n';

	std::cout << "b -> a  = " << lea2.impl(e2, e1) << '\n';
	std::cout << "a'" << lea2.orthosupplement(e1) << '\n';
	std::cout << "a^b" << lea2.inf(e1, e2) << '\n';
	std::cout << "a->b" << lea2.oplus({0,6}, {1,3}) << '\n';
	std::cout << "a -> b  = " << lea2.impl(e1, e2) << '\n';

	std::cout << "{0,6} = ";
        lea2.blocks[0].disp_elem_with_ids(lea2.blocks[0].get_elem_from_id(6));
	std::cout << '\n';

	std::cout << "{1,3} = ";
        lea2.blocks[1].disp_elem_with_ids(lea2.blocks[1].get_elem_from_id(3));
	std::cout << '\n';

	std::cout << "{1,5} = ";
        lea2.blocks[1].disp_elem_with_ids(lea2.blocks[1].get_elem_from_id(5));
	std::cout << '\n';

	std::cout << "{0,6} <= {1,5}: " << lea2.are_leq({0,6}, {1,5}) << '\n';


    // display elements
	for (int i1 = 0; i1 < (int)lea2.blocks.size(); i1++)
	for (int i2 = 0; i2 < lea2.blocks[i1].get_max_id(); i2++)
	{
		auto p = lea2.get_cannonical({i1,i2});
		if (p.block_id == i1 && p.elem_id == i2)
        {
			std::cout << '{' << p.block_id << ',' << p.elem_id << '}' << ' ';

            auto p_c = lea2.orthosupplement({i1, i2});

            std::cout << " -> ";
			std::cout << '{' << p_c.block_id << ',' << p_c.elem_id << '}' << ' ';
            std::cout << '\n';

        }
	}
	std::cout << '\n';


	std::vector<int> remove_elems;

	for (int i1 = 0; i1 < (int)lea2.blocks.size(); i1++)
	for (int i2 = 0; i2 < lea2.blocks[i1].get_max_id(); i2++)
	{
		lea2.blocks[i1].disp_elem_with_ids(lea2.blocks[i1].get_elem_from_id(i2));
		std::cout << '\n';
	}

	for (int i1 = 0; i1 < (int)lea2.blocks.size(); i1++)
	for (int i2 = 0; i2 < lea2.blocks[i1].get_max_id(); i2++)
	{
		for (int j1 = 0; j1 < (int)lea2.blocks.size(); j1++)
		for (int j2 = 0; j2 < lea2.blocks[j1].get_max_id(); j2++)
		{
			bool same = lea2.are_same({i1,i2}, {j1,j2}); 	
			if (same)
			{
				if (j2 == lea2.blocks[j1].get_max_id() - 1)
				{
					if (lea2.get_id_total({i1,i2}) < lea2.get_id_total({j1,j2}))
						remove_elems.push_back(lea2.get_id_total({i1,i2}));
				}
				else
				{
				if (lea2.get_id_total({i1,i2}) < lea2.get_id_total({j1,j2}))
					remove_elems.push_back(lea2.get_id_total({j1,j2}));
				}
			}
			std::cout << same << " ";
		}
		std::cout << '\n';
	}

	std::sort(remove_elems.begin(), remove_elems.end());
	for (const auto &v: remove_elems)
		std:: cout << v << " ";

	std::cout << '\n';


	std::cout << '\n';
	for (int i1 = 0; i1 < (int)lea2.blocks.size(); i1++)
	for (int i2 = 0; i2 < lea2.blocks[i1].get_max_id(); i2++)
	{
		if (std::find(remove_elems.begin(), remove_elems.end(), lea2.get_id_total({i1, i2})) == remove_elems.end())
		{
			for (int j1 = 0; j1 < (int)lea2.blocks.size(); j1++)
			for (int j2 = 0; j2 < lea2.blocks[j1].get_max_id(); j2++)
			{
				if (std::find(remove_elems.begin(), remove_elems.end(), lea2.get_id_total({j1, j2})) == remove_elems.end())
					std::cout << lea2.are_leq({i1,i2}, {j1,j2}) << " ";
			}
			std::cout << '\n';
		}
	}

	// define map_compress
	std::map<std::pair<int,int>, int> compress_map;
	int id = 0;

	for (int i1 = 0; i1 < (int)lea2.blocks.size(); i1++)
	for (int i2 = 0; i2 < lea2.blocks[i1].get_max_id(); i2++)
	{
		if (std::find(remove_elems.begin(), remove_elems.end(), lea2.get_id_total({i1, i2})) == remove_elems.end())
		{
			compress_map.insert({{i1,i2}, id});
			id++;

		}
	}

	// oplus
	std::cout << '\n';
	for (int i1 = 0; i1 < (int)lea2.blocks.size(); i1++)
	for (int i2 = 0; i2 < lea2.blocks[i1].get_max_id(); i2++)
	{
		if (std::find(remove_elems.begin(), remove_elems.end(), lea2.get_id_total({i1, i2})) == remove_elems.end())
		{
			for (int j1 = 0; j1 < (int)lea2.blocks.size(); j1++)
			for (int j2 = 0; j2 < lea2.blocks[j1].get_max_id(); j2++)
			{
				if (std::find(remove_elems.begin(), remove_elems.end(), lea2.get_id_total({j1, j2})) == remove_elems.end())
				{
					auto res = lea2.oplus({i1,i2}, {j1,j2});
#if 0
					if (res.block_id == -1)
						std::cout << "{}, ";
					else
						std::cout << '{' << res.block_id << "," <<res.elem_id << "},";
#endif
					if (res.block_id == -1)
						std::cout << std::setw(4) << -1;
					else
						std::cout << std::setw(4) << compress_map[{res.block_id, res.elem_id}];
				}

			}
			std::cout << '\n';lattice_ea
		}
	}

	std::cout << "\nInf: \n";
	for (int b1 = 0; b1 < (int)lea2.blocks.size(); b1++)
	for (int e1 = 0; e1 < lea2.blocks[b1].get_max_id(); e1++)
		if (lea2.is_cannonical({b1, e1}))
		{
			for (int b2 = 0; b2 < (int)lea2.blocks.size(); b2++)
			for (int e2 = 0; e2 < lea2.blocks[b2].get_max_id(); e2++)
				if (lea2.is_cannonical({b2, e2}))
				{
					auto p = lea2.inf({b1,e1}, {b2,e2});
					std::cout << '{' << p.block_id << ", " << p.elem_id << "}, ";
				}
			std::cout << '\n';

		}

	std::cout << "\nSup: \n";
	for (int b1 = 0; b1 < (int)lea2.blocks.size(); b1++)
	for (int e1 = 0; e1 < lea2.blocks[b1].get_max_id(); e1++)
		if (lea2.is_cannonical({b1, e1}))
		{
			for (int b2 = 0; b2 < (int)lea2.blocks.size(); b2++)
			for (int e2 = 0; e2 < lea2.blocks[b2].get_max_id(); e2++)
				if (lea2.is_cannonical({b2, e2}))
				{
					auto p = lea2.sup({b1,e1}, {b2,e2});
					std::cout << '{' << p.block_id << ", " << p.elem_id << "}, ";
				}
			std::cout << '\n';

		}
#endif

#if 0
	// test b->a <= (((a->b)->b)->a)
	std::cout << '\n';
	for (int b1 = 0; b1 < (int)lea2.blocks.size(); b1++)
		for (int e1 = 0; e1 < lea2.blocks[b1].get_max_id(); e1++)
			if (lea2.is_cannonical({b1, e1}))
			{
				for (int b2 = 0; b2 < (int)lea2.blocks.size(); b2++)
					for (int e2 = 0; e2 < lea2.blocks[b2].get_max_id(); e2++)
						if (lea2.is_cannonical({b2, e2}))
						{
							std::cout << "a = " << b1 << "," << e1 << '\n';
							std::cout << "b = " << b2 << "," << e2 << '\n';

							// b->a
							auto p1 = lea2.impl({b2, e2}, {b1, e1});
							/*				std::cout << '{' << p1.block_id << ", " << p1.elem_id << "}, ";*/
							auto p2 = lea2.impl(lea2.impl(lea2.impl({b1, e1}, {b2, e2}), {b2, e2}), {b1, e1});
							/*					std::cout << '{' << p2.block_id << ", " << p2.elem_id << "}, \n"; */

							if (!lea2.are_leq(p1, p2))
								std::cout << "Problem: a = " << '{' << p1.block_id << ", " << p1.elem_id << "}, "
										  << "b = " << '{' << p2.block_id << ", " << p2.elem_id << "}\n";
						}
				/*			std::cout << '\n';*/
			}
#endif

#if 0
// inf test
	std::cout << '\n';
	for (int b1 = 0; b1 < (int)lea2.blocks.size(); b1++)
		for (int e1 = 1; e1 < lea2.blocks[b1].get_max_id() - 1 ; e1++)
			if (lea2.is_cannonical({b1, e1}))
			{
				for (int b2 = 0; b2 < (int)lea2.blocks.size(); b2++)
					for (int e2 = 1; e2 < lea2.blocks[b2].get_max_id() - 1; e2++)
						if (lea2.is_cannonical({b2, e2}))
						{
							if (b1 != b2)
							{
								std::cout << "a = " << b1 << "," << e1 << '\n';
								std::cout << "b = " << b2 << "," << e2 << '\n';

								// b->a
								auto p = lea2.inf({b2, e2}, {b1, e1});
								std::cout << "inf = " << p << " %=% ";
								lea2.blocks[p.block_id].disp_elem_with_ids(
									lea2.blocks[p.block_id].get_elem_from_id(p.elem_id)
								);
								std::cout << '\n';
							}
						}
				/*			std::cout << '\n';*/
			}

#endif
    return 0;
}
