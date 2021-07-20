// *****************************************************************
// ASTAR
// *****************************************************************
namespace ASTAR
{
	static const int SIZEX = IMAGESCALE;
	static const int SIZEY = IMAGESCALE;

	static const int MAX_RAD = 150;		// 最大辐射阈值
	static real blendfactor = 0.5;		// 场与距离H函数的混合值

	// nearest path
	using mov_t = std::pair<char, char>;
	struct node_t
	{
		int parent = -1;
		point_t p;
		real g_field;
		real g;
		real h_field;
		real h;
	};

	PNG* png = 0, * field_png = 0;
	static byte hismap[SIZEX][SIZEY] = { 0 };
	std::vector< node_t> openlist;
	std::vector< node_t> closelist;
	node_t work_node;
	std::vector< point_t> bestway;

	// -----------------------------------------------------------------
	// 碰撞与辐射场测试
	// -----------------------------------------------------------------
	bool(*ptestcolor)(int x, int y);
	bool(*ptestfieldcolor)(int x, int y);
	real(*pgetfield)(int x, int y);

	bool testcolor(int x, int y)
	{
		return png && png->getcolor(x, y) == 0;
	}
	bool testfieldcolor(int x, int y)
	{
		return field_png && field_png->getcolor(x, y) > MAX_RAD;
	}
	real getfield(int x, int y)
	{
		if (!field_png)
			return 0;
		return 0.1 * field_png->getcolor(x, y);
	}
	// -----------------------------------------------------------------
	inline real getdis(int x, int y, const point_t& B)
	{
		return point_t::dis(point_t(x, y), B);
	}
	inline real getH_field(int x, int y, const point_t& B)
	{
		if (!pgetfield)
			return 0;
		real sum = 0;
		for (int i = x; i < B.x; i++)
			for (int j = y; j < B.y; j++)
			{
				sum += (*pgetfield)(i, j);
			}
		int w = abs(B.x - x);
		int h = abs(B.y - y);
		real k = h + w - 2 * _MIN(w, h) + 1.414 * _MIN(w, h);

		real s = abs(w + 1) * abs(h + 1);

		return k * sum / s;
	}
	inline real getH(int x, int y, const point_t& B)
	{
		real dis = getdis(x, y, B);
		return dis;
	}
	inline real sumH(real H_field, real H)
	{
		//return H_field;
		return blend(H, H_field, blendfactor);
		//return sqrt(H_field * H_field + H * H);
	}
	inline real sumG(real G_field, real G)
	{
		//return G_field;
		return blend(G, G_field, blendfactor);
		//return sqrt(G_field * G_field + G * G);
	}
	int findmin()
	{
		int minpos = -1;
		real min = 1e5;
		//PRINT("openlist=" << openlist.size());
		for (int i = 0; i < openlist.size(); i++)
		{
			const node_t& nd = openlist[i];
			real sum = sumH(nd.h_field, nd.h) + sumG(nd.g_field, nd.g);
			//PRINT("sum=" << sum);
			if (sum < min)
			{
				min = sum;
				minpos = i;
			}
		}
		return minpos;
	}

	bool path(node_t& node, const point_t& B, int depth)
	{
		work_node.g = -1;

		if (depth > 1500)
		{
			PRINT("depth!");
			work_node = node;
			return false;
		}

		//pointi(node.p.x, node.p.y, 1, 0xFFFFFFFF);

		if (getdis(node.p.x, node.p.y, B) < 1)
		{
			PRINT("GOAL!");
			return true;
		}

		static point_t movdirs[] =
		{
			point_t(0, 1),
			point_t(0, -1),
			point_t(1, 0),
			point_t(-1, 0),

			point_t(-1, -1),
			point_t(-1, 1),
			point_t(1, -1),
			point_t(1, 1)
		};

		for (auto it : movdirs)
		{
			static point_t np;
			np.x = node.p.x + it.x, np.y = node.p.y + it.y;
			if (np.x < 0 || np.x > SIZEX - 1 || np.y < 0 || np.y > SIZEY - 1)
			{
				continue;
			}
			{
				if (hismap[np.x][np.y] != 0)
					continue;
				hismap[np.x][np.y] = 1;
			}
			if ((*ptestcolor)(np.x, np.y) || (*ptestfieldcolor)(np.x, np.y))
			{
				PRINT("HIT!");
				continue;
			}

			static node_t n;
			{
				n.p = np;
				n.g = node.g + it.len();
				n.g_field = node.g_field + (*pgetfield)(np.x, np.y);
				n.h = getH(n.p.x, n.p.y, B);
				n.h_field = getH_field(np.x, np.y, B);

				n.parent = closelist.size() - 1;
			}
			{
				//bool in = false;
				//{ // check in openlist
				//	for (int i = 0; i < openlist.size(); i++)
				//	{
				//		if (openlist[i].p == n.p)
				//		{
				//			openlist[i].g = n.g;
				//			openlist[i].g_field = n.g_field;
				//			in = true;
				//			break;
				//		}
				//	}
				//}
				//if (!in)
				{
					pointi(n.p.x, n.p.y, 1, 0xFFFFFFFF);
					openlist.push_back(n);
				}
			}
		}
		int bestind = findmin();
		if (bestind == -1)
		{
			//PRINT("bestind=" << bestind)
			PRINT("findpath failed!");
			return false;
		}
		else
		{
			node_t best = openlist[bestind];
			pointi(best.p.x, best.p.y, 1, 0xFF808080);

			openlist.erase(openlist.begin() + bestind);
			closelist.push_back(best);

			return path(best, B, depth + 1);
		}
	}

	void findpath(const point_t& a, const point_t& b)
	{
		PRINT("-----------findpath------------")
			memset(hismap, 0, sizeof(hismap));
		node_t A;
		A.p = a;
		A.g = 0;
		A.g_field = 0;
		A.h = getH(A.p.x, A.p.y, b);
		A.h_field = getH_field(A.p.x, A.p.y, b);
		point(A.p, 4, 0xFF00FF00);

		node_t B;
		B.p = b;
		point(B.p, 4, 0xFF00FFFF);

		openlist.clear();
		closelist.clear();

		openlist.push_back(A);
		work_node = A;
		if (path(work_node, B.p, 0))
		{
			if (!closelist.empty())
			{
				int it = closelist.size() - 1;
				while (closelist[it].parent != -1)
				{
					pointi(closelist[it].p.x, closelist[it].p.y, 1, 0xFF00ff00);
					bestway.push_back(closelist[it].p);
					it = closelist[it].parent;
				}
			}
		}
	}
	bool init()
	{
		if (!png)
		{
			png = new PNG();
			if (!(*png).load("C:\\Users\\18858\\Documents\\LAB\\ZEXE/test.png"))
			{
				PRINT("!LoadFromFileA(C:\\Users\\18858\\Documents\\LAB\\ZEXE/test.png)");
				return false;
			}
		}/*
		if (!field_png)
		{
			field_png = new PNG();
			if (!(*field_png).load("C:\\Users\\18858\\Documents\\LAB\\ZEXE/field.png"))
			{
				PRINT("!LoadFromFileA(C:\\Users\\18858\\Documents\\LAB\\ZEXE/field.png)");
				return false;
			}
		}*/

		ptestcolor = testcolor;
		ptestfieldcolor = testfieldcolor;
		pgetfield = getfield;

		bestway.clear();

		blendfactor = luaparam[0];

		return true;
	}
	void test()
	{
		if (!init())
		{
			return;
		}

		findpath(pnt(62, 52), pnt(152, 119));
	}
}
