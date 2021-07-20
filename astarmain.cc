// ---------------------------------------------------------------------
// 改进型A* 辐射场寻路算法
// ---------------------------------------------------------------------
namespace ASTAR
{
	static const int MAXDEPTH = 6000;		// 最大深度
	static const int SIZEX = IMAGESCALE;	// 分辨率
	static const int SIZEY = IMAGESCALE;
	static const int MAX_RAD = 50;			// 最大辐射阈值

	real blendfactor = 0.5;					// 场与距离H函数的混合值 0-忽略辐射场 1-只考虑辐射场

	// -----------------------------------------------------------------
	// nearest path
	using mov_t = std::pair<char, char>;
	struct node_t
	{
		int parent = -1;
		PM::point_t p;
		double g_field = 0;					// 累计辐射值
		double g = 0;						// 累计距离
		double h_field = 0;					// 估计辐射值	
		double h = 0;						// 估计距离
	};

	BYTE hismap[SIZEX][SIZEY] = { 0 };		// 遍历过的点图
	std::vector< node_t> openlist;			// 开放列表，有待遍历的
	std::vector< node_t> closelist;			// 关闭列表
	node_t A;								// 起始跟目标点
	node_t B;
	std::vector< PM::point_t> bestway;		// 最优路径

	// -----------------------------------------------------------------
	// 碰撞与辐射场测试
	bool(*ptestcolor)(int x, int y);
	bool(*ptestfieldcolor)(int x, int y);
	double(*pgetfield)(int x, int y);

	// -----------------------------------------------------------------
	inline real getdis(int x, int y, const PM::point_t& B)
	{
		return PM::point_t::dis(PM::point_t(x, y), B);
	}
	// 辐射量估计值
	inline double getH_field(int x, int y, const PM::point_t& B)
	{
		double sum = 0;
		for (int i = x; i < B.x; i++)
			for (int j = y; j < B.y; j++)
			{
				sum += (*pgetfield)(i, j);
			}
		int w = abs(B.x - x);
		int h = abs(B.y - y);
		double k = h + w - 2 * _MIN(w, h) + 1.414 * _MIN(w, h);

		double s = abs(w + 1) * abs(h + 1);

		return k * sum / s;
	}
	inline real getH(int x, int y, const PM::point_t& B)
	{
		real dis = getdis(x, y, B);
		return dis;
	}
	// 辐射量跟距离混合
	inline double sumH(double H_field, double H)
	{
		return blendd(H, H_field, blendfactor);
	}
	inline double sumG(double G_field, double G)
	{
		return blendd(G, G_field, blendfactor);
	}
	// 查询最优位置 g+h最小
	int findmin()
	{
		int minpos = -1;
		double min = 1e6;
		for (int i = 0; i < openlist.size(); i++)
		{
			const node_t& nd = openlist[i];
			double sum = sumH(nd.h_field, nd.h) + sumG(nd.g_field, nd.g);
			if (sum < min)
			{
				min = sum;
				minpos = i;
			}
		}
		return minpos;
	}
	// 九宫格各个方向
	static PM::point_t movdirs[] =
	{
		PM::point_t(0, 1),
		PM::point_t(0, -1),
		PM::point_t(1, 0),
		PM::point_t(-1, 0),

		PM::point_t(-1, -1),
		PM::point_t(-1, 1),
		PM::point_t(1, -1),
		PM::point_t(1, 1)
	};

	bool path(node_t& node, short depth)
	{
		if (depth > MAXDEPTH)
		{
			PRINT("!ERROR: depth > MAXDEPTH");
			return false;
		}
		if (getdis(node.p.x, node.p.y, B.p) < 1)
		{
			B.g = node.g;
			B.g_field = node.g_field;
			PRINT("GOAL!");
			return true;
		}
		for (short i = 0; i < 8; i++)
		{
			static PM::point_t np; np.x = node.p.x + movdirs[i].x, np.y = node.p.y + movdirs[i].y;

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
				//PRINT("HIT!");
				continue;
			}
			{
				static node_t n;
				n.p = np;
				n.g = node.g + movdirs[i].flen();
				n.g_field = node.g_field + (*pgetfield)(np.x, np.y);
				n.h = getH(n.p.x, n.p.y, B.p);
				n.h_field = getH_field(np.x, np.y, B.p);

				n.parent = closelist.size() - 1;
				openlist.push_back(n);
			}
		}
		static int bestind; bestind = findmin();
		if (bestind == -1)
		{
			PRINT("findpath failed！bestind=-1 depth=" << depth);
			closelist.clear();
			return false;
		}
		{
			static node_t best; best = openlist[bestind];

			openlist.erase(openlist.begin() + bestind);
			closelist.push_back(best);

			return path(best, depth + 1); // 尾部优化
		}
	}

	void findpath(const PM::point_t& a, const PM::point_t& b)
	{
		PRINT("astar findpath a=(" << a.x << "," << a.y << ") b=(" << b.x << "," << b.y << ")");
		memset(hismap, 0, sizeof(hismap));

		A.p = a;
		A.g = 0;
		A.g_field = 0;
		A.h = getH(A.p.x, A.p.y, b);
		A.h_field = getH_field(A.p.x, A.p.y, b);
		if ((*ptestcolor)(a.x, a.y) || (*ptestfieldcolor)(a.x, a.y))
		{
			PRINT("a点不可达！");
			return;
		}
		if ((*ptestcolor)(b.x, b.y) || (*ptestfieldcolor)(b.x, b.y))
		{
			PRINT("b点不可达！");
			return;
		}
		B.p = b;

		openlist.clear();
		closelist.clear();

		openlist.push_back(A);

		if (path(A, 0))
		{
			if (!closelist.empty())
			{
				int it = closelist.size() - 1;
				while (closelist[it].parent != -1)
				{
					bestway.push_back(closelist[it].p);
					it = closelist[it].parent;
				}
			}
		}
	}
}
