#ifndef comparator_h
#define comparator_h

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#include <sstream>
#include "geom.h"

struct SliceMatch
{
	int cell_id;
	int z;
	long long my_area;
	int nm;
	int other_id;
	long long other_area;
	long long ovl;
	SliceMatch() {}
	SliceMatch(int _cell_id, int _z, long long _my_area) : cell_id(_cell_id), z(_z), my_area(_my_area) {
		nm = 0;
		other_id = -1;
		other_area = 0;
		ovl = 0;
	}
};

struct CellMatch
{
	int cell_id;
	long long my_vol;
	int nm;
	int other_id;
	int other_midx;
	long long other_vol;
	long long ovl;
	CellMatch() {}
	CellMatch(int _cell_id, long long _my_vol) : cell_id(_cell_id), my_vol(_my_vol) {
		nm = 0;
		other_id = other_midx = -1;
		other_vol = ovl = 0ll;
	}
};

class SegmentationComparator3D
{
public:
	int w, h, d;
	const char *base_csv;
	const char *cmp_csv;
	
	std::vector<Particle3D> base_cells;
	std::vector<Particle3D> cmp_cells;
	
	int num_base_slices;
	int num_cmp_slices;
	
	// 2D stats
	int f_pos, f_neg, fragm, fused;
	int pct_match[100];
	
	// 3D stats
	int num_valid_base_cells;
	int num_valid_cmp_cells;
	int f_pos_3d, f_neg_3d, fragm_3d, fused_3d;
	int pct_match_3d[100];
	
	SegmentationComparator3D(int _w, int _h, int _d, const char *_base_csv, const char *_cmp_csv) :
		w(_w), h(_h), d(_d), base_csv(_base_csv), cmp_csv(_cmp_csv)
	{
		read_cell_data(base_csv, base_cells, w, h, d);
		read_cell_data(cmp_csv, cmp_cells, w, h, d);
		
		num_base_slices = count_slices(base_cells);
		num_cmp_slices = count_slices(cmp_cells);
		
		clear_stats();
	}
	
	void clear_stats() {
		f_pos = f_neg = fragm = fused = 0;
		memset(pct_match, 0, sizeof pct_match);
		num_valid_base_cells = num_valid_cmp_cells = 0;
		f_pos_3d = f_neg_3d = fragm_3d = fused_3d = 0;
		memset(pct_match_3d, 0, sizeof pct_match_3d);
	}
	
	static int count_slices(std::vector<Particle3D>& cells) {
		int rc = 0;
		for (Particle3D& cell : cells) {
			for (std::vector<HSeg> &fill : cell.fills) {
				if (!fill.empty()) ++rc;
			}
		}
		return rc;
	}
	
	static std::vector<SliceMatch> init_matches(std::vector<Particle3D>& cells, int z);

	void compare_slices_z(int z);
	void compare_slices() {
		for (int z=0; z<d; z++)
			compare_slices_z(z);
	}
	
	static std::vector<CellMatch> init_cell_matches(std::vector<Particle3D>& cells, int min_height, double min_score);
	void compare_cells(int min_height=3, double min_score=3.);
	
};

#endif
