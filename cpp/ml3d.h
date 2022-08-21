#ifndef ml3d_h
#define ml3d_h

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#include <algorithm>
#include "raster.h"

class SliceML : public Slice
{
public:
	int z, idx;
	int ptid;
	double sc;
	SliceML() : Slice(), z(-1), idx(-1), ptid(-1), sc(0.) {}
	bool valid() { return area>0; }
	bool taken() { return ptid>=0; }
};

struct SortScoreML
{
	int ptid;
	int j;
	double sc;
	SortScoreML() : ptid(-1), j(-1), sc(0.) {}
	SortScoreML(int _ptid, double _sc) : ptid(_ptid), j(-1), sc(_sc) {}
};

struct SeedML
{
	int ptid;
	double sc;
	std::vector<std::vector<int>> zmap;
	std::vector<std::vector<SortScoreML>> lnkmap;
	SeedML() {}
	SeedML(int d) : zmap(d), lnkmap(d) {}
	void init(int d) { zmap.resize(d); lnkmap.resize(d); }
	int bottom() {
		for (int i=0; size_t(i)<zmap.size(); i++)
			if (!zmap[i].empty()) return i;
		return -1;
	}
	int top() {
		for (int i=int(zmap.size())-1; i>=0; i--)
			if (!zmap[i].empty()) return i;
		return -1;
	}
};

class AssemblerML
{
public:
	double GOOD_IOU = 0.75;
	double ACCEPTABLE_IOU = 0.65;
	int MAX_GAP = 3;
	int MIN_HEIGHT = 3;
	double DNA_PHANTOM_FILL = 0.6;

	int w, h, d;
	Raster3D mstack;
	Raster8 msk;
	std::vector<std::vector<SliceML>> all_slices;
	std::vector<Particle3D> cells;
	std::vector<SortScoreML> sorted_cells;
	
	std::vector<SeedML> seeds;
	
	AssemblerML(unsigned char *mask3d, int zm3d, int hm3d, int wm3d) :
		w(wm3d), h(hm3d), d(zm3d),
		mstack(wm3d, hm3d, zm3d, mask3d),
		msk(wm3d, hm3d, NULL),
		all_slices(zm3d)
	{
		
	}

	void load_flat_slices(std::vector<std::vector<int>> & all_flat_particles);
	void load_particles_3d(std::vector<std::vector<std::vector<int>>> & particles_3d);
	
	void find_seeds();
	void consolidate_seeds();
	void cells_from_seeds();
	void paint_final();
	
	void fill_gaps();
	void fix_borders_actin();
	void fix_borders_dna();
	
protected:
	void smooth_slice(SliceML& pt, int z=0);
	void find_z_match(std::vector<int>& res, SliceML& pt, int z, double min_iou, bool any=false);
	void assign_ptids_to_seed(SeedML& seed, std::vector<int>& ptids, int z);
	int assign_below(SeedML& seed, int z, double min_iou);
	void dedupe_seed(SeedML& seed);
	SortScoreML best_lnk_at(SeedML& seed, int z, int idx);
	int find_seed_match(SliceML& pt, int z, int *pptid);
	void attach_seed_slice(SeedML& seed, int z, int idx);
	bool merge_seed(SeedML& seed, SeedML& seed1);
	
	void interpolate_frame(Particle3D& cell, int zmiss);
	
};

#endif
