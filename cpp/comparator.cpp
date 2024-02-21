
#include "comparator.h"

std::vector<SliceMatch> SegmentationComparator3D::init_matches(std::vector<Particle3D>& cells, int z)
{
	std::vector<SliceMatch> res;
	
	for (int cell_id=0; size_t(cell_id)<cells.size(); cell_id++) {
		std::vector<HSeg> &fill = cells[cell_id].fills[z];
		if (fill.empty()) continue;
		long long area = fill_area(fill);
		if (area >= 10) {
			res.push_back(SliceMatch(cell_id, z, area));
		}
	}
	
	return res;
}

std::vector<CellMatch> SegmentationComparator3D::init_cell_matches(
		std::vector<Particle3D>& cells, int min_height, double min_score)
{
	std::vector<CellMatch> res;
	
	for (int cell_id=0; size_t(cell_id)<cells.size(); cell_id++) {
		Particle3D& cell = cells[cell_id];
		long long vol = cell.update_from_fill();
		int height = cell.bnd.zmax - cell.bnd.zmin + 1;
		if (height < min_height) continue;
		double sc = cell.iou_score(3);
		if (sc >= min_score) {
			res.push_back(CellMatch(cell_id, vol));
		}
	}
	
	return res;
}

static int count_back_matches(std::vector<SliceMatch>& sms, int my_id)
{
	int cnt = 0;
	for (SliceMatch& sm : sms)
		if (sm.other_id == my_id) ++cnt;
	return cnt;
}

void SegmentationComparator3D::compare_slices_z(int z)
{
	std::vector<SliceMatch> bsms = init_matches(base_cells, z);
	std::vector<SliceMatch> csms = init_matches(cmp_cells, z);
	for (int b_id=0; size_t(b_id)<bsms.size(); b_id++) {
		SliceMatch& bsm = bsms[b_id];
		std::vector<HSeg> &fill = base_cells[bsm.cell_id].fills[z];
		for (int c_id=0; size_t(c_id)<csms.size(); c_id++) {
			SliceMatch& csm = csms[c_id];
			std::vector<HSeg> &cfill = cmp_cells[csm.cell_id].fills[z];
			long long ovl = fill_overlay_area(fill, cfill);
			long long ovl2 = ovl + ovl;
			if (ovl2 < bsm.my_area && ovl2 < csm.my_area) continue;
			
			++bsm.nm;
			bsm.other_id = c_id;
			bsm.other_area += csm.my_area;
			bsm.ovl += ovl;
			
			++csm.nm;
			csm.other_id = b_id;
			csm.other_area += bsm.my_area;
			csm.ovl += ovl;
		}
	}

	for (int b_id=0; size_t(b_id)<bsms.size(); b_id++) {
		SliceMatch& bsm = bsms[b_id];
		if (bsm.nm == 0) {
			// Not detected -- false negative
			++f_neg;
		} else if (bsm.nm > 1) {
			// Multiple matches -- fused
			++fused;
		} else {
			int back_nm = count_back_matches(csms, b_id);
			if (back_nm > 1) {
				// Maps to multiple -- fragmented
				++fragm;
			} else {
				// 1-1 match
				long long uni = bsm.my_area + bsm.other_area - bsm.ovl;
				int idx = int(double(bsm.ovl) * 100. / double(uni));
				if (idx >= 100) idx = 99;
				//else if (idx < 50) idx = 50;
				++pct_match[idx];
			}
		}
	}

	for (int c_id=0; size_t(c_id)<csms.size(); c_id++) {
		SliceMatch& csm = csms[c_id];
		if (csm.nm == 0) {
			// Not matched -- false positive
			++f_pos;
		}
	}
}

static int count_back_matches_3d(std::vector<CellMatch>& cms, int my_id)
{
	int cnt = 0;
	for (CellMatch& cm : cms) {
		if (cm.other_id == my_id)
			++cnt;
	}
	return cnt;
}

void SegmentationComparator3D::compare_cells(int min_height, double min_score)
{
	std::vector<CellMatch> base_cell_matches = init_cell_matches(base_cells, min_height, min_score);
	std::vector<CellMatch> cmp_cell_matches = init_cell_matches(cmp_cells, min_height, min_score);
	
	num_valid_base_cells = int(base_cell_matches.size());
	num_valid_cmp_cells = int(cmp_cell_matches.size());
	
	for (int base_midx=0; size_t(base_midx)<base_cell_matches.size(); base_midx++) {
		CellMatch& bcm = base_cell_matches[base_midx];
		Particle3D& bcell = base_cells[bcm.cell_id];
		for (int cmp_midx=0; size_t(cmp_midx)<cmp_cell_matches.size(); cmp_midx++) {
			CellMatch& ccm = cmp_cell_matches[cmp_midx];
			Particle3D& ccell = cmp_cells[ccm.cell_id];
			if (!bcell.bnd.intersects(ccell.bnd)) continue;
			long long ivol = bcell.overlay_volume(ccell);
			long long ivol2 = ivol + ivol;
			if (ivol2 < bcm.my_vol && ivol2 < ccm.my_vol) continue;
			
			++bcm.nm;
			bcm.other_id = ccm.cell_id;
			bcm.other_vol += ccm.my_vol;
			bcm.other_midx = cmp_midx;
			bcm.ovl += ivol;
			
			++ccm.nm;
			ccm.other_id = bcm.cell_id;
			ccm.other_vol += bcm.my_vol;
			ccm.other_midx = base_midx;
			ccm.ovl += ivol;
		}
	}
	
	for (CellMatch& bcm : base_cell_matches) {
		if (bcm.nm == 0) {
			// Not detected -- false negative
			++f_neg_3d;
		} else if (bcm.nm > 1) {
			// Multiple matches -- fused
			++fused_3d;
		} else {
			int back_nm = count_back_matches_3d(cmp_cell_matches, bcm.cell_id);
			if (back_nm > 1) {
				// Maps to multiple -- fragmented
				++fragm_3d;
			} else {
				// 1-1 match
				CellMatch& ccm = cmp_cell_matches[bcm.other_midx];
				long long uni = bcm.my_vol + ccm.my_vol - bcm.ovl;
				int idx = int(double(bcm.ovl) * 100. / double(uni));
				if (idx >= 100) idx = 99;
				//else if (idx < 50) idx = 50;
				++pct_match_3d[idx];
			}
		}
	}
	
	for (CellMatch& ccm : cmp_cell_matches) {
		if (ccm.nm == 0) {
			// Not matched -- false positive
			++f_pos_3d;
		}
	}
}




