#include "ml3d.h"

void AssemblerML::load_flat_slices(std::vector<std::vector<int>> & all_flat_particles)
{
	for (int z=0; z<d; z++) {
		if (size_t(z) >= all_flat_particles.size()) break;
		std::vector<int> & flat_particles = all_flat_particles[z];
		std::vector<SliceML> & particles = all_slices[z];
		
		size_t iStart = 0;
		while (iStart < flat_particles.size()) {
			int hlen = flat_particles[iStart++];
			particles.resize(particles.size()+1);
			SliceML& ptc = particles[particles.size()-1];
			for (int i=0; i<hlen; i++) {
				ptc.fill.push_back(HSeg(flat_particles[iStart], flat_particles[iStart+1], flat_particles[iStart+2]));
				iStart += 3;
			}
			ptc.area = ptc.update_from_fill();
			if (ptc.area < 50) {
				particles.resize(particles.size()-1);
			} else {
				ptc.z = z;
				ptc.idx = int(particles.size()-1);
			}
		}
		
		// std::cout << "z=" << z << "; #slices=" << particles.size() << std::endl;
	}

}

void AssemblerML::load_particles_3d(std::vector<std::vector<std::vector<int>>> & particles_3d)
{
	for (int z=0; z<d; z++) {
		if (size_t(z) >= particles_3d.size()) break;
		std::vector<std::vector<int>>& particles_2d = particles_3d[z];
		std::vector<SliceML> & slices = all_slices[z];
		
		for (std::vector<int>& pt2d : particles_2d) {
			slices.resize(slices.size()+1);
			SliceML& ptc = slices[slices.size()-1];
			for (size_t i=0; i<pt2d.size(); i+=3) {
				if (i+2 >= pt2d.size()) break;
				ptc.fill.push_back(HSeg(pt2d[i], pt2d[i+1], pt2d[i+2]));
			}
			ptc.area = ptc.update_from_fill();
			if (ptc.area < 50) {
				slices.resize(slices.size()-1);
			} else {
				ptc.z = z;
				ptc.idx = int(slices.size()-1);
			}
		}
	}
}

void AssemblerML::find_seeds()
{
	msk.fill(0x10);
	for (int z=0; z<d; z++) {
		for (SliceML& pt : all_slices[z]) {
			smooth_slice(pt, z);
			pt.ptid = -1;
		}
	}
	
	for (int z=d-1; z>0; z--) {
		std::vector<SliceML>& slices = all_slices[z];
		
		for (int idx=0; size_t(idx)<slices.size(); idx++) {
			SliceML& pt = slices[idx];
			if (pt.taken()) continue;
			std::vector<int> ptids;
			find_z_match(ptids, pt, z-1, GOOD_IOU, false);
			if (ptids.empty()) continue;
			std::vector<int> ptids0;
			find_z_match(ptids0, pt, z, GOOD_IOU, false);
			int sid = int(seeds.size());
			seeds.resize(seeds.size()+1);
			SeedML& seed = seeds[seeds.size()-1];
			seed.init(d);
			seed.ptid = sid;
			assign_ptids_to_seed(seed, ptids0, z);
			assign_ptids_to_seed(seed, ptids, z-1);
			for (int zz=z; zz>0; zz--) {
				int cnt = assign_below(seed, zz, GOOD_IOU);
				if (cnt == 0) break;
			}
			dedupe_seed(seed);
		}
	}
}

static bool compare_cell_sort(SortScoreML& ss1, SortScoreML& ss2)
{
	return ss1.sc > ss2.sc;
}

void AssemblerML::consolidate_seeds()
{
	sorted_cells.resize(seeds.size());
	for (int idx=0; size_t(idx)<seeds.size(); idx++) {
		sorted_cells[idx].ptid = idx;
		sorted_cells[idx].sc = seeds[idx].sc;
	}
	std::sort(sorted_cells.begin(), sorted_cells.end(), compare_cell_sort);
	
	for (SortScoreML& ss : sorted_cells) {
		SeedML& seed = seeds[ss.ptid];
		if (seed.ptid < 0) continue;

		while (true) {
			int z0 = seed.bottom();
			SliceML& pt = all_slices[z0][seed.zmap[z0][0]];
			int sid=-1, ptid=-1, z1=-1;
			for (int z=z0-1; z>=0; z--) {
				if (z0 - z > MAX_GAP) break;
				sid = find_seed_match(pt, z, &ptid);
				if (ptid >= 0) {
					z1 = z;
					break;
				}
			}
			if (ptid >= 0) {
				if (sid >= 0) {
					SeedML& seed1 = seeds[sid];
					if (!merge_seed(seed, seed1)) break;
				} else {
					SliceML& pt1 = all_slices[z1][ptid];
					attach_seed_slice(seed, z1, ptid);
				}
				
			} else break;
			if (seed.bottom() == z0) break;
		}
		
		while (true) {
			int z0 = seed.top();
			SliceML& pt = all_slices[z0][seed.zmap[z0][0]];
			int sid=-1, ptid=-1, z1=-1;
			for (int z=z0+1; z<d; z++) {
				if (z - z0 > MAX_GAP) break;
				sid = find_seed_match(pt, z, &ptid);
				if (ptid >= 0) {
					z1 = z;
					break;
				}
			}
			if (ptid >= 0) {
				if (sid >= 0) {
					SeedML& seed1 = seeds[sid];
					int sid1 = seed1.ptid;
					if (!merge_seed(seed, seed1)) break;
				} else {
					SliceML& pt1 = all_slices[z1][ptid];
					attach_seed_slice(seed, z1, ptid);
				}
				
			} else break;
			if (seed.top() == z0) break;
		}
	}

	for (SeedML& seed : seeds) {
		if (seed.ptid < 0) continue;
		if (seed.top() - seed.bottom() + 1 >= MIN_HEIGHT) continue;
		for (int z=seed.bottom(); z<=seed.top(); z++) {
			if (seed.zmap[z].empty()) continue;
			SliceML& ptc = all_slices[z][seed.zmap[z][0]];
			ptc.ptid = -1;
			seed.zmap[z].clear();
		}
		seed.ptid = -1;
	}
}

void AssemblerML::cells_from_seeds()
{
	for (SeedML& seed : seeds) {
		if (seed.ptid < 0) continue;
		cells.resize(cells.size()+1);
		Particle3D& cell = cells[cells.size()-1];
		cell.fills.resize(d);
		for (int z=0; z<d; z++) {
			if (seed.zmap[z].empty()) continue;
			SliceML& ptc = all_slices[z][seed.zmap[z][0]];
			std::vector<HSeg>& my_fill = cell.fills[z];
			for (HSeg& hs : ptc.fill)
				my_fill.push_back(hs);
		}
		cell.update_from_fill();
	}
	
	sorted_cells.resize(cells.size());
	for (int idx=0; size_t(idx)<cells.size(); idx++) {
		Particle3D& cell = cells[idx];
		sorted_cells[idx].ptid = idx;
		sorted_cells[idx].sc = cell.iou_score(MAX_GAP);
	}
	std::sort(sorted_cells.begin(), sorted_cells.end(), compare_cell_sort);
}

static void frame_seq(std::vector<int> &res, int zlo, int zhi)
{
	int zmid = (zlo + zhi) / 2;
	if (zmid <= zlo || zmid >= zhi) return;
	res.push_back(zmid);
	if (zmid - zlo > 1) {
		frame_seq(res, zlo, zmid);
	}
	if (zhi - zmid > 1) {
		frame_seq(res, zmid, zhi);
	}
}

void AssemblerML::fill_gaps()
{
	mstack.fill(0);
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x20, 1);
	}
	for (SortScoreML& ss : sorted_cells) {
		Particle3D& cell = cells[ss.ptid];
		
		int zlo = -1, zhi = -1;
		for (int z=cell.bnd.zmin; z<=cell.bnd.zmax; z++) {
			std::vector<HSeg>& fill = cell.fills[z];
			if (fill.empty()) continue;
			Raster8 msk = mstack.getPlane(z);
			msk.paintParticleFillInto(fill, 0x60, 0);
			Boundary bnd = fill_boundary(fill);
			msk.rescanParticleFill(bnd, fill, 0x60);
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.expandBordersInto(bnd, 0x60, 0, 0x50, HOOD_SIZE_MOORE);
			msk.replaceColor(bnd, 0x60, 0x40);
			msk.paintParticleFill(fill, 0x80);
			
			if (zlo < 0) zlo = z;
			zhi = z;
			if (zhi - zlo > 1) {
				std::vector<int> mf;
				frame_seq(mf, zlo, zhi);
				for (int zmiss : mf) {
					interpolate_frame(cell, zmiss);
				}
			}
			zlo = zhi;
		}
		
	}
}

void AssemblerML::fix_borders_actin()
{
	mstack.fill(0);

	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x10, 1);
		int npt = 0;
		for (Particle3D& cell : cells) {
			if (cell.fills[z].empty()) continue;
			++npt;
			msk.paintParticleFill(cell.fills[z], 0xC0);
		}
		if (!npt) continue;
		
		msk.expandBorders(0xC0, 0x5, HOOD_SIZE_MOORE, 0);
		msk.expandBorders(0, 0x15, HOOD_SIZE_MOORE, 0x5);
		msk.replaceColor(0x15, 0);
		msk.replaceColor(0x5, 0x40);
		
		msk.expandBorders(0xC0, 0x40, HOOD_SIZE_FATCROSS, 0);
		msk.expandBorders(0x40, 0x30, HOOD_SIZE_FATCROSS, 0);
		msk.replaceColor(0x30, 0x40);
		msk.filterParticles(0, 0x5, 300, 0x80);
		msk.replaceColor(0x5, 0);
		msk.replaceColor(0x80, 0x40);
		msk.expandBorders(0, 0x60, HOOD_SIZE_FATCROSS, 0x40);
		msk.replaceColor(0x60, 0);
		msk.expandBorders(0, 0x60, HOOD_SIZE_NEUMANN, 0x40);
		msk.replaceColor(0x60, 0);

		for (int pass=0; pass<8; pass++) {
			int nbsz = (pass & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE;
			for (SortScoreML& ss : sorted_cells) {
				Particle3D& cell = cells[ss.ptid];
				std::vector<HSeg>& fill = cell.fills[z];
				if (fill.empty()) continue;
				Boundary bnd = fill_boundary(fill);
				bnd.expand(1);
				msk.clip(bnd, 1);
				msk.paintParticleFill(fill, 0xD0);
				msk.bordersAround(bnd, 0xD0, 0x40, 0xA0, nbsz);
				
				for (int y=bnd.ymin; y<=bnd.ymax; y++) {
					unsigned char *p = msk.scanLine(y);
					for (int x=bnd.xmin; x<=bnd.xmax; x++) {
						if (p[x] != 0xA0) continue;
						unsigned char xc = 0xD0;
						for (int j=1; j<HOOD_SIZE_MOORE; j++) {
							unsigned char c = msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy);
							if (c == 0xC0) {
								xc = 0;
								break;
							}
						}
						p[x] = xc;
					}
				}

				msk.rescanParticleFill(bnd, fill, 0xD0);
				msk.paintParticleFill(fill, 0xC0);
			}
		}
	}
}

void AssemblerML::fix_borders_dna()
{
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x10, 1);
		
		int npt = 0;
		// Filter out "phantom slices"
		for (Particle3D& cell : cells) {
			std::vector<HSeg>& fill = cell.fills[z];
			if (fill.empty()) continue;
			double ff = double(msk.countColor(fill, 0xFF)) / fill_area(fill);
			if (ff < DNA_PHANTOM_FILL) {
				fill.clear();
			} else {
				++npt;
				msk.paintParticleFill(fill, 0xC0);
			}
		}
		if (!npt) continue;
		
		for (int pass=0; pass<2; pass++) {
			int nbsz = (pass & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE;
			for (SortScoreML& ss : sorted_cells) {
				Particle3D& cell = cells[ss.ptid];
				std::vector<HSeg>& fill = cell.fills[z];
				if (fill.empty()) continue;
				Boundary bnd = fill_boundary(fill);
				bnd.expand(1);
				msk.clip(bnd, 1);
				msk.paintParticleFill(fill, 0xD0);
				msk.bordersAround(bnd, 0xD0, 0xFF, 0xA0, nbsz);
				
				for (int y=bnd.ymin; y<=bnd.ymax; y++) {
					unsigned char *p = msk.scanLine(y);
					for (int x=bnd.xmin; x<=bnd.xmax; x++) {
						if (p[x] != 0xA0) continue;
						unsigned char xc = 0xD0;
						for (int j=1; j<HOOD_SIZE_MOORE; j++) {
							unsigned char c = msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy);
							if (c == 0xC0) {
								xc = 0;
								break;
							}
						}
						p[x] = xc;
					}
				}

				msk.rescanParticleFill(bnd, fill, 0xD0);
				msk.paintParticleFill(fill, 0xC0);
			}
		}
	}

	for (Particle3D& cell : cells)
		cell.update_from_fill();
}

void AssemblerML::paint_final()
{
	Boundary3D bnd = mstack.getBoundary();
	mstack.fill(0x80);
	for (Particle3D& cell : cells) {
		if (cell.empty()) continue;
		mstack.paintParticle(cell, 0xE1);
		
		for (int z=0; size_t(z)<cell.fills.size(); z++) {
			if (cell.fills[z].empty()) continue;
			for (HSeg &hs : cell.fills[z]) {
				unsigned char *p = mstack.scanLine(hs.y, z);
				for (int x=hs.xl; x<=hs.xr; x++) {
					for (int j=1; j<HOOD3D_26; j++) {
						int x0 = x + hood3d_pts[j].dx;
						int y0 = hs.y + hood3d_pts[j].dy;
						int z0 = z + hood3d_pts[j].dz;
						if (!bnd.IsInside(x0, y0, z0)) continue;
						unsigned char v = mstack.value(x0, y0, z0);
						if (v != 0xE0 && v != 0xE1) {
							p[x] = 0xE0;
							break;
						}
					}
				}
			}
		}

		mstack.paintParticleADD(cell, 0x1F);
	}
}


// Protected

static void shrink_particle(Raster8& msk, Particle& pt, unsigned char fg, unsigned char bord)
{
	for (HSeg& hs: pt.fill) {
		unsigned char *p = msk.scanLine(hs.y);
		for (int x0 = hs.xl; x0 <= hs.xr; x0++) {
			for (int j=1; j<HOOD_SIZE_FATCROSS; j++) {
				unsigned char c = msk.value(x0+hood_pts[j].dx, hs.y+hood_pts[j].dy);
				if (c != fg && c != bord) {
					p[x0] = bord;
					break;
				}
			}
		}
	}
}

static void expand_particle(Raster8& msk, Particle& pt, unsigned char fg, unsigned char bord, unsigned char tmpc, int nbsz)
{
	for (HSeg& hs: pt.fill) {
		unsigned char *p = msk.scanLine(hs.y);
		for (int x0 = hs.xl; x0 <= hs.xr; x0++) {
			if (p[x0] != bord) continue;
			for (int j=1; j<nbsz; j++) {
				if (msk.value(x0+hood_pts[j].dx, hs.y+hood_pts[j].dy) == fg) {
					p[x0] = tmpc;
					break;
				}
			}
		}
	}
	msk.paintParticleInto(pt, fg, tmpc);
}

static void biggest_particle(Raster8& msk, Particle& ptc, unsigned char fg, unsigned char bord, unsigned char tmpc)
{
	std::vector<Slice> ptcs;
	int a = 0;
	int idx = -1;
	for (HSeg& hs: ptc.fill) {
		unsigned char *p = msk.scanLine(hs.y);
		for (int x0 = hs.xl; x0 <= hs.xr; x0++) {
			if (p[x0] != fg) continue;
			int _idx = int(ptcs.size());
			ptcs.resize(ptcs.size()+1);
			Slice& pt = ptcs[ptcs.size()-1];
			pt.area = msk.detectParticle(pt, x0, hs.y, tmpc);
			if (pt.area > a) {
				a = pt.area;
				idx = _idx;
			}
		}
	}
	for (int i=0; size_t(i)<ptcs.size(); i++) {
		msk.paintParticle(ptcs[i], (i == idx) ? fg : bord);
	}
}

void AssemblerML::smooth_slice(SliceML& pt, int z)
{
	Boundary bnd = pt.bnd;
	bnd = pt.bnd;
	bnd.expand(3);
	msk.clip(bnd, 2);
	msk.fill(bnd, 0);
	msk.paintParticle(pt, 0x80);

	shrink_particle(msk, pt, 0x80, 0x60);
	
	biggest_particle(msk, pt, 0x80, 0x60, 0x50);
	
	expand_particle(msk, pt, 0x80, 0x60, 0x50, HOOD_SIZE_NEUMANN);
	expand_particle(msk, pt, 0x80, 0x60, 0x50, HOOD_SIZE_MOORE);
	expand_particle(msk, pt, 0x80, 0x60, 0x50, HOOD_SIZE_NEUMANN);
	
	msk.replaceColor(0x60, 0);
	
	bnd.expand(-2);
	msk.expandBordersInto(bnd, 0x80, 0, 0x40, HOOD_SIZE_FATCROSS, true);
	msk.paintBoundary(bnd, 0);
	msk.fillParticle(bnd.xmin, bnd.ymin, 0x20);
	msk.replaceColor(bnd, 0, 0x80);
	msk.replaceColor(bnd, 0x20, 0);
	
	msk.expandBordersInto(bnd, 0, 0x40, 0x50, HOOD_SIZE_NEUMANN);
	msk.expandBordersInto(bnd, 0, 0x40, 0x50, HOOD_SIZE_MOORE);
	msk.expandBordersInto(bnd, 0, 0x40, 0x50, HOOD_SIZE_NEUMANN);
	msk.replaceColor(0x40, 0x80);
	
	msk.rescanParticleFill(bnd, pt.fill, 0x80);
	
	pt.area = pt.update_from_fill();
	bnd.expand(2);
	msk.fill(bnd, 0x10);
}

void AssemblerML::find_z_match(std::vector<int>& res, SliceML& pt, int z, double min_iou, bool any)
{
	res.resize(0);
	
	std::vector<SliceML> & slices = all_slices[z];
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		SliceML& ptc = slices[idx];
		if (!ptc.valid() || (!any && ptc.taken())) continue;
		if (!ptc.bnd.intersects(pt.bnd)) continue;
		int o = ptc.overlay_area(pt);
		if (o == 0) continue;
		double iou = double(o) / (ptc.area + pt.area - o);
		if (iou >= min_iou) res.push_back(idx);
	}
}

void AssemblerML::assign_ptids_to_seed(SeedML& seed, std::vector<int>& ptids, int z)
{
	std::vector<SliceML> & slices = all_slices[z];
	std::vector<int>& szids = seed.zmap[z];
	for (int idx : ptids) {
		SliceML& pt = slices[idx];
		szids.push_back(idx);
		pt.ptid = seed.ptid;
	}
}

int AssemblerML::assign_below(SeedML& seed, int z, double min_iou)
{
	if (z < 1) return 0;
	std::vector<SliceML> & slices = all_slices[z];
	if (slices.empty()) return 0;
	for (int idx : seed.zmap[z]) {
		SliceML& pt = slices[idx];
		std::vector<int> ptids;
		find_z_match(ptids, pt, z-1, min_iou, false);
		if (!ptids.empty())
			assign_ptids_to_seed(seed, ptids, z-1);
	}
	return int(seed.zmap[z-1].size());
}

void AssemblerML::dedupe_seed(SeedML& seed)
{
	int z;
	for (z=0; z<d; z++) {
		if (seed.zmap[z].empty()) continue;
		seed.lnkmap[z].resize(seed.zmap[z].size());
		for (int idx=0; size_t(idx)<seed.zmap[z].size(); idx++) {
			seed.lnkmap[z][idx] = best_lnk_at(seed, z, seed.zmap[z][idx]);
		}
	}
	
	int cur_idx = -1;
	int next_idx = -1;
	int next_j = -1;
	for (z=d-1; z>=0; z--) {
		if (seed.zmap[z].empty()) continue;
		double best_sc = 0.;
		for (int j=0; size_t(j) < seed.zmap[z].size(); j++) {
			if (cur_idx < 0 || best_sc < seed.lnkmap[z][j].sc) {
				best_sc = seed.lnkmap[z][j].sc;
				cur_idx = seed.zmap[z][j];
				next_j = seed.lnkmap[z][j].j;
				next_idx = seed.lnkmap[z][j].ptid;
			}
		}
		seed.sc = best_sc;
		break;
	}
	
	for (; z>=0; z--) {
		if (seed.zmap[z].empty()) break;
		if (cur_idx >= 0) {
			seed.zmap[z].resize(1);
			seed.zmap[z][0] = cur_idx;
		} else break;
		if (z == 0) break;
		cur_idx = next_idx;
		if (next_j >= 0) {
			next_idx = seed.lnkmap[z-1][next_j].ptid;
			next_j = seed.lnkmap[z-1][next_j].j;
		} else {
			next_j = next_idx = -1;
		}
	}
}

SortScoreML AssemblerML::best_lnk_at(SeedML& seed, int z, int idx)
{
	//std::cout << "best_lnk_at: " << seed.ptid << " Z=" << z << " idx=" << idx << std::endl;
	SortScoreML ss;
	if (z < 1) return ss;
	int zz = z - 1;
	if (seed.zmap[zz].empty()) return ss;
	SliceML& ptc = all_slices[z][idx];
	for (int j=0; size_t(j)<seed.zmap[zz].size(); j++) {
		int jdx = seed.zmap[zz][j];
		SliceML& pt = all_slices[zz][jdx];
		SortScoreML& jss = seed.lnkmap[zz][j];
		int o = ptc.overlay_area(pt);
		double sc = double(o) / double(pt.area + ptc.area - o) + jss.sc;
		if (ss.ptid < 0 || ss.sc < sc) {
			ss.sc = sc;
			ss.ptid = jdx;
			ss.j = j;
		}
	}
	return ss;
}

int AssemblerML::find_seed_match(SliceML& pt, int z, int *pptid)
{
	std::vector<int> ptids;
	find_z_match(ptids, pt, z, ACCEPTABLE_IOU, true);
	int b_ptid = -1;
	double b_sc = 0.;
	for (int ptid : ptids) {
		SliceML& ptc = all_slices[z][ptid];
		if (ptc.ptid >= 0) {
			*pptid = ptid;
			return ptc.ptid;
		}
		int o = ptc.overlay_area(pt);
		double iou = double(o) / (ptc.area + pt.area - o);
		if (b_ptid < 0 || b_sc < iou) {
			b_ptid = ptid;
			b_sc = iou;
		}
	}
	*pptid = b_ptid;
	return -1;
}

void AssemblerML::attach_seed_slice(SeedML& seed, int z, int idx)
{
	SliceML& pt0 = all_slices[z][idx];
	pt0.ptid = seed.ptid;
	pt0.sc = -1.;
	seed.zmap[z].resize(1);
	seed.zmap[z][0] = idx;
	
	Boundary bnd = pt0.bnd;
	msk.paintParticle(pt0, 0x40);
	int ptcnt = 0;
	for (int dz=1; dz<10; dz++) {
		int za = z + dz;
		if (za < d) {
			if (!seed.zmap[za].empty()) {
				SliceML& pt = all_slices[za][seed.zmap[za][0]];
				msk.paintParticleInto(pt, 0x80, 0x40);
				++ptcnt;
			}
		}
		int zb = z - dz;
		if (zb >= 0) {
			if (!seed.zmap[zb].empty()) {
				SliceML& pt = all_slices[zb][seed.zmap[zb][0]];
				msk.paintParticleInto(pt, 0x80, 0x40);
				++ptcnt;
			}
		}
		if (ptcnt >= 2) break;
	}
	bnd.expand(1);
	msk.clip(bnd, 1);
	msk.expandBordersInto(bnd, 0x40, 0x80, 0x50, HOOD_SIZE_MOORE, false);
	msk.rescanParticleFill(bnd, pt0.fill, 0x80);
	pt0.area = pt0.update_from_fill();
	msk.fill(bnd, 0x10);
	smooth_slice(pt0, z);
}

bool AssemblerML::merge_seed(SeedML& seed, SeedML& seed1)
{
	for (int z=0; z<d; z++) {
		if (!seed1.zmap[z].empty() && !seed.zmap[z].empty())
			return false;
	}
	for (int z=0; z<d; z++) {
		if (seed1.zmap[z].empty()) continue;
		int idx = seed1.zmap[z][0];
		seed1.zmap[z].clear();
		SliceML& pt = all_slices[z][idx];
		seed.zmap[z].push_back(idx);
		pt.ptid = seed.ptid;
		pt.sc = -2.;
	}
	seed1.ptid = -1;
	seed1.sc = 0.;
	return true;
}

void AssemblerML::interpolate_frame(Particle3D& cell, int zmiss)
{
	Raster8 msk = mstack.getPlane(zmiss);
	int zlo, zhi;
	for (zlo = zmiss-1; zlo>=0; zlo--) {
		if (!cell.fills[zlo].empty()) break;
	}
	for (zhi = zmiss+1; zhi<d; zhi++) {
		if (!cell.fills[zhi].empty()) break;
	}
	if (zlo < 0 || zhi >= d) return;
	
	std::vector<HSeg>& fill_lo = cell.fills[zlo];
	std::vector<HSeg>& fill = cell.fills[zmiss];
	std::vector<HSeg>& fill_hi = cell.fills[zhi];
	Boundary bnd = fill_boundary(fill_lo);
	Boundary bnd_hi = fill_boundary(fill_hi);
	bnd.combo(bnd_hi);
	
	msk.paintParticleFillInto(fill_lo, 0x44, 0);
	msk.paintParticleFillInto(fill_hi, 0x60, 0x44);
	msk.paintParticleFillInto(fill_hi, 0x44, 0);
	
	for (int y0=bnd.ymin+1; y0<=bnd.ymax-1; y0++) {
		unsigned char *b = msk.scanLine(y0);
		for (int x0=bnd.xmin+1; x0<=bnd.xmax-1; x0++) {
			if (b[x0] != 0x44) continue;
			int nc = 0, nb = 0;
			for (int j=1; j<HOOD_SIZE_FATCROSS; j++) {
				unsigned char c = msk.value(x0+hood_pts[j].dx, y0+hood_pts[j].dy);
				if (c == 0x60) {
					++nc;
					break;
				}
				if (c != 0x44 && c != 0x48 && c != 0x30) {
					++nb;
				}
			}
			b[x0] = (nc > 0 || nb == 0) ? 0x48 : 0x30;
		}
	}
	msk.replaceColor(0x30, 0);
	msk.replaceColor(0x48, 0x60);
	
	msk.rescanParticleFill(bnd, fill, 0x60);
	bnd.expand(1);
	msk.clip(bnd, 1);
	msk.expandBordersInto(bnd, 0x60, 0, 0x50, HOOD_SIZE_MOORE);
	msk.replaceColor(bnd, 0x60, 0x40);
	msk.paintParticleFill(fill, 0xFF);

}





