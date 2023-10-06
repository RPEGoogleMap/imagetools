
#include "geom.h"
#include "raster.h"
#include "ml3d.h"
#include "comparator.h"
#include "multifit.h"
#define __MAIN__
#include "imagetools.h"

#include "version.h"

const char *version() { return __version__; }

static void paint_border_particles(Raster8& msk, unsigned char c)
{
	int x0, y0;
	unsigned char *p;
	
	y0 = 1;
	p = msk.scanLine(y0);
	for (x0=1; x0<msk.w-1; x0++) {
		if (p[x0] == 0) {
			Particle pt;
			msk.detectParticle(pt, x0, y0, c);
		}
	}
	y0 = msk.h - 2;
	p = msk.scanLine(y0);
	for (x0=1; x0<msk.w-1; x0++) {
		if (p[x0] == 0) {
			Particle pt;
			msk.detectParticle(pt, x0, y0, c);
		}
	}
	for (y0=1; y0<msk.h-1; y0++) {
		x0 = 1;
		if (msk.value(x0, y0) == 0) {
			Particle pt;
			msk.detectParticle(pt, x0, y0, c);
		}
		x0 = msk.w - 2;
		if (msk.value(x0, y0) == 0) {
			Particle pt;
			msk.detectParticle(pt, x0, y0, c);
		}
	}
}

template <class T>
static std::vector<T> detect_particles_2d(Raster8& msk)
{
	// Exclude pixels at the border
	msk.fillBorder(0x10, 1);
	paint_border_particles(msk, 0x10);
	
	std::vector<T> particles;
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *p = msk.scanLine(y0);
		for (int x0=1; x0<msk.w-1; x0++) {
			if (p[x0] == 0) {
				particles.resize(particles.size()+1);
				T& ptc = particles[particles.size()-1];
				ptc.area = msk.detectParticle(ptc, x0, y0, 0x80);
			}
		}
	}
	
	return particles;
}
		
int detect_particles_csv(unsigned char *mask, int hm, int wm, const char *csvfile)
{
	Raster8 msk(wm, hm, mask);
	std::vector<Slice> particles = detect_particles_2d<Slice>(msk);
	write_particle_data(particles, csvfile);
	
	return int(particles.size());
}

// std::vector<std::vector<int>> res;
std::vector<std::vector<int>> detect_particles(unsigned char *mask, int hm, int wm)
{
	Raster8 msk(wm, hm, mask);
	std::vector<Slice> particles = detect_particles_2d<Slice>(msk);
	std::vector<std::vector<int>> res(particles.size());
	
	for (size_t i=0; i<particles.size(); i++) {
		Slice& ptc = particles[i];
		std::vector<int>& ptfill = res[i];
		for (HSeg & hs : ptc.fill) {
			ptfill.push_back(hs.y);
			ptfill.push_back(hs.xl);
			ptfill.push_back(hs.xr);
		}
	}

	return res;
}

void postprocess_particle_borders(unsigned char *mask, int hm, int wm, bool expand)
{
	Raster8 out(wm, hm, mask);

	out.fillBorder(0x10, 1);
	out.expandBorders(0xFF, 0x80, HOOD_SIZE_FATCROSS, 0);
	out.expandBorders(0x80, 0x78, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 20, 0x30);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x30, 0x80);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x78, 0x70, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x70, 0x68, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x68, 0x60, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x60, 0x58, HOOD_SIZE_FATCROSS, 0);

	out.replaceColor(0x40, 0);
	out.replaceColor(0xC0, 0);
	
	out.filterParticles(0, 0x20, 50, 0x80);
	out.replaceColor(0x20, 0);
	
	// Try to rejoin small chips to bigger ones nearby
	out.filterParticles(0, 0x20, 200, 0xA0);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x58, 0x60);
	for (int pass=0; pass<5; pass++) {
		out.expandBorders(0xA0, 0x30, HOOD_SIZE_FATCROSS, 0x60);
		out.replaceColor(0x30, 0xA0);
	}
	out.replaceColor(0xA0, 0);
	
	std::vector<Slice> particles;
	for (int y0=1; y0<out.h-1; y0++) {
		unsigned char *p = out.scanLine(y0);
		for (int x0=1; x0<out.w-1; x0++) {
			if (p[x0] != 0) continue;
			particles.resize(particles.size()+1);
			Slice &ptc = particles[particles.size()-1];
			ptc.area = out.detectParticle(ptc, x0, y0, 0x30);
		}
	}
	for (Slice &ptc : particles)
		out.paintParticle(ptc, 0);
	int npass = 8;
	for (int lc=0x60; lc<=0x80; lc+=0x8) {
		unsigned char curc = (unsigned char)(lc);
		unsigned char nextc = curc + 0x8;
		for (int pass=0; pass<npass; pass++) {
			int nbsz = (pass & 1) ? HOOD_SIZE_MOORE : HOOD_SIZE_NEUMANN;
			for (Slice &ptc : particles)
				out.expandParticle(ptc, 0, 0x30, curc, nbsz);
		}
		if (curc == 0x80) break;
		out.replaceColor(curc, nextc);
		npass = 4;
	}
	
	// Perform a few more erosion iterations to cut into deep corners
	out.replaceColor(0xFF, 0x80);
	for (int pass=0; pass<8; pass++) {
		for (Slice &ptc : particles)
			out.expandParticle(ptc, 0, 0x30, 0x80, HOOD_SIZE_NEUMANN);
	}

	out.replaceColor(out.forcedc, 0x80);
	out.expandBorders(0, 0xFF, HOOD_SIZE_MOORE, 0x80);
	
	if (expand) {
		out.expandBorders(0xFF, 0xC0, HOOD_SIZE_MOORE, 0);
		out.replaceColor(0xC0, 0xFF);
	}

	out.fillBorder(0, 1);
}


// -- called from masks_to_particles() --
// Process single detected mask; if it consists of several disconnected areas, paint small ones (<minarea)
// with out_fgd color, the rest -- with in_fgd.
// Return bounding box containing all areas (big and small), expanded by 1 pixel in each direction.
static Boundary filterParticlesIn(Raster8& msk, int minarea, unsigned char fgd, unsigned char in_fgd, unsigned char out_fgd)
{
	Boundary bnd(msk.w, msk.h, -1, -1);
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *p = msk.scanLine(y0);
		for (int x0=0; x0<msk.w-1; x0++) {
			if (p[x0] != fgd) continue;
			std::vector<HSeg> fill;
			msk.findParticleFill(fill, x0, y0, in_fgd);
			if (fill_area(fill) < minarea)
				msk.paintParticleFill(fill, out_fgd);
			Boundary ptbnd = fill_boundary(fill);
			bnd.combo(ptbnd);
		}
	}
	bnd.expand(1);
	msk.clip(bnd, 1);
	return bnd;
}

std::vector<std::vector<int>> masks_to_particles(unsigned char *ptmask, int npts, int hptm, int wptm,
		int x_orig, int y_orig, int minarea)
{
	std::vector<std::vector<int>> res;
	
	//std::cout << "npts=" << npts << " w=" << wptm << " h=" << hptm << std::endl;
	
	size_t fr_size = size_t(hptm) * size_t(wptm);
	for (int ipt=0; ipt<npts; ipt++) {
		Raster8 msk(wptm, hptm, ptmask + (ipt * fr_size));
		Boundary bnd = filterParticlesIn(msk, minarea, 0xFF, 0xA0, 0x80);
		// Smooth masked object boundaries a little
		msk.bordersAround(bnd, 0xA0, 0, 0x80, HOOD_SIZE_NEUMANN);
		msk.bordersAround(bnd, 0, 0x80, 0x40, HOOD_SIZE_MOORE);
		msk.replaceColor(bnd, 0x40, 0);
		msk.replaceColor(bnd, 0x80, 0xA0);
		Slice ptc;
		ptc.bnd = bnd;
		ptc.area = msk.rescanParticle(ptc, 0xA0);
		if (ptc.area == 0) continue;	// if the object disappeared after filtering/smoothing
		res.resize(res.size() + 1);
		std::vector<int>& ptfill = res[res.size()-1];
		for (HSeg & hs : ptc.fill) {
			ptfill.push_back(hs.y + y_orig);
			ptfill.push_back(hs.xl + x_orig);
			ptfill.push_back(hs.xr + x_orig);
		}
	}
	
	//std::cout << "Found " << res.size() << " masks." << std::endl;

	return res;
}

void assemble_ml(std::vector<std::vector<std::vector<int>>> particles_3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *csvfile, int postproc,
		double good_iou, double ok_iou)
{
	AssemblerML aml(mask3d, zm3d, hm3d, wm3d);
	aml.GOOD_IOU = good_iou;
	aml.ACCEPTABLE_IOU = ok_iou;
	
	aml.load_particles_3d(particles_3d);
	
	aml.find_seeds();
	aml.consolidate_seeds();
	aml.cells_from_seeds();
	
	if (postproc & POSTPROC_DNA) {
		// std::cout << "fix_borders_dna" << std::endl;
		aml.fix_borders_dna();
	} else if (postproc & POSTPROC_ACTIN) {
		// std::cout << "fix_borders_actin" << std::endl;
		aml.fix_borders_actin();
	}
	aml.fill_gaps();
	
	if (postproc & POSTPROC_SANDPAPER) {
		aml.mstack.sandPaperCells(aml.cells);
	}
	
	aml.paint_final();
	
	std::cout << "Write " << aml.cells.size() << " objects to " << csvfile << std::endl;
	write_cell_data(aml.cells, csvfile);
}

void assemble_2d(std::vector<std::vector<int>> particles_2d, std::vector<double> scores,
		unsigned short *data, int hd, int wd, const char *csvfile, int postproc)
{
	Raster16 dat(wd, hd, data);
	Raster8 msk(wd, hd, NULL);
	
	std::vector<Slice> slices(particles_2d.size());
	for (int idx=0; size_t(idx)<particles_2d.size(); idx++) {
		std::vector<int>& pt2d = particles_2d[idx];
		Slice& ptc = slices[idx];
		for (size_t i=0; i<pt2d.size(); i+=3) {
			if (i+2 >= pt2d.size()) break;
			ptc.fill.push_back(HSeg(pt2d[i], pt2d[i+1], pt2d[i+2]));
		}
		ptc.area = ptc.update_from_fill();
		
		for (int j=0; j<idx; j++) {
			Slice& ptc0 = slices[j];
			if (ptc0.fill.empty()) continue;
			if (!ptc0.bnd.intersects(ptc.bnd)) continue;
			int ovl = ptc0.overlay_area(ptc);
			double iou = double(ovl) / double(ptc.area + ptc0.area - ovl);
			if (iou >= 0.25) {
				if (scores[j] >= scores[idx]) {
					ptc.clear();
					break;
				} else {
					ptc0.clear();
				}
			}
		}
	}
	
	msk.fill(0);
	msk.fillBorder(0x10);
	for (Slice& ptc : slices) {
		if (ptc.fill.empty()) continue;
		msk.paintParticleInto(ptc, 0x60, 0);
		ptc.area = msk.rescanParticle(ptc, 0x60);
		msk.paintParticle(ptc, 0x80);
	}
	
	// Post-processing
	if (postproc & POSTPROC_DNA) {
		for (long long i=0; i<msk.len; i++) {
			if (msk.buf[i] == 0 && dat.buf[i] != 0)
				msk.buf[i] = 0x40;
		}
		for (Slice& ptc : slices) {
			if (ptc.fill.empty()) continue;
			msk.paintParticle(ptc, 0xC0);
			ptc.bnd.expand(2);
			msk.clip(ptc.bnd, 1);
			msk.expandBordersInto(ptc.bnd, 0xC0, 0x40, 0x55, HOOD_SIZE_MOORE);
			msk.expandBordersInto(ptc.bnd, 0xC0, 0x40, 0x55, HOOD_SIZE_NEUMANN);
			ptc.area = msk.rescanParticle(ptc, 0xC0);
			msk.paintParticle(ptc, 0x80);
		}
		msk.replaceColor(0x40, 0);
	} else if (postproc & POSTPROC_ACTIN) {
		for (Slice& ptc : slices) {
			if (ptc.fill.empty()) continue;
			msk.paintParticle(ptc, 0xC0);
			ptc.bnd.expand(2);
			msk.clip(ptc.bnd, 1);
			msk.expandBordersInto(ptc.bnd, 0xC0, 0, 0x55, HOOD_SIZE_MOORE);
			msk.expandBordersInto(ptc.bnd, 0xC0, 0, 0x55, HOOD_SIZE_NEUMANN);
			ptc.area = msk.rescanParticle(ptc, 0xC0);
			msk.paintParticle(ptc, 0x80);
		}
	}
	
	// Enforce 1-pixel separation between particles
	for (Slice& ptc : slices) {
		if (ptc.fill.empty()) continue;
		msk.paintParticle(ptc, 0xC0);
		for (HSeg& hs : ptc.fill) {
			unsigned char *p = msk.scanLine(hs.y);
			for (int x=hs.xl; x<=hs.xr; x++) {
				for (int j=1; j<HOOD_SIZE_MOORE; j++) {
					if (msk.value(x+hood_pts[j].dx, hs.y+hood_pts[j].dy) == 0x80) {
						p[x] = 0x40;
						break;
					}
				}
				
			}
		}
		ptc.area = msk.rescanParticle(ptc, 0xC0);
		msk.paintParticle(ptc, 0x80);
	}
	
	// Output
	dat.fill(0);
	int id = 1;
	for (Slice& ptc : slices) {
		if (ptc.fill.empty()) continue;
		dat.paintParticle(ptc, (unsigned short)id);
		++id;
	}
	std::cout << "Write " << id << " particles to " << csvfile << std::endl;
	write_particle_data(slices, csvfile);
}

Compare3dResult compare_3d_annotations(int w, int h, int d, const char *base_csv, const char *cmp_csv)
{
	Compare3dResult res;
	
	SegmentationComparator3D cmp(w, h, d, base_csv, cmp_csv);
	
	res.base_slices = cmp.num_base_slices;
	res.cmp_slices = cmp.num_cmp_slices;
	
	cmp.compare_slices();
	
	res.f_pos = cmp.f_pos;
	res.f_neg = cmp.f_neg;
	res.fragm = cmp.fragm;
	res.fused = cmp.fused;
	for (int i=0; i<100; i++) {
		res.pct_match.push_back(cmp.pct_match[i]);
	}
	
	cmp.compare_cells();

	res.base_cells = cmp.num_valid_base_cells;
	res.cmp_cells = cmp.num_valid_cmp_cells;

	res.f_pos_3d = cmp.f_pos_3d;
	res.f_neg_3d = cmp.f_neg_3d;
	res.fragm_3d = cmp.fragm_3d;
	res.fused_3d = cmp.fused_3d;
	for (int i=0; i<100; i++) {
		res.pct_match_3d.push_back(cmp.pct_match_3d[i]);
	}
	
	return res;
}

std::vector<std::vector<std::vector<int>>> border_pixels(int w, int h, int d, const char *csvfile)
{
	std::vector<std::vector<std::vector<int>>> res;
	
	std::vector<Particle3D> cells;
	if (read_cell_data(csvfile, cells, w, h, d) < 0)
		return res;
	if (cells.empty())
		return res;

	Raster3D mstack(w, h, d, NULL);
	mstack.fill(0x80);
	for (Particle3D& cell : cells)
		mstack.paintParticle(cell, 0);
	mstack.expandBorders(0x80, 0, 0xFF);
	
	res.resize(cells.size());
	for (int i=0; size_t(i)<cells.size(); i++) {
		Particle3D& cell = cells[i];
		std::vector<std::vector<int>>& cres = res[i];
		for (int z=0; z<size_t(cell.fills.size()); z++) {
			for (HSeg& hs : cell.fills[z]) {
				for (int x=hs.xl; x<=hs.xr; x++) {
					if (mstack.value(x, hs.y, z) != 0xFF) continue;
					cres.resize(cres.size()+1);
					std::vector<int>& cpt = cres[cres.size() - 1];
					cpt.push_back(x);
					cpt.push_back(hs.y);
					cpt.push_back(z);
				}
			}
		}
	}
	
	return res;
}

static void extrapolate_particle_3d(Raster3D& mstack, Particle3D& pt, int z0, int z1)
{
	Raster8 msk = mstack.getPlane(z1);
	Boundary bnd = pt.bnd.boundary2d();
	bnd.expand(1);
	msk.clip(bnd, 1);
	msk.paintParticleFillInto(pt.fills[z0], 0x50, 0);
	msk.paintParticleFillInto(pt.fills[z0], 0x50, 0xFF);
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (msk.value(x, y) != 0x50) continue;
			for (int j=1; j<HOOD_SIZE_MOORE; j++) {
				unsigned char c = msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy);
				if (c != 0x50 && c != 0x51) {
					msk.setValue(x, y, 0x51);
					break;
				}
			}
		}
	}
	msk.replaceColor(0x51, 0);
	msk.rescanParticleFill(bnd, pt.fills[z1], 0x50);
}

static void expand_borders_3d(
	Raster3D& mstack,
	unsigned char oldc, unsigned char newc, unsigned char bordc,
	int nbsz=HOOD3D_18)
{
	for (int z=0; z<mstack.d; z++) {
		for (int y=1; y<mstack.h-1; y++) {
			for (int x=1; x<mstack.w-1; x++) {
				if (mstack.value(x, y, z) != oldc) continue;
				
				for (int j=1; j<nbsz; j++) {
					int z0 = z + hood3d_pts[j].dz;
					if (z0 < 0 || z0 >= mstack.d) continue;
					int y0 = y + hood3d_pts[j].dy;
					int x0 = x + hood3d_pts[j].dx;
					if (mstack.value(x0, y0, z0) == bordc) {
						mstack.setValue(x, y, z, newc);
						break;
					}
				}
			}
		}
	}
}

void rs_tops_bottoms(unsigned char *mask3d, int zm3d, int hm3d, int wm3d, int xcv)
{
	Raster3D src_mstack(wm3d, hm3d, zm3d, mask3d);
	AssemblerML aml(NULL, zm3d, hm3d, wm3d);
	aml.GOOD_IOU = 0.75;
	aml.ACCEPTABLE_IOU = 0.6;
	
	for (int z=0; z<src_mstack.d; z++) {
		Raster8 msk = src_mstack.getPlane(z);
		aml.all_slices[z] = detect_particles_2d<SliceML>(msk);
		src_mstack.replaceColorInv(0xFF, 0);
	}
	
	aml.find_seeds();
	aml.consolidate_seeds();
	aml.cells_from_seeds();
	aml.fill_gaps();
	
	memcpy(aml.mstack.buf, src_mstack.buf, aml.mstack.len);
	
//	int painted = 0;
	for (Particle3D& pt : aml.cells) {
		if (pt.realHeight() > 3) {
			src_mstack.paintParticle(pt, 0x40);
//			++painted;
		} else {
			pt.clear();
		}
	}
//	std::cout << "Painted " << painted << " of " << aml.cells.size() << std::endl;
	
	Particle3D xpt;
	xpt.fills.resize(src_mstack.d);
	
	for (size_t idx=0; idx<aml.cells.size(); idx++) {
		Particle3D& pt = aml.cells[idx];
		src_mstack.paintParticle(pt, 0x50);
		
		int zl=-1, zh;
		for (int z=0; size_t(z) < pt.fills.size(); z++) {
			if (pt.fills[z].empty()) continue;
			if (zl < 0) zl = z;
			zh = z;
		}
		if (zl < 0) continue;
		
		if (zl == 0) {
			xpt.add_fill(pt.fills[zl], zl);
		}
		if (zh+1 == src_mstack.d) {
			xpt.add_fill(pt.fills[zh], zh);
		}
		if (zl > 0) {
			extrapolate_particle_3d(src_mstack, pt, zl, zl-1);
			--zl;
			if (zl == 0)
				xpt.add_fill(pt.fills[zl], zl);
		}
		if (zl > 0) {
			extrapolate_particle_3d(src_mstack, pt, zl, zl-1);
			--zl;
			xpt.add_fill(pt.fills[zl], zl);
		}
		if (zh+1 < src_mstack.d) {
			extrapolate_particle_3d(src_mstack, pt, zh, zh+1);
			++zh;
			xpt.add_fill(pt.fills[zh], zh);
		}
		
		src_mstack.paintParticle(pt, 0x40);
	}
	
	expand_borders_3d(src_mstack, 0, 0x80, 0x40, HOOD3D_18);

	for (int z=0; z<src_mstack.d; z++) {
		Raster8 msk = src_mstack.getPlane(z);
		msk.expandBorders(0x80, 0xA0, HOOD_SIZE_MOORE, 0);
		msk.replaceColor(0x80, 0xFF);
		msk.replaceColor(0xA0, 0xFF);
	}
	
	src_mstack.replaceColorInv(0xFF, 0);

	if (xcv != 0xFF) {
		src_mstack.paintParticle(xpt, 0x40);
		expand_borders_3d(src_mstack, 0xFF, 0x7F, 0x40, HOOD3D_18);

		for (int z=0; z<src_mstack.d; z++) {
			Raster8 msk = src_mstack.getPlane(z);
			msk.expandBorders(0x7F, 0xA0, HOOD_SIZE_FATCROSS, 0xFF);
			msk.replaceColor(0xA0, 0x7F);
		}
		
		src_mstack.replaceColor(0x40, 0);
		if (xcv != 0x7F)
			src_mstack.replaceColor(0x7F, (unsigned char)xcv);
	}
}

void sandpaper_cells(int w, int h, int d, const char *in_csv, const char *out_csv)
{
	std::vector<Particle3D> cells;
	read_cell_data(in_csv, cells, w, h, d);
	
	Raster3D mstack(w, h, d, NULL);
	mstack.sandPaperCells(cells);

	std::cout << "Write " << cells.size() << " objects to " << out_csv << std::endl;
	write_cell_data(cells, out_csv);
}

void paint_cells(unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *in_csv)
{
	std::vector<Particle3D> cells;
	read_cell_data(in_csv, cells, wm3d, hm3d, zm3d);
	
	Raster3D mstack(wm3d, hm3d, zm3d, mask3d);
	mstack.fill(0x80);
	
	for (Particle3D& cell : cells)
		mstack.paintParticle(cell, 0);
	mstack.expandBorders(0x80, 0, 0xFF, HOOD3D_26);
}

void filter_particles(unsigned char *mask, int hm, int wm, int minarea)
{
	Raster8 msk(wm, hm, mask);
	msk.fillBorder(0x40);
	msk.filterParticles(0xFF, 0xC0, minarea, 0);
	msk.replaceColor(0xC0, 0xFF);
	msk.fillBorder(0);
}

/*

class Raster16SampleIterator2D : public SampleIterator2D
{
protected:
	Raster16& dat;
	Boundary bnd;
	int x=0, y=0;
public:
	Raster16SampleIterator2D(Raster16& _dat, Boundary _bnd) : dat(_dat), bnd(_bnd) {}
	virtual int GetNumSamples() override
	{
		return int(bnd.area());
	}
	virtual void reset() override
	{
		y = bnd.ymin;
		x = bnd.xmin;
	}
	virtual bool next(double *px, double *py, double *pz) override
	{
		if (!bnd.IsInside(x, y)) return false;
		*px = double(x);
		*py = double(y);
		*pz = double(dat.value(x,y));
		++x;
		if (x >= bnd.xmax) {
			x = bnd.xmin;
			++y;
		}
		return true;
	}
};

class Raster16SampleIterator2D_p8 : public SampleIterator2D
{
protected:
	Raster16& dat;
	Boundary bnd;
	int x=0, y=0;
public:
	Raster16SampleIterator2D_p8(Raster16& _dat, Boundary _bnd) : dat(_dat), bnd(_bnd) {}
	virtual int GetNumSamples() override
	{
		return (bnd.xmax-bnd.xmin+1)*(bnd.ymax-bnd.ymin+1)/64;
	}
	virtual void reset() override
	{
		y = bnd.ymin;
		x = bnd.xmin;
	}
	virtual bool next(double *px, double *py, double *pz) override
	{
		if (x < bnd.xmin || x+7 > bnd.xmax || y < bnd.ymin || y+7 > bnd.ymax)
			return false;
		
		double v = 0.;
		for (int j=0; j<8; j++) {
			unsigned short *buf = dat.scanLine(y+j) + x;
			for (int i=0; i<8; i++)
				v += double(buf[i]);
		}
		*px = x + 4.;
		*py = y + 4.;
		*pz = v / 64.;
		
		x += 8;
		if (x+7 > bnd.xmax) {
			x = bnd.xmin;
			y += 8;
		}
		return true;
	}
};

double polyfit(unsigned short *data, int hd, int wd, int ord, int tilesz, int overlap)
{
	Raster16 dat(wd, hd, data);
	
	std::vector<PolynomialFit2D *> fits;
	std::vector<Boundary> bnds;
	
	for (int y=0; y<dat.h; y+=tilesz) {
		int ymin = y;
		int ymax = y + tilesz - 1;
		if (ymax >= dat.h) {
			ymax = dat.h - 1;
			ymin = ymax - tilesz + 1;
			if (ymin < 0) ymin = 0;
		}
		std::cout << "Fit: " << ymin << " till " << dat.h << std::endl;
		
		for (int x=0; x<dat.w; x+=tilesz) {
			int xmin = x;
			int xmax = x + tilesz - 1;
			if (xmax >= dat.w) {
				xmax = dat.w - 1;
				xmin = xmax - tilesz + 1;
				if (xmin < 0) xmin = 0;
			}
			Boundary inner(xmin, ymin, xmax, ymax);
			Boundary outer(xmin-overlap,ymin-overlap, xmax+overlap,ymax+overlap);
			dat.clip(outer);
			
			PolynomialFit2D* fit = new PolynomialFit2D(ord);
			Raster16SampleIterator2D_p8 iter(dat, outer);
			fit->fit(iter);
			fits.push_back(fit);
			bnds.push_back(inner);
		}
	}
	
	std::cout << "Predict " << fits.size() << std::endl;
	double dev = 0.;
	for (size_t j=0; j<fits.size(); j++) {
		PolynomialFit2D* fit = fits[j];
		Boundary bnd = bnds[j];
		
		for (int y=bnd.ymin; y<=bnd.ymax; y++) {
			unsigned short *buf = dat.scanLine(y);
			for (int x=bnd.xmin; x<=bnd.xmax; x++) {
				double v = fit->value(double(x), double(y));
				if (v < 0.) v = 0.;
				else if (v > 65535.) v = 65535.;
				double dv = v - buf[x];
				dev += dv*dv;
				buf[x] = (unsigned short)v;
			}
		}
		
		delete fit;
	}

	std::cout << "All done." << std::endl;
	
	return sqrt(dev / dat.len);
}
*/

void artimask(unsigned short *data, int hd, int wd, unsigned char *mask, int hm, int wm, double cutoff)
{
	Raster16 dat(wd, hd, data);
	Raster8 msk(wm, hm, mask);
	
	int h8 = hd / 8;
	int w8 = wd / 8;
	
	Raster16 hst(w8, h8, NULL);
	hst.fill(0);
	for (int y=0; y<h8; y++) {
		unsigned short *acc = hst.scanLine(y);
		for (int x=0; x<w8; x++) {
			double v = 0.;
			for (int j=0; j<8; j++) {
				unsigned short *buf = dat.scanLine(y*8+j) + (x*8);
				for (int i=0; i<8; i++) v += double(buf[i]);
			}
			acc[x] = (unsigned short)(v/64.);
		}
	}
	double stdev;
	double mean = hst.mean_std(&stdev);
	//std::cout << "mean=" << mean << " stdev=" << stdev << std::endl;
	
	unsigned short cut = (unsigned short)(mean + stdev*cutoff + 1.);
	//std::cout << "cut=" << cut << std::endl;

	for (long long i=0; i<msk.len; i++) {
		msk.buf[i] = (dat.buf[i] > cut) ? 0xFF : 0;
	}
	
	msk.expandBorders(0, 0x40, HOOD_SIZE_FATCROSS, 0xFF);
	msk.replaceColor(0x40, 0);
	
	msk.filterParticles(0xFF, 0xFC, 500, 0x80);
	msk.replaceColor(0x80, 0);
	msk.replaceColor(0xFC, 0xFF);
	
	msk.expandBorders(0xFF, 0x80, HOOD_SIZE_RAD3, 0);
	msk.replaceColor(0x80, 0xFF);

	msk.expandBorders(0xFF, 0x80, HOOD_SIZE_RAD3, 0);
	msk.replaceColor(0x80, 0xFF);

}


