
#include "geom.h"
#include "raster.h"
#include "ml3d.h"
#include "comparator.h"
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

static std::vector<Slice> detect_particles_2d(Raster8& msk)
{
	// Exclude pixels at the border
	msk.fillBorder(0x10, 1);
	paint_border_particles(msk, 0x10);
	
	std::vector<Slice> particles;
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *p = msk.scanLine(y0);
		for (int x0=1; x0<msk.w-1; x0++) {
			if (p[x0] == 0) {
				particles.resize(particles.size()+1);
				Slice& ptc = particles[particles.size()-1];
				ptc.area = msk.detectParticle(ptc, x0, y0, 0x80);
			}
		}
	}
	
	return particles;
}
		
int detect_particles_csv(unsigned char *mask, int hm, int wm, const char *csvfile)
{
	Raster8 msk(wm, hm, mask);
	std::vector<Slice> particles = detect_particles_2d(msk);
	write_particle_data(particles, csvfile);
	
	return int(particles.size());
}

// std::vector<std::vector<int>> res;
std::vector<std::vector<int>> detect_particles(unsigned char *mask, int hm, int wm)
{
	Raster8 msk(wm, hm, mask);
	std::vector<Slice> particles = detect_particles_2d(msk);
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

void postprocess_particle_borders(unsigned char *mask, int hm, int wm)
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
	for (int pass=0; pass<5; pass++) {
		for (Slice &ptc : particles)
			out.expandParticle(ptc, 0, 0x30, 0x80, HOOD_SIZE_NEUMANN);
	}

	out.replaceColor(out.forcedc, 0x80);
	out.expandBorders(0, 0xFF, HOOD_SIZE_MOORE, 0x80);

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
	
	if (postproc == POSTPROC_DNA) {
		// std::cout << "fix_borders_dna" << std::endl;
		aml.fix_borders_dna();
	} else if (postproc == POSTPROC_ACTIN) {
		// std::cout << "fix_borders_actin" << std::endl;
		aml.fix_borders_actin();
	}
	aml.fill_gaps();
	
	aml.paint_final();
	
	std::cout << "Write " << aml.cells.size() << " objects to " << csvfile << std::endl;
	write_cell_data(aml.cells, csvfile);

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
