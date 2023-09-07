#include "dsp.h"

//--- DSP_Detector: public

static void average_dat16(Raster16_3D& dstack, Raster16& davg, int z)
{
	Raster16 dat1 = dstack.getPlane(z-1);
	Raster16 dat2 = dstack.getPlane(z);
	Raster16 dat3 = dstack.getPlane(z+1);
	for (int y=0; y<dstack.h; y++) {
		unsigned short *p1 = dat1.scanLine(y);
		unsigned short *p2 = dat2.scanLine(y);
		unsigned short *p3 = dat3.scanLine(y);
		unsigned short *p = davg.scanLine(y);
		for (int x=0; x<dstack.w; x++) {
			p[x] = (unsigned short)((int(p1[x]) + int(p2[x]) + int(p3[x]))/3);
		}
	}
}

static bool comp_dsp_seeds(DSP_Seed& a, DSP_Seed& b) {
	return a.score > b.score;
}

static bool comp_dsp_scores(DSP_Score& a, DSP_Score& b) {
	return a.score > b.score;
}

void DSP_Detector::detect_upd()
{
	Raster16_3D dstack(w, h, d, NULL);
	Histogram hist = normalize_data(dstack);
	unsigned short otsu = hist.otsu16();
	unsigned short cut = hist.pix_level(uint64(hist.getCount()*0.99));
	
	std::cout << "Count: " << hist.getCount() << " OTSU: " << otsu << " x0.99: " << cut << std::endl;
	
	for (int z=0; z<d; z++) {
		Raster16 dat = dstack.getPlane(z);
		Raster8 msk = mstack.getPlane(z);
		for (int y=2; y<h-2; y++) {
			unsigned short *dp = dat.scanLine(y);
			unsigned char *dm = msk.scanLine(y);
			for (int x=2; x<w-2; x++) {
				if (dp[x] < otsu) continue;
				if (msk.countColors(x, y, HOOD_SIZE_FATCROSS, 0xFF) > 0) continue;
				dm[x] = (dp[x] >= cut && dat.is_local_max(x, y, HOOD_SIZE_FATCROSS)) ? 0xFF : 0x40;
			}
		}
	}
	
	double sc = 32/otsu;
	for (long long i=0; i<dstack.len; i++) {
		double v = sc * dstack.buf[i];
		if (v > 255.) v = 255.;
		srcstack.buf[i] = (unsigned short)v;
	}
}

std::vector<Particle3D> DSP_Detector::detect_particles()
{
	Raster16_3D dstack(w, h, d, NULL);
	Histogram hist = normalize_data(dstack);
	unsigned short otsu = hist.otsu16();
	unsigned short cut = hist.pix_level(uint64(hist.getCount()*0.99));

	for (int z=0; z<d; z++) {
		Raster16 dat = dstack.getPlane(z);
		Raster8 msk = mstack.getPlane(z);
		for (int y=2; y<h-2; y++) {
			unsigned short *dp = dat.scanLine(y);
			unsigned char *dm = msk.scanLine(y);
			for (int x=2; x<w-2; x++) {
				if (dp[x] < otsu) continue;
				if (msk.countColors(x, y, HOOD_SIZE_FATCROSS, 0xFF) > 0) continue;
				dm[x] = (dp[x] >= cut && dat.is_local_max(x, y, HOOD_SIZE_FATCROSS)) ? 0xFF : 0x40;
			}
		}
	}
	
	std::vector<DSP_Seed> dspseeds;
	
	for (int z=1; z<mstack.d-1; z++) {
		Raster8 msk = mstack.getPlane(z);
		for (int y=3; y<msk.h-3; y++) {
			unsigned char *p = msk.scanLine(y);
			for (int x=3; x<mstack.w-3; x++) {
				if (p[x] != 0xFF) continue;
				
				dspseeds.resize(dspseeds.size()+1);
				DSP_Seed& dsd = dspseeds[dspseeds.size()-1];
				dsd.pts.push_back(Point3D(x, y, z));
				p[x] = 0x40;
				int gap = 0;
				int z1 = z;
				int x0 = x;
				int y0 = y;
				while (gap < 3) {
					++z1;
					if (z1 >= mstack.d) break;
					Raster8 msk1 = mstack.getPlane(z1);
					for (int j=0; j<HOOD_SIZE_MOORE; j++) {
						int x1 = x0 + hood_pts[j].dx;
						int y1 = y0 + hood_pts[j].dy;
						if (msk1.value(x1, y1) == 0xFF) {
							dsd.pts.push_back(Point3D(x1, y1, z1));
							msk1.setValue(x1, y1, 0x40);
							x0 = x1;
							y0 = y1;
							gap = -1;
							break;
						}
					}
					++gap;
				}
				
				if (dsd.pts.size() < 3)
					dspseeds.resize(dspseeds.size()-1);
			}
		}
	}

	for (DSP_Seed& dsd : dspseeds) {
		for (Point3D& pt : dsd.pts) {
			Raster16 dat = dstack.getPlane(pt.z);
			for (int j=0; j<HOOD_SIZE_FATCROSS; j++) {
				dsd.score += double(dat.value(pt.x+hood_pts[j].dx, pt.y+hood_pts[j].dy));
			}
		}
	}
	
	std::sort(dspseeds.begin(), dspseeds.end(), comp_dsp_seeds);

	mstack.fill(0);
	for (DSP_Seed& dsd : dspseeds) {
		for (Point3D& pt : dsd.pts) {
			mstack.setValue(pt.x, pt.y, pt.z, 0xFF);
		}
	}

	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		Raster16 dat = dstack.getPlane(z);
		msk.expandBorders(0xFF, 0x40, HOOD_SIZE_RAD3, 0);
		for (int y=0; y<dat.h; y++) {
			hist.add_row16(dat.scanLine(y), dat.w, msk.scanLine(y));
		}
	}
	cut = (unsigned short)(hist.otsu16() / sensitivity);
	
	// std::cout << "Cut=" << cut << " OTSU=" << hist.otsu16() << std::endl;
	
	for (long long i=0; i<dstack.len; i++) {
		mstack.buf[i] = (dstack.buf[i] >= cut) ? 0x40 : 0;
	}
	
	std::vector<Particle3D> particles(dspseeds.size());
	for (size_t j=0; j<particles.size(); j++) {
		DSP_Seed& dsd = dspseeds[j];
		Particle3D& cell = particles[j];
		cell.fills.resize(mstack.d);
		for (Point3D& pt : dsd.pts) {
			cell.fills[pt.z].push_back(HSeg(pt.y, pt.x, pt.x));
		}
		cell.update_from_fill();
	}
	
	for (int pass=0; pass<3; pass++) {
		for (int z=0; z<mstack.d; z++) {
			Raster8 msk = mstack.getPlane(z);
			for (Particle3D& cell : particles) {
				std::vector<HSeg>& fill = cell.fills[z];
				if (fill.empty()) continue;
				msk.paintParticleFill(fill, 0xA0);
				Boundary bnd = cell.bnd.boundary2d();
				bnd.expand(1);
				msk.clip(bnd, 1);
				msk.expandBordersInto(bnd, 0xA0, 0x40, 0x51, HOOD_SIZE_NEUMANN, false);
				msk.rescanParticleFill(bnd, fill, 0xA0);
				msk.paintParticleFill(fill, 0xFF);
			}
		}
		for (Particle3D& cell : particles)
			cell.update_from_fill();
	}

	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		std::vector<int> zs;
		if (z > 0) zs.push_back(z-1);
		if (z < mstack.d-1) zs.push_back(z+1);
		for (Particle3D& cell : particles) {
			if (!cell.fills[z].empty()) continue;
			bool fixed = false;
			for (int z1 : zs) {
				if (!cell.fills[z1].empty()) {
					msk.paintParticleFillInto(cell.fills[z1], 0xA0, 0x40);
					fixed = true;
				}
			}
			if (!fixed) continue;
			Boundary bnd = cell.bnd.boundary2d();
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.rescanParticleFill(bnd, cell.fills[z], 0xA0);
			msk.paintParticleFill(cell.fills[z], 0xFF);
		}
	}
	
	int valid = 0;
	for (Particle3D& cell : particles) {
		long long vol = cell.update_from_fill();
		if (vol < min_volume) {
			//mstack.paintParticle(cell, 0x40);
			cell.clear();
			continue;
		}
		++valid;
	}
	
	mstack.fill(0);
	for (Particle3D& cell : particles) {
		if (cell.empty()) continue;
		mstack.paintParticle(cell, 0xFF);
	}
	
	for (size_t idx=0; idx<particles.size(); idx++) {
		Particle3D& cell = particles[idx];
		if (cell.empty()) continue;
		double sc = particle_snr(cell);
		if (sc < min_snr) {
			mstack.paintParticle(cell, 0);
			cell.clear();
			--valid;
		}
	}
	
	// std::cout << "Valid DSP objects: " << valid << " of " << particles.size() << std::endl;
	
	return particles;
}


/*
std::vector<Particle3D> DSP_Detector::detect_particles()
{
	Raster16_3D dstack(w, h, d, NULL);
	Raster16 davg(w, h, NULL);
	
	Histogram hist = normalize_data(dstack);
	unsigned short otsu = hist.otsu16();
	
	double sc_thresh = 1.;
	for (int pass=0; pass<3; pass++)
	{
		unsigned short cut = (unsigned short)(otsu * sc_thresh);
		sc_thresh /= sqrt(2.);
		
		for (int z=1; z<mstack.d-1; z++) {
			average_dat16(dstack, davg, z);
			Raster8 msk = mstack.getPlane(z);
			for (int y=3; y<mstack.h-3; y++) {
				unsigned short *dpt = davg.scanLine(y);
				unsigned char *mpt = msk.scanLine(y);
				for (int x=3; x<mstack.w-3; x++) {
					if (mpt[x] == 0 && dpt[x] >= cut) {
						if (msk.countColors(x, y, HOOD_SIZE_FATCROSS, 0xFF) > 0) {
							mpt[x] = 0x40;
							continue;
						}
						
						double xc, yc;
						davg.centerMass(x, y, HOOD_SIZE_RAD3, &xc, &yc);
						xc -= x;
						yc -= y;
						mpt[x] = (sqrt(xc*xc + yc*yc) < 0.75) ? 0xFF : 0x40;
					}
				}
			}
		}
	}
	
	std::vector<DSP_Seed> dspseeds;
	
	for (int z=1; z<mstack.d-1; z++) {
		Raster8 msk = mstack.getPlane(z);
		for (int y=3; y<msk.h-3; y++) {
			unsigned char *p = msk.scanLine(y);
			for (int x=3; x<mstack.w-3; x++) {
				if (p[x] != 0xFF) continue;
				
				dspseeds.resize(dspseeds.size()+1);
				DSP_Seed& dsd = dspseeds[dspseeds.size()-1];
				dsd.pts.push_back(Point3D(x, y, z));
				p[x] = 0x40;
				int gap = 0;
				int z1 = z;
				int x0 = x;
				int y0 = y;
				while (gap < 3) {
					++z1;
					if (z1 >= mstack.d) break;
					Raster8 msk1 = mstack.getPlane(z1);
					for (int j=0; j<HOOD_SIZE_MOORE; j++) {
						int x1 = x0 + hood_pts[j].dx;
						int y1 = y0 + hood_pts[j].dy;
						if (msk1.value(x1, y1) == 0xFF) {
							dsd.pts.push_back(Point3D(x1, y1, z1));
							msk1.setValue(x1, y1, 0x40);
							x0 = x1;
							y0 = y1;
							gap = -1;
							break;
						}
					}
					++gap;
				}
				
				if (dsd.pts.size() < 3)
					dspseeds.resize(dspseeds.size()-1);
			}
		}
	}

	for (DSP_Seed& dsd : dspseeds) {
		for (Point3D& pt : dsd.pts) {
			Raster16 dat = dstack.getPlane(pt.z);
			for (int j=0; j<HOOD_SIZE_FATCROSS; j++) {
				dsd.score += double(dat.value(pt.x+hood_pts[j].dx, pt.y+hood_pts[j].dy));
			}
		}
	}
	
	std::sort(dspseeds.begin(), dspseeds.end(), comp_dsp_seeds);

	mstack.fill(0);
	for (DSP_Seed& dsd : dspseeds) {
		for (Point3D& pt : dsd.pts) {
			mstack.setValue(pt.x, pt.y, pt.z, 0xFF);
		}
	}

	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		Raster16 dat = dstack.getPlane(z);
		msk.expandBorders(0xFF, 0x40, HOOD_SIZE_RAD3, 0);
		for (int y=0; y<dat.h; y++) {
			hist.add_row16(dat.scanLine(y), dat.w, msk.scanLine(y));
		}
	}
	unsigned short cut = (unsigned short)(hist.otsu16() / sensitivity);
	
	// std::cout << "Cut=" << cut << " OTSU=" << hist.otsu16() << std::endl;
	
	for (long long i=0; i<dstack.len; i++) {
		mstack.buf[i] = (dstack.buf[i] >= cut) ? 0x40 : 0;
	}
	
	std::vector<Particle3D> particles(dspseeds.size());
	for (size_t j=0; j<particles.size(); j++) {
		DSP_Seed& dsd = dspseeds[j];
		Particle3D& cell = particles[j];
		cell.fills.resize(mstack.d);
		for (Point3D& pt : dsd.pts) {
			cell.fills[pt.z].push_back(HSeg(pt.y, pt.x, pt.x));
		}
		cell.update_from_fill();
	}
	
	for (int pass=0; pass<3; pass++) {
		for (int z=0; z<mstack.d; z++) {
			Raster8 msk = mstack.getPlane(z);
			for (Particle3D& cell : particles) {
				std::vector<HSeg>& fill = cell.fills[z];
				if (fill.empty()) continue;
				msk.paintParticleFill(fill, 0xA0);
				Boundary bnd = cell.bnd.boundary2d();
				bnd.expand(1);
				msk.clip(bnd, 1);
				msk.expandBordersInto(bnd, 0xA0, 0x40, 0x51, HOOD_SIZE_NEUMANN, false);
				msk.rescanParticleFill(bnd, fill, 0xA0);
				msk.paintParticleFill(fill, 0xFF);
			}
		}
		for (Particle3D& cell : particles)
			cell.update_from_fill();
	}

	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		std::vector<int> zs;
		if (z > 0) zs.push_back(z-1);
		if (z < mstack.d-1) zs.push_back(z+1);
		for (Particle3D& cell : particles) {
			if (!cell.fills[z].empty()) continue;
			bool fixed = false;
			for (int z1 : zs) {
				if (!cell.fills[z1].empty()) {
					msk.paintParticleFillInto(cell.fills[z1], 0xA0, 0x40);
					fixed = true;
				}
			}
			if (!fixed) continue;
			Boundary bnd = cell.bnd.boundary2d();
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.rescanParticleFill(bnd, cell.fills[z], 0xA0);
			msk.paintParticleFill(cell.fills[z], 0xFF);
		}
	}
	
	int valid = 0;
	for (Particle3D& cell : particles) {
		long long vol = cell.update_from_fill();
		if (vol < min_volume) {
			//mstack.paintParticle(cell, 0x40);
			cell.clear();
			continue;
		}
		++valid;
	}
	
	mstack.fill(0);
	for (Particle3D& cell : particles) {
		if (cell.empty()) continue;
		mstack.paintParticle(cell, 0xFF);
	}
	
	for (size_t idx=0; idx<particles.size(); idx++) {
		Particle3D& cell = particles[idx];
		if (cell.empty()) continue;
		double sc = particle_snr(cell);
		if (sc < min_snr) {
			mstack.paintParticle(cell, 0);
			cell.clear();
			--valid;
		}
	}
	
	// std::cout << "Valid DSP objects: " << valid << " of " << particles.size() << std::endl;
	
	return particles;
}
*/


//--- DSP_Detector: protected

static double savgol_5_4[25] = {
	0.041632653061224545, -0.08081632653061216, 0.07836734693877527, -0.08081632653061206, 0.04163265306122459,
	-0.0808163265306123, -0.0195918367346938, 0.20081632653061213, -0.019591836734694154, -0.08081632653061227,
	0.07836734693877541, 0.20081632653061215, 0.4416326530612249, 0.20081632653061207, 0.07836734693877563,
	-0.08081632653061208, -0.01959183673469368, 0.20081632653061227, -0.01959183673469402, -0.08081632653061217,
	0.04163265306122452, -0.08081632653061209, 0.07836734693877559, -0.0808163265306121, 0.041632653061224385
};
static double savgol_5_2[25] = {
	-0.07428571428571419, 0.011428571428571392, 0.039999999999999994, 0.01142857142857142, -0.0742857142857143,
	0.011428571428571416, 0.09714285714285713, 0.12571428571428575, 0.09714285714285717, 0.011428571428571425,
	0.03999999999999998, 0.12571428571428572, 0.15428571428571436, 0.12571428571428578, 0.040000000000000036,
	0.011428571428571408, 0.09714285714285714, 0.12571428571428572, 0.09714285714285717, 0.011428571428571437,
	-0.07428571428571433, 0.011428571428571413, 0.04000000000000002, 0.011428571428571413, -0.0742857142857143
};

Histogram DSP_Detector::normalize_data(Raster16_3D& dstack)
{
	dstack.fill(0);
	for (int z=0; z<dstack.d; z++) {
		Raster16 sdat = srcstack.getPlane(z);
		Raster16 tdat = dstack.getPlane(z);
		for (int y0=3; y0<dstack.h-3; y0++) {
			unsigned short *tp = tdat.scanLine(y0);
			for (int x0=3; x0<dstack.w-3; x0++) {
				double rough = 0.;
				double smooth = 0.;
				
				for (int j=0; j<HOOD_SIZE_RAD3; j++)
					rough += double(sdat.value(x0+hood_pts[j].dx, y0+hood_pts[j].dy));
				rough /= HOOD_SIZE_RAD3;
				
				double *sgpt = savgol_5_4;
				for (int y=y0-2; y<=y0+2; y++) {
					for (int x=x0-2; x<=x0+2; x++) {
						smooth += sdat.value(x,y) * (*sgpt++);
					}
				}
				
				smooth = (smooth - rough)*0.5;
				if (smooth < 0.) smooth = 0.;
				else if (smooth > 65530.) smooth = 65530.;
				tp[x0] = (unsigned short)smooth;
			}
		}
	}

	Histogram hist(0x4000);
	for (int z=0; z<dstack.d; z++) {
		for (int y=0; y<dstack.h; y++) {
			hist.add_row16(dstack.scanLine(y, z), dstack.w);
		}
	}
	return hist;
}

double DSP_Detector::particle_snr(Particle3D& cell)
{
	Boundary3D bnd3d = cell.bnd;
	bnd3d.expand(3);
	mstack.clip(bnd3d, 3);
	Boundary bnd = bnd3d.boundary2d();
	
	std::vector<unsigned short> bkg_pixels;
	double fg_mean = 0.;
	long long fg_cnt = 0;

	mstack.paintParticle(cell, 0xA0);
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		Raster16 dat = srcstack.getPlane(z);
		
		std::vector<HSeg>& fill = cell.fills[z];
		if (fill.empty()) continue;
		for (HSeg& hs : fill) {
			for (int x=hs.xl; x<=hs.xr; x++) {
				fg_mean += double(dat.value(x, hs.y));
				++fg_cnt;
			}
		}
		
		for (int y=bnd.ymin; y<=bnd.ymax; y++) {
			unsigned char *p = msk.scanLine(y);
			for (int x=bnd.xmin; x<=bnd.xmax; x++) {
				if (p[x] != 0) continue;
				for (int j=1; j<HOOD_SIZE_RAD3; j++) {
					if (msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy) == 0xA0) {
						bkg_pixels.push_back(dat.value(x, y));
						break;
					}
				}
			}
		}
	}
	mstack.paintParticle(cell, 0xFF);
	
	if (fg_cnt == 0 || bkg_pixels.size() == 0) return 0.;
	fg_mean /= fg_cnt;

	double bg_mean = 0.;
	for (unsigned short v : bkg_pixels)
		bg_mean += double(v);
	bg_mean /= bkg_pixels.size();

	double bg_stddev = 0.;
	for (unsigned short v : bkg_pixels) {
		double dv = double(v) - bg_mean;
		bg_stddev += dv*dv;
	}
	bg_stddev = sqrt(bg_stddev/bkg_pixels.size() + 1.);

	return (fg_mean - bg_mean) / bg_stddev;
}


// DSP_Analyzer: protected

long long DSP_Analyzer::pix_distances(std::vector<long long>& pix_dist_inside,
		std::vector<long long>& pix_dist_outside, int maxdist)
{
	long long tot = 0;
	
	pix_dist_inside.resize(maxdist);
	pix_dist_outside.resize(maxdist);
	for (size_t i=0; i<pix_dist_inside.size(); i++) {
		pix_dist_inside[i] = pix_dist_outside[i] = 0;
	}
	
	Raster8 msk0 = mstack.getPlane(0);
	for (Particle3D& dsp : dsps) {
		tot += dsp.volume();
		Boundary bnd = dsp.bnd.boundary2d();
		bnd.expand(maxdist);
		msk0.clip(bnd);
		std::vector<int> idxs;
		for (int idx=0; size_t(idx)<cells.size(); idx++) {
			Boundary b = cells[idx].bnd.boundary2d();
			if (bnd.intersects(b))
				idxs.push_back(idx);
		}
		if (idxs.empty()) continue;
		
		for (int z=0; z<d; z++) {
			std::vector<HSeg>& dspfill = dsp.fills[z];
			if (dspfill.empty()) continue;
			Raster8 msk = mstack.getPlane(z);
			for (HSeg& hs : dspfill) {
				for (int x=hs.xl; x<=hs.xr; x++) {
					bool is_inside = msk.value(x, hs.y) != 0;
					double bsf = 1000000.;
					for (int idx : idxs) {
						CellBorder& brd = borders[idx];
						for (Point& pt : brd.bpixels[z]) {
							double dist = pt.dist(x, hs.y);
							if (dist < bsf) bsf = dist;
						}
					}
					int dist = int(bsf);
					if (dist >= maxdist) continue;
					if (is_inside)
						++pix_dist_inside[dist];
					else
						++pix_dist_outside[dist];
				}
			}
		}
	}
	
	return tot;
}

std::vector<DSPtoCell> DSP_Analyzer::dsp_to_cell_match(int maxdist, double close_by)
{
	std::vector<DSPtoCell> res(dsps.size());
	for (int idx=0; size_t(idx)<dsps.size(); idx++) {
		DSPtoCell& dc = res[idx];
		dc = find_closest_cell(idx, maxdist);
		if (dc.border_dist >= close_by && dc.border_dist < 0.)
			dc.border_dist = 0.;
	}
	return res;
}

void DSP_Analyzer::paint_dsp_mask(std::vector<DSPtoCell>& dsptocells)
{
	mstack.fill(0);
	for (Particle3D& cell : cells)
		mstack.paintParticle(cell, 0x40);
	for (CellBorder& brd : borders) {
		for (int z=0; z<d; z++) {
			for (Point& pt : brd.bpixels[z]) {
				mstack.setValue(pt.x, pt.y, z, 0x60);
			}
		}
	}

	int rad = 3;
	double rsq = double(rad) * rad;
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		Boundary bnd0 = msk.getBoundary();
		for (DSPtoCell& dc : dsptocells) {
			if (z < dc.zmin || z > dc.zmax) continue;
			unsigned char c = (dc.border_dist < -1.5) ? 0x80 : 0xFF;
			int xc = int(dc.xc);
			int yc = int(dc.yc);
			Boundary bnd(xc-rad, yc-rad, xc+rad, yc+rad);
			msk.clip(bnd);
			for (int y=bnd.ymin; y<=bnd.ymax; y++) {
				double dysq = double(y - yc); dysq *= dysq;
				for (int x=bnd.xmin; x<=bnd.xmax; x++) {
					double dxsq = double(x - xc); dxsq *= dxsq;
					if (fabs(dxsq+dysq - rsq) < 2.)
						msk.setValue(x, y, c);
				}
			}
			if (z == int(dc.zc)) {
				msk.setValue(xc, yc, c);
			}
			
			if (dc.cellidx < 0) continue;
			Point3D& cc = centroids[dc.cellidx];
			int dx = cc.x - xc;
			int dy = cc.y - yc;
			if (dx==0 && dy==0) continue;
			double slx = 1., sly = 1.;
			double ax = 0., ay = 0.;
			if (dx < 0) {
				slx = -1.;
				dx = -dx;
				ax = 0.5;
			}
			if (dy < 0) {
				sly = -1.;
				dy = -dy;
				ay = 0.5;
			}
			if (dy > dx) {
				slx *= double(dx) / double(dy);
			} else {
				sly *= double(dy) / double(dx);
			}
			for (int r=3; r<6; r++) {
				Point p(int(r*slx+ax+xc), int(r*sly+ay+yc));
				if (!bnd0.IsInside(p)) break;
				msk.setValue(p.x, p.y, c);
			}
		}
	}
}

// DSP_Analyzer: protected

void DSP_Analyzer::detect_cell_borders()
{
	mstack.fill(0);
	borders.clear();
	borders.resize(cells.size());
	centroids.resize(cells.size());
	for (int idx=0; size_t(idx)<cells.size(); idx++) {
		Particle3D& cell = cells[idx];
		centroids[idx] = cell.center_mass();
		CellBorder& brd = borders[idx];
		brd.idx = idx;
		brd.bpixels.resize(d);
		for (int z=0; z<d; z++) {
			std::vector<HSeg>& fill = cell.fills[z];
			if (fill.empty()) continue;
			std::vector<Point>& bpix = brd.bpixels[z];
			Raster8 msk = mstack.getPlane(z);
			msk.paintParticleFill(fill, 0x80);
			for (HSeg& hs : fill) {
				for (int x=hs.xl; x<=hs.xr; x++) {
					for (int j=1; j<HOOD_SIZE_NEUMANN; j++) {
						if (msk.value(x+hood_pts[j].dx, hs.y+hood_pts[j].dy) != 0x80) {
							bpix.push_back(Point(x, hs.y));
							break;
						}
					}
				}
			}
			msk.paintParticleFill(fill, 0x40);
		}
	}
}

static bool comp_dsp_to_cell_cands(DSPtoCell& a, DSPtoCell& b) {
	if (a.overlay > b.overlay) return true;
	if (a.overlay < b.overlay) return false;
	return a.border_dist > b.border_dist;
}

static double closest_point_2d(std::vector<Point>& bpts, double xc, double yc) {
	double bsf = 1000000000.;
	for (Point& p : bpts) {
		double dx = xc - p.x;
		double dy = yc - p.y;
		double dsq = dx*dx + dy*dy;
		if (dsq < bsf) bsf = dsq;
	}
	return sqrt(bsf);
}

DSPtoCell DSP_Analyzer::find_closest_cell(int dspidx, int maxdist)
{
	DSPtoCell res;
	res.dspidx = dspidx;
	Particle3D& dsp = dsps[dspidx];
	res.zmin = dsp.bnd.zmin;
	res.zmax = dsp.bnd.zmax;
	res.value = srcstack.centerMass(dsp, &res.xc, &res.yc, &res.zc);
	res.volume = dsp.volume();
	res.cellidx = -1;
	res.border_dist = -1000000.;
	
	Point3D pc(int(res.xc+0.5), int(res.yc+0.5), int(res.zc+0.5));

	Boundary3D bnd = dsp.bnd;
	bnd.expand(maxdist);
	bnd.zmin = dsp.bnd.zmin-2;
	bnd.zmax = dsp.bnd.zmax+2;
	mstack.clip(bnd);
	
	std::vector<DSPtoCell> candidates;
	for (int cellidx=0; size_t(cellidx)<cells.size(); cellidx++) {
		Particle3D& cell = cells[cellidx];
		if (!bnd.intersects(cell.bnd)) continue;
		if (cell.contains(pc)) {
			res.cellidx = cellidx;
			res.overlay = cell.overlay_volume(dsp);
			res.border_dist = borders[cellidx].closest(res.xc, res.yc, res.zc);
			return res;
		}
		DSPtoCell cand = res;
		cand.cellidx = cellidx;
		cand.overlay = cell.overlay_volume(dsp);
		cand.border_dist = -borders[cellidx].closest(res.xc, res.yc, res.zc);
		candidates.push_back(cand);
	}
	
	if (candidates.size() > 0) {
		std::sort(candidates.begin(), candidates.end(), comp_dsp_to_cell_cands);
		DSPtoCell& c1 = candidates[0];
		if (c1.border_dist < -1.) {
			for (int z=0; z<d; z++) {
				std::vector<HSeg> fill = dsp.fills[z];
				if (fill.empty()) continue;
				Raster16 dat = srcstack.getPlane(z);
				double xc, yc;
				dat.centerMass(fill, &xc, &yc);
				Point3D ptc(int(xc+0.5), int(yc+0.5), z);
				for (DSPtoCell& dc : candidates) {
					Particle3D& cell = cells[dc.cellidx];
					if (cell.contains(ptc)) {
						res = dc;
						res.border_dist = 0.;
						return res;
					}
					double dist = -closest_point_2d(borders[dc.cellidx].bpixels[z], xc, yc);
					if (dist > dc.border_dist) dc.border_dist = dist;
				}
			}
			std::sort(candidates.begin(), candidates.end(), comp_dsp_to_cell_cands);
		}
	}
	
	if (!candidates.empty()) {
		res = candidates[0];
	}
	return res;
}






