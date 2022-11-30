#include "raster.h"
#include "laminb.h"

static double calc_rad_3d(Particle3D& cell, double xyscale, double zscale, double* pcx, double* pcy, double* pcz)
{
	double cx = 0., cy = 0., cz = 0., r = 0.;
	long long npix = 0;
	
	for (int z=0; size_t(z)<cell.fills.size(); z++) {
		double scz = z * zscale;
		for (HSeg& hs : cell.fills[z]) {
			double scy = hs.y * xyscale;
			for (int x=hs.xl; x<=hs.xr; x++) {
				cx += x*xyscale;
				cy += scy;
				cz += scz;
				++npix;
			}
		}
	}
	if (npix == 0) ++npix;
	cx /= npix;
	cy /= npix;
	cz /= npix;
	
	for (int z=0; size_t(z)<cell.fills.size(); z++) {
		double dz = z*zscale - cz;
		for (HSeg& hs : cell.fills[z]) {
			double dy = hs.y*xyscale - cy;
			for (int x=hs.xl; x<=hs.xr; x++) {
				double dx = x*xyscale - cx;
				r += (dx*dx + dy*dy + dz*dz);
			}
		}
	}
	r = sqrt(r/npix) * pow(2., 1./3.);
	*pcx = cx;
	*pcy = cy;
	*pcz = cz;
	return r;
}

static double calc_lmn_rad3d(Raster16_3D& dstack, Particle3D& cell, double thresh, double xyscale, double zscale,
	double cx, double cy, double cz)
{
	unsigned short ts = (unsigned short)thresh;
	double r = 0.;
	double tot = 0.;

	for (int z=0; size_t(z)<cell.fills.size(); z++) {
		double dz = z*zscale - cz;
		for (HSeg& hs : cell.fills[z]) {
			double dy = hs.y*xyscale - cy;
			unsigned short *b = dstack.scanLine(hs.y, z);
			for (int x=hs.xl; x<=hs.xr; x++) {
				if (b[x] < ts) continue;
				double wt = double(b[x]) / thresh;
				double dx = x*xyscale - cx;
				
				tot += wt;
				r += (dx*dx + dy*dy + dz*dz) * wt;
			}
		}
	}
	
	if (tot < 1.) tot = 1.;
	r = sqrt(r/tot);
	
	return r;
}

static Particle3D expand_particle3d(Raster8 &msk, Particle3D& cell)
{
	Particle3D res;
	res.fills.resize(cell.fills.size());
	
	Boundary bnd = cell.bnd.boundary2d();
	bnd.expand(2);
	msk.clip(bnd, 1);
	
	for (int z=0; size_t(z)<cell.fills.size(); z++) {
		if (cell.fills[z].empty()) continue;
		msk.fill(bnd, 0);
		msk.paintParticleFill(cell.fills[z], 0xC0);
		msk.expandBordersInto(bnd, 0xC0, 0, 0x50, HOOD_SIZE_MOORE);
		msk.rescanParticleFill(bnd, res.fills[z], 0xC0);
	}
	
	long long rvol = res.update_from_fill();
	//long long svol = cell.update_from_fill();
	
	//std::cout << "Src Vol: " << svol << " -> Res Vol: " << rvol << std::endl; 
	return res;
}

LaminB_Result lamin_b_analysis(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *csvfile, double thresh, double xyscale, double zscale)
{
	LaminB_Result res;
	
	std::vector<Particle3D> cells;
	if (read_cell_data(csvfile, cells, wd3d, hd3d, zd3d) < 0)
		return res;
	if (cells.empty())
		return res;
	std::cout << "Analyzing Lamin B1 of " << cells.size() << " cells from " << csvfile << std::endl;
	
	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);
	double pixvol = xyscale * xyscale * zscale;
	
	Raster8 msk(wd3d, hd3d, NULL);
	msk.fill(0);
	
	for (Particle3D& cell : cells) {
		int ch = 0;
		for (int z=0; size_t(z)<cell.fills.size(); z++)
			if (!cell.fills[z].empty()) ++ch;
		if (ch < 5) continue;
		double cx, cy, cz;
		double dnar = calc_rad_3d(cell, xyscale, zscale, &cx, &cy, &cz);
		res.dnarad.push_back(dnar);
		
		Particle3D xcell = expand_particle3d(msk, cell);
		
		double lmnr = calc_lmn_rad3d(dstack, xcell, thresh, xyscale, zscale, cx, cy, cz);
		res.lmnrad.push_back(lmnr);
		long long vol = 0;
		for (int z=0; size_t(z)<cell.fills.size(); z++) {
			vol += fill_area(cell.fills[z]);
		}
		res.volume.push_back(vol * pixvol);
	}

	
	return res;
}


template <class T>
class RasterF3D {
public:
	int w, h, d;
	T* buf;
	long long len;
	long long plane_len;
	bool is_temp;
	RasterF3D(int w, int h, int d, T* buf)
	{
		this->w = w;
		this->h = h;
		this->d = d;
		plane_len = (long long)(this->w) * this->h;
		len = plane_len * this->d;
		if (!buf) {
			this->buf = new T[len];
			is_temp = true;
		} else {
			this->buf = buf;
			is_temp = false;
		}
	}
	~RasterF3D() { if (is_temp) delete[] buf; }
	void fill(T c) {
		for (long long i=0; i<len; i++) buf[i] = c;
	}
	T *scanLine(int y, int z) { return buf + (z*plane_len + (long long)(y)*w); }
	T value(int x, int y, int z) { return *(scanLine(y, z) + x); }
	void setValue(int x, int y, int z, T v) { *(scanLine(y, z) + x) = v; }
	void fillBorderXY(T c, int bsz)
	{
		int x, y, z;
		for (z=0; z<d; z++) {
			for (y=0; y<bsz; y++) {
				T* ptop = scanLine(y, z);
				T* pbot = scanLine(h-1-y, z);
				for (x=0; x<w; x++) {
					ptop[x] = pbot[x] = c;
				}
			}
			for (y=bsz; y<h-bsz; y++) {
				T* p = scanLine(y, z);
				for (x=0; x<bsz; x++) {
					p[w-1-x] = p[x] = c;
				}
			}
		}
	}
	void paintParticle(Particle3D &pt, T c)
	{
		for (int z=0; size_t(z)<pt.fills.size(); z++) {
			for (HSeg &hs : pt.fills[z]) {
				T *p = scanLine(hs.y, z);
				for (int x=hs.xl; x<=hs.xr; x++)
					p[x] = c;
			}
		}
	}
};

struct NbrDist3D
{
	int dx, dy, dz;
	float dist;
	NbrDist3D() {}
	NbrDist3D(int _dx, int _dy, int _dz, float _dist) : dx(_dx), dy(_dy), dz(_dz), dist(_dist) {}
};

static bool comp_nbrdist3d(const NbrDist3D& a, const NbrDist3D& b)
{
	if (a.dz < b.dz) return true;
	if (a.dz > b.dz) return false;
	if (a.dist < b.dist) return true;
	if (a.dist > b.dist) return false;
	if (a.dy < b.dy) return true;
	if (a.dy > b.dy) return false;
	if (a.dx < b.dx) return true;
	return false;
};

static std::vector<NbrDist3D> calc_dist_3d(double xyscale, double zscale, double maxr)
{
	std::vector<NbrDist3D> nbd;
	for (int z=-3; z<=3; z++) {
		double dz = z * zscale;
		for (int y=-3; y<=3; y++) {
			double dy = y * xyscale;
			for (int x=-3; x<=3; x++) {
				double dx = x * xyscale;
				if (dx == 0 && dy == 0 && dz == 0) continue;
				double dist = sqrt(dz*dz + dy*dy + dx*dx);
				if (dist <= maxr)
					nbd.push_back(NbrDist3D(x, y, z, float(dist)));
			}
		}
	}
	std::sort(nbd.begin(), nbd.end(), comp_nbrdist3d);
	
	return nbd;
}

LaminB_Result_2 lamin_b_analysis_2(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *csvfile, double thresh, double xyscale, double zscale)
{
	LaminB_Result_2 res;

	std::vector<Particle3D> cells;
	if (read_cell_data(csvfile, cells, wd3d, hd3d, zd3d) < 0)
		return res;
	if (cells.empty())
		return res;
	std::cout << "Analyzing Lamin B1 of " << cells.size() << " cells from " << csvfile << std::endl;
	
	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);
	double pixvol = xyscale * xyscale * zscale;

	double maxr = xyscale * 2.95;
	double maxr_z = zscale * 2.95;
	if (maxr > maxr_z) maxr = maxr_z;
	
	std::vector<NbrDist3D> nbd = calc_dist_3d(xyscale, zscale, maxr);
	
	RasterF3D<float> fstack(wd3d, hd3d, zd3d, NULL);
	RasterF3D<float>* gstack = new RasterF3D<float>(wd3d, hd3d, zd3d, NULL);
	fstack.fill(0.);
	for (Particle3D& cell : cells)
		fstack.paintParticle(cell, 1000000.);
	fstack.fillBorderXY(0., 2);
	
	int iter;
	for (iter=0; iter<100; iter++) {
		//std::cout << "iter=" << iter << std::endl;
		memcpy(gstack->buf, fstack.buf, gstack->len * sizeof(float));
		long long nmod = 0;
		for (int z=0; z<fstack.d; z++) {
			int j0=-1, j1=-1;
			for (int j=0; size_t(j)<nbd.size(); j++) {
				int z1 = z + nbd[j].dz;
				if (j0<0 && z1>=0) j0 = j;
				if (z1 >= dstack.d) break;
				j1 = j;
			}
			
			for (int y=2; y<fstack.h-2; y++) {
				for (int x=2; x<fstack.w-2; x++) {
					float sdist = fstack.value(x, y, z);
					float cdist = sdist;
					if (sdist < 0.1) continue;
					for (int j=j0; j<=j1; j++) {
						NbrDist3D& nb = nbd[j];
						float _cdist = fstack.value(x+nb.dx, y+nb.dy, z+nb.dz) + nb.dist;
						if (_cdist < cdist) cdist = _cdist;
					}
					if (cdist < sdist) {
						gstack->setValue(x, y, z, cdist);
						++nmod;
					}
				}
			}
		}
		// std::cout << "iter=" << iter << " nmod=" << nmod << std::endl;
		if (nmod == 0) break;
		memcpy(fstack.buf, gstack->buf, gstack->len * sizeof(float));
	}
	
	delete gstack;
	
	Raster8 msk(wd3d, hd3d, NULL);
	msk.fill(0);
	
	unsigned short ts = (unsigned short)thresh;
	for (Particle3D& cell : cells) {
		int ch = 0;
		for (int z=0; size_t(z)<cell.fills.size(); z++)
			if (!cell.fills[z].empty()) ++ch;
		if (ch < 5) continue;
		
		double cx, cy, cz;
		double dnar = calc_rad_3d(cell, xyscale, zscale, &cx, &cy, &cz);
		res.dnarad.push_back(dnar);
		
		Particle3D xcell = expand_particle3d(msk, cell);

		long long vol = 0;
		double sdist = 0.;
		double sval = 0.;
		for (int z=0; size_t(z)<cell.fills.size(); z++) {
			if (cell.fills[z].empty()) continue;
			for (HSeg& hs : cell.fills[z]) {
				vol += (hs.xr - hs.xl +1);
			}
			for (HSeg& hs : xcell.fills[z]) {
				for (int x=hs.xl; x<=hs.xr; x++) {
					if (dstack.value(x, hs.y, z) < ts) continue;
					double val = double(dstack.value(x, hs.y, z));
					double dist = double(fstack.value(x, hs.y, z));
					sdist += dist*dist*val;
					sval += val;
				}
			}
		}
		if (sval > 0.1) sdist = sqrt(sdist/sval);
		res.lmndev.push_back(sdist);
		res.volume.push_back(vol * pixvol);
	}
	
/*
	for (long long i=0; i<fstack.len; i++) {
		float v = fstack.buf[i] * 50.;
		if (v > 255.) v = 255.;
		dstack.buf[i] = (unsigned short)v;
	}
*/

	std::cout << "Done in " << iter << " iterations." << std::endl;
	
	return res;
}

