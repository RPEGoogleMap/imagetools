#ifndef raster_h
#define raster_h

#include "geom.h"

const double SEGM_EPS = 0.00001;

class Histogram
{
protected:
	int nbins;
	uint64* bins;
	double* weights;
	double* jweights;
	uint64 ccount;
	int hilvl;
	int state;

	void init();
	void finish();
public:
	Histogram(int _nbins);
	virtual ~Histogram();
	int getNbins() { return nbins; }
	uint64 getCount() { return ccount; }
	unsigned short high_level() { return (unsigned short)(uint64(hilvl) * 0x10000L / uint64(nbins)); }
	unsigned short otsu16() { return (unsigned short)(uint64(otsu_bin()) * 0x10000L / uint64(nbins)); }
	unsigned char otsu8() { return (unsigned char)(uint64(otsu_bin()) * 0x100L / uint64(nbins)); }

	void add_row16(unsigned short *buf, int len, unsigned char *mask=NULL);
	void add_row8(unsigned char *buf, int len);
	int otsu_bin();
	// 16-bit level at which pixel count <pcount> is reached
	unsigned short pix_level(uint64 pcount);
};

template <class T>
class SummedArea {
protected:
	int w, h;
	T* buf;
	size_t len;
	//
	double sum_area(Boundary& bnd) {
		if (bnd.xmin == 0 || bnd.ymin == 0) {
			double s = double(value(bnd.xmax, bnd.ymax));
			if (bnd.xmin > 0) s -= double(value(bnd.xmin - 1, bnd.ymax));
			if (bnd.ymin > 0) s -= double(value(bnd.xmax, bnd.ymin - 1));
			return s;
		}
		return double(value(bnd.xmax, bnd.ymax)) - double(value(bnd.xmax, bnd.ymin - 1))
			+ double(value(bnd.xmin - 1, bnd.ymin - 1)) - double(value(bnd.xmin - 1, bnd.ymax));
	}
	double avg_area(Boundary& bnd) { return sum_area(bnd) / double(bnd.area()); }
	//
public:
	SummedArea(int w, int h) {
		this->w = w;
		this->h = h;
		len = size_t(this->w) * size_t(this->h);
		buf = new T[len];
	}
	~SummedArea() {
		delete[] buf;
	}
	T* scanLine(int y) { return buf + (size_t(w) * y); }
	T value(int x, int y) { return *(scanLine(y) + x); }
	// Add a row of pixels. Must be called sequentially, y=0,1,2, etc.
	template<class P>
	void add_line(int y, P* s) {
		T* t = scanLine(y);
		if (y == 0) {
			double csum = 0;
			for (int x = 0; x < w; x++) {
				csum += double(s[x]);
				t[x] = T(csum);
			}
		}
		else {
			T* prev = scanLine(y - 1);
			t[0] = T(double(prev[0]) + double(s[0]));
			for (int x = 1; x < w; x++) {
				t[x] = T(double(t[x - 1]) - double(prev[x - 1]) + double(prev[x]) + double(s[x]));
			}
		}
	}
	double sum_around(int x, int y, int dist) {
		Boundary bnd = boundary_around(x, y, dist, w, h);
		return sum_area(bnd);
	}
	double avg_around(int x, int y, int dist) {
		Boundary bnd = boundary_around(x, y, dist, w, h);
		return avg_area(bnd);
	}
};

struct Raster16 {
	int w, h;
	unsigned short *buf;
	long long len;
	bool is_temp;
	Raster16(int w, int h, unsigned short *buf)
	{
		this->w = w;
		this->h = h;
		this->len = (long long)(this->w) * this->h;
		if (!buf) {
			this->buf = new unsigned short[len];
			is_temp = true;
		} else {
			this->buf = buf;
			is_temp = false;
		}
	}
	~Raster16() { if (is_temp) delete[] buf; }
	Boundary getBoundary() { return Boundary(0, 0, w-1, h-1); }
	void fill(unsigned short c) {
		for (long long i=0; i<len; i++) buf[i] = c;
	}
	unsigned short *scanLine(int y) { return buf + ((long long)(w) * y); }
	unsigned short value(int x, int y) { return *(scanLine(y) + x); }
	void setValue(int x, int y, unsigned short v) { *(scanLine(y) + x) = v; }
	void setValue(const Point& p, unsigned short v) { *(scanLine(p.y) + p.x) = v; }
	void clip(Boundary &bnd, int bord=0) {
		if (bnd.xmin < bord) bnd.xmin = bord;
		if (bnd.xmax >= w-bord) bnd.xmax = w-1-bord;
		if (bnd.ymin < bord) bnd.ymin = bord;
		if (bnd.ymax >= h-bord) bnd.ymax = h-1-bord;
	}
	void paintParticle(Particle &ptc, unsigned short c) {
		paintParticleFill(ptc.fill, c);
	}
	void countColors(std::vector<ColorCounter> &res, Particle &ptc) {
		countColors(res, ptc.fill);
	}
	
	double mean_std(double *p_std=NULL);
	unsigned short avg_value(int x, int y, int nbsz=HOOD_SIZE_MOORE);
	void centerMass(int x0, int y0, int nbsz, double *p_xc, double *p_yc);
	double centerMass(std::vector<HSeg> &fill, double *p_xc, double *p_yc);
	bool is_local_max(int x, int y, int nbsz=HOOD_SIZE_FATCROSS);
	int local_rank(int x, int y, int nbsz=HOOD_SIZE_RAD3);
	void paintParticleFill(std::vector<HSeg> &fill, unsigned short c);
	void countColors(std::vector<ColorCounter> &res, std::vector<HSeg> &fill);
	void fill(Boundary &bnd, unsigned short c);
	void paintParticleADD(Particle &ptc, int incr);
	int rescanParticleThresh(Particle &ptc, unsigned short tc);
	void paintParticleInto(Particle &ptc, unsigned short c, unsigned short bkg);

};

class Raster8 {
public:
	int w, h;
	unsigned char *buf;
	long long len;
	bool is_temp;
	unsigned char forcedc;
	Raster8(int w, int h, unsigned char *buf)
	{
		this->w = w;
		this->h = h;
		this->len = (long long)(this->w) * this->h;
		if (!buf) {
			this->buf = new unsigned char[len];
			is_temp = true;
		} else {
			this->buf = buf;
			is_temp = false;
		}
		forcedc = 0xFE;
	}
	~Raster8() { if (is_temp) delete[] buf; }
	Boundary getBoundary() { return Boundary(0, 0, w-1, h-1); }
	unsigned char *scanLine(int y) { return buf + ((long long)(w) * y); }
	void fill(unsigned char c) { memset(buf, c, len); }
	unsigned char value(int x, int y) { return *(scanLine(y) + x); }
	void setValue(int x, int y, unsigned char c) { *(scanLine(y)+x) = c; }
	void setValue(const Point& p, unsigned char c) { *(scanLine(p.y) + p.x) = c; }
	void setMaxValue(int x, int y, unsigned char c) {
		unsigned char *p = scanLine(y) + x;
		if (*p < c) *p = c;
	}
	void clip(Boundary &bnd, int bord=0) {
		if (bnd.xmin < bord) bnd.xmin = bord;
		if (bnd.xmax >= w-bord) bnd.xmax = w-1-bord;
		if (bnd.ymin < bord) bnd.ymin = bord;
		if (bnd.ymax >= h-bord) bnd.ymax = h-1-bord;
	}
	long long fillParticle(int x0, int y0, unsigned char newc) {
		std::vector<HSeg> fill;
		findParticleFill(fill, x0, y0, newc);
		return fill_area(fill);
	}
	void paintParticle(Particle &ptc, unsigned char c) { paintParticleFill(ptc.fill, c); }
	void paintParticleInto(Particle &ptc, unsigned char c, unsigned char bkg) {
		paintParticleFillInto(ptc.fill, c, bkg);
	}
	void paintParticleADD(Particle &ptc, int incr) {
		paintParticleFillADD(ptc.fill, incr);
	}

	void fill(Boundary &bnd, unsigned char c);
	void fillBorder(unsigned char c, int bsz=1);
	void paintBoundary(const Boundary& bnd, unsigned char c);
	void replaceColor(unsigned char oldc, unsigned char newc);
	void replaceColor(Boundary &b, unsigned char oldc, unsigned char newc);
	void paintSegmentOut(Point pt1, Point pt2, unsigned char c, unsigned char cb);
	void paintPath(std::vector<Point> &path, unsigned char c, unsigned char cb);
	void paintParticleFill(std::vector<HSeg> &fill, unsigned char c);
	void paintParticleFillInto(std::vector<HSeg> &fill, unsigned char c, unsigned char bkg);
	void paintParticleFillADD(std::vector<HSeg> &fill, int incr);
	void paintParticlePat(Particle &ptc, unsigned char fg);
	void paintSegment(Point pt1, Point pt2, unsigned char c);
	void fillContour(Contour &cont, unsigned char fill_c, unsigned char bord_c);
	void paintContourBorder(Contour &cont, unsigned char bord_c);
	// replace c0 with newc if next to oldc
	int expandBorders(unsigned char oldc, unsigned char newc,
			int nbsz=HOOD_SIZE_RAD3, unsigned char c0=0);
	void expandBordersInto(Boundary &bnd, unsigned char fg, unsigned char bk, unsigned char tmpc,
			int nbsz=HOOD_SIZE_NEUMANN, bool keeptmpc=false);
	void filterParticles(unsigned char oldc, unsigned char newc, int minsz, unsigned char bk);
	void filterParticles(Boundary& bnd, unsigned char oldc, unsigned char newc, int minsz, unsigned char bk);
	bool touches(int x, int y, unsigned char c, int nbsz=HOOD_SIZE_NEUMANN);
	int countColors(int x0, int y0, int nbsz, unsigned char c);
	int countColorsNb(std::vector<Point> path, unsigned char c, int nbsz=1);
	int detectParticle(Particle &ptc, int x0, int y0, unsigned char newc);
	void shrinkParticle(Particle &ptc, unsigned char fillc, unsigned char tmpc,
			unsigned char bordc, int nbsz=HOOD_SIZE_MOORE);
	void expandParticle(Particle &ptc, unsigned char fillc, unsigned char tmpc,
			unsigned char bordc, int nbsz=HOOD_SIZE_NEUMANN);
	int countColor(std::vector<HSeg>& fill, unsigned char c);
	int countColor(Particle &ptc, unsigned char c) { return countColor(ptc.fill, c); }
	void rescanParticleFill(Boundary &bnd, std::vector<HSeg> &fill, unsigned char c);
	int rescanParticle(Particle &ptc, unsigned char c);
	int removeFalseWalls(Particle &ptc, unsigned char fg, unsigned char bk, unsigned char tmpc);
	int mergeParticles(Particle &ptc, Particle &other,
			unsigned char fgd, unsigned char bkg, unsigned char tmpc);

	void findParticleFill(std::vector<HSeg> &fill, int x0, int y0, unsigned char newc);
	void bordersAround(Boundary &bnd, unsigned char fg, unsigned char bk, unsigned char bordc,
			int hoodsz=HOOD_SIZE_NEUMANN);

protected:
	inline void paintHSeg(HSeg &hs, unsigned char c) {
		unsigned char *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) p[x] = c;
	}
	void findFillSegments(std::vector<HSeg> &fill, int y, int xl, int xr,
			unsigned char oldc, unsigned char newc);
};

class Raster3D
{
public:
	int w, h, d;
	unsigned char *buf;
	long long len;
	long long plane_len;
	bool is_temp;
	//
	Raster3D(int _w, int _h, int _d, unsigned char *_buf) :
			w(_w), h(_h), d(_d), buf(_buf) {
		plane_len = (long long)(w) * h;
		len = plane_len * d;
		if (!buf) {
			buf = new unsigned char[len];
			is_temp = true;
		} else {
			is_temp = false;
		}
	}
	virtual ~Raster3D() { if (is_temp) delete [] buf; }
	//
	Raster8 getPlane(int z) {
		return Raster8(w, h, buf + (plane_len *z));
	}
	void fill(unsigned char c) { memset(buf, c, len); }
	unsigned char value(int x, int y, int z) { return *(buf +(plane_len*z + (long long)(w)*y + x)); }
	void setValue(int x, int y, int z, unsigned char c) {
		*(buf +(plane_len*z + (long long)(w)*y + x)) = c;
	}
	unsigned char *scanLine(int y, int z) {
		return buf + (plane_len*z + (long long)(w)*y);
	}
	Boundary3D getBoundary() { return Boundary3D(0, 0, 0, w-1, h-1, d-1); }
	void clip(Boundary3D &bnd, int bord=0) {
		if (bnd.xmin < bord) bnd.xmin = bord;
		if (bnd.xmax >= w-bord) bnd.xmax = w-1-bord;
		if (bnd.ymin < bord) bnd.ymin = bord;
		if (bnd.ymax >= h-bord) bnd.ymax = h-1-bord;
		if (bnd.zmin < bord) bnd.zmin = bord;
		if (bnd.zmax >= d-bord) bnd.zmax = d-1-bord;
	}
	void replaceColor(unsigned char oldc, unsigned char newc) {
		for (long long i=0; i<len; i++)
			if (buf[i] == oldc) buf[i] = newc;
	}
	void replaceColorInv(unsigned char oldc, unsigned char newc) {
		for (long long i=0; i<len; i++)
			if (buf[i] != oldc) buf[i] = newc;
	}
	//
	void replaceColor(Boundary3D& bnd, unsigned char oldc, unsigned char newc);
	void expandBorders(Boundary3D& bnd, unsigned char fg, unsigned char bk, unsigned char bordc, int nbsz=HOOD3D_26);
	void expandBorders(unsigned char fg, unsigned char bk, unsigned char bordc, int nbsz=HOOD3D_26)
	{
		Boundary3D bnd = getBoundary();
		expandBorders(bnd, fg, bk, bordc, nbsz);
	}
	void chopBorders(unsigned char fg, unsigned char bk, unsigned char bordc, int nbsz=HOOD3D_6);
	void paintParticle(Particle3D &pt, unsigned char c);
	void paintParticleADD(Particle3D &pt, int incr);
	void paintParticleInto(Particle3D &pt, unsigned char newc, unsigned char oldc);
	void forceSeparation(Particle3D &pt, unsigned char fg, unsigned char bk, unsigned char tmpc);
	void forceBorder(Particle3D& pt, unsigned char fg, unsigned char bordc,
			unsigned char tmp_fg, unsigned char tmp_bordc);
	void rescanShrunkParticles(std::vector<Particle3D> &cells, unsigned char fg, unsigned char tmpc);
	void rescanShrunkCells(std::vector<Cell> &cells, unsigned char fg, unsigned char tmpc);
	void sandPaperCells(std::vector<Particle3D> &cells);
	long long rescanParticle(Particle3D& cell, unsigned char c);
	long long detectParticle(Particle3D& cell, int x0, int y0, int z0, unsigned char newc);
protected:
	long long findParticleFillsZ(std::vector<std::pair<int, std::vector<HSeg>>>& zfills,
			int x0, int y0, int z0, unsigned char newc);
	void findFillZSlices(std::vector<std::pair<int, std::vector<HSeg>>>& zfills,
			int z0, std::vector<HSeg>& fill, unsigned char oldc, unsigned char newc);
	Particle3D findBiggestCell(Boundary3D& bnd, unsigned char fg, unsigned char tmpc);
};

class Raster16_3D
{
public:
	int w, h, d;
	unsigned short *buf;
	long long len;
	long long plane_len;
	bool is_temp;
	//
	Raster16_3D(int _w, int _h, int _d, unsigned short *_buf) :
			w(_w), h(_h), d(_d), buf(_buf) {
		plane_len = (long long)(w) * h;
		len = plane_len * d;
		if (!buf) {
			buf = new unsigned short[len];
			is_temp = true;
		} else {
			is_temp = false;
		}
	}
	virtual ~Raster16_3D() { if (is_temp) delete [] buf; }
	//
	void fill(unsigned short c) {
		for (long long i=0; i<len; i++) buf[i] = c;
	}
	Raster16 getPlane(int z) {
		return Raster16(w, h, buf + (plane_len *z));
	}
	unsigned short value(int x, int y, int z) { return *(buf +(plane_len*z + (long long)(w)*y + x)); }
	void setValue(int x, int y, int z, unsigned short c) {
		*(buf +(plane_len*z + (long long)(w)*y + x)) = c;
	}
	unsigned short *scanLine(int y, int z) {
		return buf +(plane_len*z + (long long)(w)*y);
	}
	unsigned short otsu();
	double centerMass(Particle3D& cell, double *p_xc, double *p_yc, double *p_zc);
};

std::vector<NbrPoint3D> hood_sorted_by_z(int nbsz=HOOD3D_26);

#endif
