#include "raster.h"

//----------------------------- Histogram --------------------------------------

Histogram::Histogram(int _nbins) : nbins(_nbins)
{
	bins = new uint64[nbins];
	weights = new double[nbins];
	jweights = new double[nbins];
	init();
}

Histogram::~Histogram()
{
	delete [] bins;
	delete [] weights;
	delete [] jweights;
}

void Histogram::init() {
	memset(bins, 0, nbins * sizeof(uint64));
	memset(weights, 0, nbins * sizeof(double));
	memset(jweights, 0, nbins * sizeof(double));
	ccount = 0;
	hilvl = 0;
	state = 0;
}

void Histogram::finish()
{
	state = 2;
	if (ccount < 1) return;
	for (int j=0; j<nbins; j++) {
		weights[j] = double(bins[j]) / ccount;
		jweights[j] = j * weights[j];
	}
}

void Histogram::add_row16(unsigned short *buf, int len)
{
	if (state != 0 && state != 1) init();
	state = 1;
	for (int i=0; i<len; i++) {
		int v = int(uint64(buf[i]) * uint64(nbins) / 0x10000L);
		if (hilvl < v) hilvl = v;
		++bins[v];
	}
	if (len > 0)
		ccount += len;
}

void Histogram::add_row8(unsigned char *buf, int len)
{
	if (state != 0 && state != 1) init();
	state = 1;
	for (int i=0; i<len; i++) {
		int v = int(uint64(buf[i]) * uint64(nbins) / 0x100L);
		if (hilvl < v) hilvl = v;
		++bins[v];
	}
	if (len > 0)
		ccount += len;
}

int Histogram::otsu_bin()
{
	if (state != 2) finish();
	if (ccount == 0) return 0;
	int thresh = 0;
	double diff_max = 0.;
	double wB, wF, uB, uF;
	for (int i = 0; i < nbins; i++) {
		wB = wF = uB = uF = 0.;
		for (int j = 0; j < nbins; j++) {
			if (j <= i) {
				wB += weights[j];
				uB += jweights[j];
			}
			else {
				wF += weights[j];
				uF += jweights[j];
			}
		}
		if (wB == 0. || wF == 0.) continue;
		uB /= wB;
		uF /= wF;
		double diff = uF - uB;
		diff = wB * wF * diff * diff;
		if (diff > diff_max) {
			diff_max = diff;
			thresh = i;
		}
	}
	return thresh;
}

unsigned short Histogram::pix_level(uint64 pcount)
{
	uint64 cnt = 0;
	int lvl;
	for (lvl = 0; lvl < nbins; lvl++) {
		cnt += bins[lvl];
		if (cnt >= pcount) break;
	}
	return (unsigned short)(uint64(lvl) * 0x10000L / uint64(nbins));
}

//----------------------------- Raster16 --------------------------------------

double Raster16::mean_std(double *p_std)
{
	double sumv = 0.;
	
	for (long long i=0; i<len; i++) sumv += double(buf[i]);
	sumv /= len;
	
	if (p_std) {
		double sumsq = 0.;
		for (long long i=0; i<len; i++) {
			double dv = double(buf[i]) - sumv;
			sumsq += dv*dv;
		}
		*p_std = sqrt(sumsq / len);
	}
	
	return sumv;
}

unsigned short Raster16::avg_value(int x, int y, int nbsz)
{
	uint64 sumv = 0;
	for (int j=0; j<nbsz; j++)
		sumv += uint64(value(x+hood_pts[j].dx, y+hood_pts[j].dy));
	return (unsigned short)(sumv / nbsz);
}

void Raster16::centerMass(int x0, int y0, int nb, double *p_xc, double *p_yc)
{
	double xc = 0., yc = 0.;
	double mass = 0.;
	for (int y=y0-nb; y<=y0+nb; y++) {
		unsigned short *p = scanLine(y) + (x0-nb);
		for (int x=x0-nb; x<=x0+nb; x++) {
			double v = double(*p++);
			mass += v;
			xc += v * (x+0.5);
			yc += v * (y+0.5);
		}
	}
	if (mass < SEGM_EPS) {
		xc = yc = 0.;
	} else {
		xc = xc / mass - 0.5;
		yc = yc / mass - 0.5;
	}
	*p_xc = xc;
	*p_yc = yc;
}

bool Raster16::is_local_max(int x, int y, int nbsz)
{
	int vmax = int(value(x, y));
	if (vmax <= 0) return false;

	for (int j=0; j<nbsz; j++) {
		NbrPoint & npt = hood_pts[j];
		int cx = x + npt.dx;
		int cy = y + npt.dy;
		if (int(value(cx, cy)) > vmax)
			return false;
	}
	return true;
}

int Raster16::local_rank(int x, int y, int nbsz)
{
	int vcent = int(value(x, y));
	int rank = 0;
	for (int j=1; j<nbsz; j++) {
		if (int(value(x+hood_pts[j].dx, y+hood_pts[j].dy)) > vcent)
			++rank;
	}
	return rank;
}

void Raster16::paintParticleFill(std::vector<HSeg> &fill, unsigned short c)
{
	for (HSeg &hs : fill) {
		unsigned short *b = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) b[x] = c;
	}
}

void Raster16::countColors(std::vector<ColorCounter> &res, std::vector<HSeg> &fill)
{
	res.resize(0);
	for (HSeg &hs : fill) {
		unsigned short *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) {
			int c = (int) p[x];
			if (c == 0xFFFF) continue;
			bool found = false;
			for (ColorCounter &cc : res) {
				if (cc.c == c) {
					++cc.cnt;
					found = true;
					break;
				}
			}
			if (!found)
				res.push_back(ColorCounter(c, 1));
		}
	}
}

void Raster16::fill(Boundary &bnd, unsigned short c)
{
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned short *p = scanLine(y);
		for (int x=bnd.xmin; x<=bnd.xmax; x++)
			p[x] = c;
	}
}

void Raster16::paintParticleADD(Particle &ptc, int incr)
{
	unsigned short ic = (unsigned short)(incr);
	for (HSeg &hs : ptc.fill) {
		unsigned short *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) {
			p[x] += ic;
		}
	}
}

int Raster16::rescanParticleThresh(Particle &ptc, unsigned short tc)
{
	ptc.fill.resize(0);
	Boundary bnd = ptc.bnd;
	HSeg hs;
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned short *p = scanLine(y);
		bool is_in = false;
		hs.y = y;
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (p[x] >= tc) {
				if (!is_in) {
					hs.xl = x;
					is_in = true;
				}
				hs.xr = x;
			} else {
				if (is_in) {
					ptc.fill.push_back(hs);
					is_in = false;
				}
			}
		}
		if (is_in) ptc.fill.push_back(hs);
	}
	if (ptc.fill.size() > 0) {
		ptc.x0 = ptc.fill[0].xl;
		ptc.y0 = ptc.fill[0].y;
	} else {
		ptc.x0 = ptc.y0 = 0;
	}
	ptc.bnd = fill_boundary(ptc.fill);
	return fill_area(ptc.fill);
}

void Raster16::paintParticleInto(Particle &ptc, unsigned short c, unsigned short bkg)
{
	for (HSeg &hs : ptc.fill) {
		unsigned short *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) {
			if (p[x] == bkg) p[x] = c;
		}
	}
}

//----------------------------- Raster8 --------------------------------------

void Raster8::fill(Boundary &bnd, unsigned char c)
{
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned char *p = scanLine(y);
		for (int x=bnd.xmin; x<=bnd.xmax; x++)
			p[x] = c;
	}
}

void Raster8::fillBorder(unsigned char c, int bsz)
{
	int x, y;
	for (y=0; y<bsz; y++) {
		memset(scanLine(y), c, w);
		memset(scanLine(h-1-y), c, w);
	}
	for (y=bsz; y<h-bsz; y++) {
		unsigned char *p = scanLine(y);
		for (x=0; x<bsz; x++) {
			p[w-1-x] = p[x] = c;
		}
	}
}

void Raster8::paintBoundary(const Boundary& bnd, unsigned char c)
{
	int y;
	int w = bnd.xmax - bnd.xmin + 1;
	memset(scanLine(bnd.ymin)+bnd.xmin, c, w);
	for (y = bnd.ymin+1; y < bnd.ymax; y++) {
		unsigned char *p = scanLine(y);
		p[bnd.xmin] = p[bnd.xmax] = c;
	}
	memset(scanLine(bnd.ymax)+bnd.xmin, c, w);
}

void Raster8::replaceColor(unsigned char oldc, unsigned char newc)
{
	for (int y=0; y<h; y++) {
		unsigned char *p = scanLine(y);
		for (int x=0; x<w; x++) {
			if (p[x] == oldc) p[x] = newc;
		}
	}
}

void Raster8::replaceColor(Boundary &b, unsigned char oldc, unsigned char newc)
{
	for (int y=b.ymin; y<=b.ymax; y++) {
		unsigned char *p = scanLine(y);
		for (int x=b.xmin; x<=b.xmax; x++) {
			if (p[x] == oldc) p[x] = newc;
		}
	}
}

void Raster8::paintSegmentOut(Point pt1, Point pt2, unsigned char c, unsigned char cb)
{
	int dx = pt2.x - pt1.x;
	if (dx < 0) dx = -dx;
	int dy = pt2.y - pt1.y;
	if (dy < 0) dy = -dy;
	if (dx >= dy) {
		if (dx == 0) dx = 1;
		else if (pt1.x > pt2.x) pt1.swap(pt2);
		dy = pt2.y - pt1.y;
		setMaxValue(pt1.x-1, pt1.y, cb);
		for (int x=pt1.x; x<=pt2.x; x++) {
			int y = pt1.y + (x - pt1.x) * dy / dx;
			setMaxValue(x, y-1, cb);
			setValue(x, y, c);
			setMaxValue(x, y+1, cb);
		}
		setMaxValue(pt2.x+1, pt2.y, cb);
	} else {
		if (pt1.y > pt2.y) pt1.swap(pt2);
		dx = pt2.x - pt1.x;
		setMaxValue(pt1.x, pt1.y-1, cb);
		for (int y=pt1.y; y<=pt2.y; y++) {
			int x = pt1.x + (y - pt1.y) * dx / dy;
			setMaxValue(x-1, y, cb);
			setValue(x, y, c);
			setMaxValue(x+1, y, cb);
		}
		setMaxValue(pt2.x, pt2.y+1, cb);
	}
}

void Raster8::paintPath(std::vector<Point> &path, unsigned char c, unsigned char cb)
{
	for (size_t j=1; j<path.size(); j++) {
		paintSegmentOut(path[j-1], path[j], c, cb);
	}
}

void Raster8::paintParticleFill(std::vector<HSeg> &fill, unsigned char c)
{
	for (HSeg &hs : fill) {
		unsigned char *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) p[x] = c;
	}
}

void Raster8::paintParticleFillInto(std::vector<HSeg> &fill, unsigned char c, unsigned char bkg)
{
	for (HSeg &hs : fill) {
		unsigned char *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++)
			if (p[x] == bkg) p[x] = c;
	}
}

void Raster8::paintParticleFillADD(std::vector<HSeg> &fill, int incr)
{
	unsigned char c = (unsigned char)(incr);
	for (HSeg &hs : fill) {
		unsigned char *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) p[x] += c;
	}
}

void Raster8::paintParticlePat(Particle &ptc, unsigned char fg)
{
	for (HSeg &hs : ptc.fill) {
		unsigned char *p = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) {
			if ((hs.y+x) & 1) p[x] = fg;
		}
	}
}

void Raster8::paintSegment(Point pt1, Point pt2, unsigned char c)
{
	int dx = pt2.x - pt1.x;
	if (dx < 0) dx = -dx;
	int dy = pt2.y - pt1.y;
	if (dy < 0) dy = -dy;
	if (dx >= dy) {
		if (dx == 0) dx = 1;
		else if (pt1.x > pt2.x) pt1.swap(pt2);
		dy = pt2.y - pt1.y;
		for (int x=pt1.x; x<=pt2.x; x++) {
			int y = pt1.y + (x - pt1.x) * dy / dx;
			setValue(x, y, c);
		}
	} else {
		if (pt1.y > pt2.y) pt1.swap(pt2);
		dx = pt2.x - pt1.x;
		for (int y=pt1.y; y<=pt2.y; y++) {
			int x = pt1.x + (y - pt1.y) * dx / dy;
			setValue(x, y, c);
		}
	}
}

void Raster8::fillContour(Contour &cont, unsigned char fill_c, unsigned char bord_c)
{
	if (cont.is_closed) {
		for (int y=cont.bnd.ymin; y<=cont.bnd.ymax; y++) {
			unsigned char *p = scanLine(y);
			for (int x=cont.bnd.xmin; x<=cont.bnd.xmax; x++) {
				if (cont.IsInside(Point(x, y)))
					p[x] = fill_c;
			}
		}
	}
	for (size_t j=1; j<cont.path.size(); j++) {
		paintSegment(cont.path[j-1], cont.path[j], bord_c);
	}
}

void Raster8::paintContourBorder(Contour &cont, unsigned char bord_c)
{
	for (size_t j=1; j<cont.path.size(); j++) {
		paintSegment(cont.path[j-1], cont.path[j], bord_c);
	}
}

int Raster8::expandBorders(unsigned char oldc, unsigned char newc, int nbsz, unsigned char c0)
{
	int cnt = 0;
	int bsz = 1;
	if (nbsz > HOOD_SIZE_FATCROSS) bsz = 3;
	else if (nbsz > HOOD_SIZE_MOORE) bsz = 2;
	for (int y0=bsz; y0<h-bsz; y0++) {
		unsigned char *p = scanLine(y0);
		for (int x0=bsz; x0<w-bsz; x0++) {
			if (p[x0] != c0) continue;
			for (int j=1; j<nbsz; j++)
			{
				int x = x0+hood_pts[j].dx;
				int y = y0+hood_pts[j].dy;
				if (value(x, y) == oldc) {
					p[x0] = newc;
					++cnt;
					break;
				}
			}
		}
	}
	return cnt;
}

void Raster8::expandBordersInto(Boundary &bnd, unsigned char fg, unsigned char bk, unsigned char tmpc,
		int nbsz, bool keeptmpc)
{
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned char *p = scanLine(y);
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (p[x] != bk) continue;
			for (int j=1; j<nbsz; j++) {
				if (value(x+hood_pts[j].dx, y+hood_pts[j].dy) == fg) {
					p[x] = tmpc;
					break;
				}
			}
		}
	}
	if (!keeptmpc)
		replaceColor(bnd, tmpc, fg);
}

void Raster8::filterParticles(unsigned char oldc, unsigned char newc, int minsz, unsigned char bk)
{
	for (int y0=2; y0<h-2; y0++) {
		for (int x0=2; x0<w-2; x0++) {
			if (value(x0, y0) != oldc) continue;
			std::vector<HSeg> fill;
			findParticleFill(fill, x0, y0, newc);
			if (fill_area(fill) < minsz)
				paintParticleFill(fill, bk);
		}
	}
}

void Raster8::filterParticles(Boundary& bnd, unsigned char oldc, unsigned char newc, int minsz, unsigned char bk)
{
	for (int y0=bnd.ymin+1; y0<bnd.ymax; y0++) {
		for (int x0=bnd.xmin+1; x0<bnd.xmax; x0++) {
			if (value(x0, y0) != oldc) continue;
			std::vector<HSeg> fill;
			findParticleFill(fill, x0, y0, newc);
			if (fill_area(fill) < minsz)
				paintParticleFill(fill, bk);
		}
	}
}

bool Raster8::touches(int x, int y, unsigned char c, int nbsz)
{
	for (int j=1; j<nbsz; j++) {
		if (value(x + hood_pts[j].dx, y + hood_pts[j].dy) == c) return true;
	}
	return false;
}

int Raster8::countColors(int x0, int y0, int nbsz, unsigned char c)
{
	int cnt = 0;
	for (int j=0; j<nbsz; j++) {
		if (value(x0+hood_pts[j].dx, y0+hood_pts[j].dy) == c)
			++cnt;
	}
	return cnt;
}

int Raster8::countColorsNb(std::vector<Point> path, unsigned char c, int nbsz)
{
	if (path.size() == 0) return 0;
	size_t n = path.size();
	if (n>2 && path[0].equals(path[n-1])) --n;
	int cnt = 0;
	for (size_t i=0; i<n; i++) {
		Point &p = path[i];
		bool found_another = false;
		for (int j=0; j<nbsz; j++) {
			if (value(p.x+hood_pts[j].dx, p.y+hood_pts[j].dy) != c) {
				found_another = true;
				break;
			}
		}
		if (!found_another) ++cnt;
	}
	return cnt;
}

int Raster8::detectParticle(Particle &ptc, int x0, int y0, unsigned char newc)
{
	ptc.fill.resize(0);
	ptc.x0 = x0;
	ptc.y0 = y0;
	findParticleFill(ptc.fill, ptc.x0, ptc.y0, newc);
	ptc.bnd = fill_boundary(ptc.fill);
	return fill_area(ptc.fill);
}

void Raster8::shrinkParticle(Particle &ptc, unsigned char fillc, unsigned char tmpc,
		unsigned char bordc, int nbsz)
{
	int bsz = 1;
	if (nbsz > HOOD_SIZE_FATCROSS) bsz = 3;
	else if (nbsz > HOOD_SIZE_MOORE) bsz = 2;
	paintParticle(ptc, tmpc);
	Boundary &bnd = ptc.bnd;
	int x, y;
	for (y=bnd.ymin; y<bnd.ymin+bsz; y++) {
		unsigned char *p = scanLine(y);
		for (x=bnd.xmin; x<=bnd.xmax; x++)
			if (p[x] == tmpc) p[x] = bordc;
	}
	for (y=bnd.ymin+bsz; y<=bnd.ymax-bsz; y++) {
		unsigned char *p = scanLine(y);
		for (x=bnd.xmin; x<bnd.xmin+bsz; x++)
			if (p[x] == tmpc) p[x] = bordc;
		for (x=bnd.xmax-bsz+1; x<=bnd.xmax; x++)
			if (p[x] == tmpc) p[x] = bordc;
	}
	for (y=bnd.ymax-bsz+1; y<=bnd.ymax; y++) {
		unsigned char *p = scanLine(y);
		for (x=bnd.xmin; x<=bnd.xmax; x++)
			if (p[x] == tmpc) p[x] = bordc;
	}
	for (y=bnd.ymin+bsz; y<=bnd.ymax-bsz; y++) {
		unsigned char *p = scanLine(y);
		for (x=bnd.xmin+bsz; x<=bnd.xmax-bsz; x++) {
			if (p[x] != tmpc) continue;
			for (int j=1; j<nbsz; j++) {
				unsigned char c = value(x+hood_pts[j].dx, y+hood_pts[j].dy);
				if (c != tmpc && c != bordc) {
					p[x] = bordc;
					break;
				}
			}
		}
	}
	rescanParticle(ptc, tmpc);
	paintParticle(ptc, fillc);
}

void Raster8::expandParticle(Particle &ptc, unsigned char fillc, unsigned char tmpc,
		unsigned char bordc, int nbsz)
{
	paintParticle(ptc, tmpc);
	Boundary &bnd = ptc.bnd;
	bnd.expand(2);
	clip(bnd);
	std::vector<Point> expts;
	for (int y=bnd.ymin+1; y<=bnd.ymax-1; y++) {
		unsigned char *p = scanLine(y);
		for (int x=bnd.xmin+1; x<=bnd.xmax-1; x++) {
			if (p[x] != bordc) continue;
			
			bool touch_this = false;
			bool touch_other = false;
			for (int j=1; j<HOOD_SIZE_MOORE; j++) {
				unsigned char c = value(x+hood_pts[j].dx, y+hood_pts[j].dy);
				if (c == fillc) touch_other = true;
				else if (c == tmpc && j<nbsz) touch_this = true;
			}
			if (!touch_this) continue;
			if (touch_other) p[x] = forcedc;
			else expts.push_back(Point(x, y));
		}
	}
	for (Point &pt : expts) setValue(pt.x, pt.y, tmpc);
	rescanParticle(ptc, tmpc);
	paintParticle(ptc, fillc);
}

int Raster8::countColor(std::vector<HSeg>& fill, unsigned char c)
{
	int cnt = 0;
	for (HSeg &hs : fill) {
		unsigned char *b = scanLine(hs.y);
		for (int x=hs.xl; x<=hs.xr; x++) {
			if (b[x] == c) ++cnt;
		}
	}
	return cnt;
}

void Raster8::rescanParticleFill(Boundary &bnd, std::vector<HSeg> &fill, unsigned char c)
{
	fill.resize(0);
	HSeg hs;
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned char *p = scanLine(y);
		bool is_in = false;
		hs.y = y;
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (p[x] == c) {
				if (!is_in) {
					hs.xl = x;
					is_in = true;
				}
				hs.xr = x;
			} else {
				if (is_in) {
					fill.push_back(hs);
					is_in = false;
				}
			}
		}
		if (is_in) fill.push_back(hs);
	}
}

int Raster8::rescanParticle(Particle &ptc, unsigned char c)
{
	rescanParticleFill(ptc.bnd, ptc.fill, c);
	if (ptc.fill.size() > 0) {
		ptc.x0 = ptc.fill[0].xl;
		ptc.y0 = ptc.fill[0].y;
	} else {
		ptc.x0 = ptc.y0 = 0;
	}
	ptc.bnd = fill_boundary(ptc.fill);
	return fill_area(ptc.fill);
}

int Raster8::removeFalseWalls(Particle &ptc, unsigned char fg, unsigned char bk, unsigned char tmpc)
{
	unsigned char tmpc1 = tmpc + 1;
	paintParticle(ptc, tmpc);
	Boundary bnd = ptc.bnd;
	bnd.expand(1);
	clip(bnd);
	expandBordersInto(bnd, tmpc, bk, tmpc1);
	expandBordersInto(bnd, bk, tmpc, tmpc1);
	ptc.bnd = bnd;
	int rc = rescanParticle(ptc, tmpc);
	paintParticle(ptc, fg);
	return rc;
}

int Raster8::mergeParticles(Particle &ptc, Particle &other,
		unsigned char fgd, unsigned char bkg, unsigned char tmpc)
{
	unsigned char tmp_fg1 = tmpc;
	unsigned char tmp_bord = tmpc + 1;
	unsigned char tmp_fg2 = tmpc + 2;

	Boundary bnd = ptc.bnd;
	bnd.combo(other.bnd);

	paintParticle(ptc, tmp_fg1);
	paintParticle(other, tmp_fg2);
	other.clear();
	
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned char *p = scanLine(y);
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (p[x] != bkg) continue;
			bool t1=false, t2=false;
			for (int j=1; j<HOOD_SIZE_MOORE; j++) {
				int y0 = y + hood_pts[j].dy;
				int x0 = x + hood_pts[j].dx;
				if (!bnd.IsInside(x0, y0)) continue;
				unsigned char c = value(x0, y0);
				if (c == tmp_fg1) {
					t1 = true;
					if (t2) break;
				} else if (c == tmp_fg2) {
					t2 = true;
					if (t1) break;
				}
			}
			if (t1 && t2) p[x] = tmp_bord;
		}
	}
	replaceColor(bnd, tmp_bord, tmp_fg1);
	replaceColor(bnd, tmp_fg2, tmp_fg1);
	
	ptc.bnd = bnd;
	int a = rescanParticle(ptc, tmp_fg1);
	paintParticle(ptc, fgd);
	return a;
}


void Raster8::findFillSegments(std::vector<HSeg> &fill, int y, int xl, int xr,
		unsigned char oldc, unsigned char newc)
{
	unsigned char *p = scanLine(y);
	if (p[xl] == oldc) {
		for (; xl>=0; xl--) if (p[xl] != oldc) break;
		++xl;
	}
	if (p[xr] == oldc) {
		for (; xr<w; xr++) if (p[xr] != oldc) break;
		--xr;
	}
	HSeg hs;
	hs.y = y;
	bool is_in = false;
	for (int x=xl; x<=xr; x++) {
		if (p[x] == oldc) {
			if (!is_in) {
				is_in = true;
				hs.xl = x;
			}
			hs.xr = x;
		} else {
			if (is_in) {
				is_in = false;
				paintHSeg(hs, newc);
				fill.push_back(hs);
			}
		}
	}
	if (is_in) {
		paintHSeg(hs, newc);
		fill.push_back(hs);
	}
}

void Raster8::findParticleFill(std::vector<HSeg> &fill, int x0, int y0, unsigned char newc)
{
	unsigned char *p = scanLine(y0);
	unsigned char oldc = p[x0];
	int xl, xr;
	for (xl=x0; xl>=0; xl--)
		if (p[xl] != oldc) break;
	++xl;
	for (xr=x0; xr<w; xr++)
		if (p[xr] != oldc) break;
	--xr;
	
	HSeg hs = HSeg(y0, xl, xr);
	paintHSeg(hs, newc);

	size_t idx = fill.size();
	fill.push_back(hs);
	while(idx < fill.size()) {
		hs = fill[idx];
		++idx;
		if (hs.y > 0) {
			findFillSegments(fill, hs.y - 1, hs.xl, hs.xr, oldc, newc);
		}
		if (hs.y < h-1) {
			findFillSegments(fill, hs.y + 1, hs.xl, hs.xr, oldc, newc);
		}
	}
}

void Raster8::bordersAround(Boundary &bnd, unsigned char fg, unsigned char bk, unsigned char bordc,
		int hoodsz)
{
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned char *p = scanLine(y);
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (p[x] != bk) continue;
			for (int j=1; j<hoodsz; j++) {
				if (value(x+hood_pts[j].dx, y+hood_pts[j].dy) == fg) {
					p[x] = bordc;
					break;
				}
			}
		}
	}
}


//----------------------------- Raster3D --------------------------------------

void Raster3D::replaceColor(Boundary3D& bnd, unsigned char oldc, unsigned char newc)
{
	Boundary b2 = bnd.boundary2d();
	for (int z=bnd.zmin; z<=bnd.zmax; z++) {
		Raster8 msk = getPlane(z);
		msk.replaceColor(b2, oldc, newc);
	}
}

void Raster3D::expandBorders(Boundary3D& bnd, unsigned char fg, unsigned char bk, unsigned char bordc, int nbsz)
{
	std::vector<NbrPoint3D> nb3d = hood_sorted_by_z(nbsz);
	for (int z=bnd.zmin; z<=bnd.zmax; z++) {
		int j0 = -1, j1;
		for (int j=0; size_t(j)<nb3d.size(); j++) {
			int zz = z + nb3d[j].dz;
			if (zz < 0 || zz >= d) continue;
			if (j0 < 0) j0 = j;
			j1 = j;
		}
		
		for (int y=bnd.ymin; y<=bnd.ymax; y++) {
			unsigned char *p = scanLine(y, z);
			for (int x=bnd.xmin; x<=bnd.xmax; x++) {
				if (p[x] != bk) continue;
				for (int j=j0; j<=j1; j++) {
					NbrPoint3D& nbp = nb3d[j];
					if (value(x + nbp.dx, y + nbp.dy, z + nbp.dz) == fg) {
						p[x] = bordc;
						break;
					}
				}
			}
		}
	}
}


void Raster3D::chopBorders(unsigned char fg, unsigned char bk, unsigned char bordc, int nbsz)
{
	Boundary3D bnd = getBoundary();
	for (int z=bnd.zmin; z<=bnd.zmax; z++) {
		for (int y=bnd.ymin; y<=bnd.ymax; y++) {
			unsigned char *p = scanLine(y, z);
			for (int x=bnd.xmin; x<=bnd.xmax; x++) {
				if (p[x] != bordc) continue;
				bool found = false;
				for (int j=1; j<nbsz; j++) {
					int x0 = x + hood3d_pts[j].dx;
					int y0 = y + hood3d_pts[j].dy;
					int z0 = z + hood3d_pts[j].dz;
					if (!bnd.IsInside(x0, y0, z0)) continue;
					if (value(x0, y0, z0) == fg) {
						found = true;
						break;
					}
				}
				if (!found) p[x] = bk;
			}
		}
	}
}

void Raster3D::paintParticle(Particle3D &pt, unsigned char c)
{
	for (int z=0; size_t(z)<pt.fills.size(); z++) {
		for (HSeg &hs : pt.fills[z]) {
			unsigned char *p = scanLine(hs.y, z);
			for (int x=hs.xl; x<=hs.xr; x++)
				p[x] = c;
		}
	}
}

void Raster3D::paintParticleADD(Particle3D &pt, int incr)
{
	for (int z=0; size_t(z)<pt.fills.size(); z++) {
		for (HSeg &hs : pt.fills[z]) {
			unsigned char *p = scanLine(hs.y, z);
			for (int x=hs.xl; x<=hs.xr; x++) {
				int v = int(p[x]) + incr;
				//if (v < 0) v = 0;
				//else if (v > 0xFF) v = 0xFF;
				p[x] = (unsigned char)(v);
			}
		}
	}
}

void Raster3D::paintParticleInto(Particle3D &pt, unsigned char newc, unsigned char oldc)
{
	for (int z=0; size_t(z)<pt.fills.size(); z++) {
		for (HSeg &hs : pt.fills[z]) {
			unsigned char *p = scanLine(hs.y, z);
			for (int x=hs.xl; x<=hs.xr; x++)
				if (p[x] == oldc) p[x] = newc;
		}
	}
}

void Raster3D::forceSeparation(Particle3D &pt, unsigned char fg, unsigned char bk, unsigned char tmpc)
{
	int nbsz = HOOD3D_26;
	Boundary3D bnd = getBoundary();
	paintParticle(pt, tmpc);
	for (int z=0; size_t(z)<pt.fills.size(); z++) {
		for (HSeg &hs : pt.fills[z]) {
			unsigned char *p = scanLine(hs.y, z);
			for (int x=hs.xl; x<=hs.xr; x++) {
				for (int j=1; j<nbsz; j++) {
					int x0 = x + hood3d_pts[j].dx;
					int y0 = hs.y + hood3d_pts[j].dy;
					int z0 = z + hood3d_pts[j].dz;
					if (!bnd.IsInside(x0, y0, z0)) continue;
					if (value(x0, y0, z0) == fg) {
						p[x] = bk;
						break;
					}
				}
			}
		}
	}
	paintParticleInto(pt, fg, tmpc);
}

void Raster3D::forceBorder(Particle3D& pt, unsigned char fg, unsigned char bordc,
		unsigned char tmp_fg, unsigned char tmp_bordc)
{
	paintParticle(pt, tmp_fg);
	for (int z0=0; z0<d; z0++) {
		if (pt.fills[z0].empty()) continue;
		Raster8 msk = getPlane(z0);
		for (int y0=pt.bnd.ymin; y0<=pt.bnd.ymax; y0++) {
			unsigned char *p = msk.scanLine(y0);
			for (int x0=pt.bnd.xmin; x0<=pt.bnd.xmax; x0++) {
				if (p[x0] != tmp_fg) continue;
				bool at_bord = false;
				for (int z=z0-1; z<=z0+1; z++) {
					if (z < 0 || z >= d) continue;
					for (int j=1; j<HOOD_SIZE_MOORE; j++) {
						unsigned char c = value(x0+hood_pts[j].dx, y0+hood_pts[j].dy, z);
						if (c!=tmp_fg && c!=tmp_bordc) {
							p[x0] = tmp_bordc;
							at_bord = true;
							break;
						}
					}
					if (at_bord) break;
				}
			}
		}
	}
	for (int z=0; z<d; z++) {
		if (pt.fills[z].empty()) continue;
		Raster8 msk = getPlane(z);
		for (int y=pt.bnd.ymin; y<=pt.bnd.ymax; y++) {
			unsigned char *p = msk.scanLine(y);
			for (int x=pt.bnd.xmin; x<=pt.bnd.xmax; x++) {
				if (p[x] == tmp_fg) p[x] = fg;
				else if (p[x] == tmp_bordc) p[x] = bordc;
			}
		}
	}
}

void Raster3D::rescanShrunkParticles(std::vector<Particle3D> &cells, unsigned char fg, unsigned char tmpc)
{
	for (int z=0; z<d; z++) {
		Raster8 msk = getPlane(z);
		for (Particle3D &cell : cells) {
			Boundary bnd = cell.bnd.boundary2d();
			std::vector<HSeg> & fill = cell.fills[z];
			if (fill.empty()) continue;
			msk.paintParticleFillInto(fill, tmpc, fg);
			msk.rescanParticleFill(bnd, fill, tmpc);
			fill.shrink_to_fit();
			msk.paintParticleFill(fill, fg);
		}
	}
	for (Particle3D &pt : cells) {
		pt.update_from_fill();
	}
}

void Raster3D::rescanShrunkCells(std::vector<Cell> &cells, unsigned char fg, unsigned char tmpc)
{
	for (int z=0; z<d; z++) {
		Raster8 msk = getPlane(z);
		for (Cell &cell : cells) {
			Boundary bnd = cell.bnd.boundary2d();
			std::vector<HSeg> & fill = cell.fills[z];
			if (fill.empty()) continue;
			msk.paintParticleFillInto(fill, tmpc, fg);
			msk.rescanParticleFill(bnd, fill, tmpc);
			fill.shrink_to_fit();
			msk.paintParticleFill(fill, fg);
		}
	}
	for (Cell &pt : cells) {
		pt.update_from_fill();
	}
}

void Raster3D::sandPaperCells(std::vector<Particle3D> &cells)
{
	fill(0x80);
	for (Particle3D& cell : cells)
		paintParticle(cell, 0);
	expandBorders(0x80, 0, 0xFF, HOOD3D_26);
	
	int j=0;
	for (int i=0; size_t(i)<cells.size(); i++) {
		Particle3D& scell = cells[i];
		Particle3D& tcell = cells[j];
		Boundary3D bnd = scell.bnd;
		paintParticleInto(scell, 0x50, 0);
		
		tcell = findBiggestCell(bnd, 0x50, 0xFF, 0x55);
		paintParticle(tcell, 0x50);
		replaceColor(bnd, 0x55, 0xFF);
		expandBorders(bnd, 0x50, 0xFF, 0x55, HOOD3D_26);
		replaceColor(bnd, 0x55, 0x50);
		tcell.bnd = bnd;
		rescanParticle(tcell, 0x50);
		paintParticle(tcell, 0);
		
		if (tcell.realHeight() >= 3) ++j;
	}
	
	cells.resize(j);
}

long long Raster3D::rescanParticle(Particle3D& cell, unsigned char c)
{
	for (std::vector<HSeg>& fill : cell.fills)
		fill.clear();
	
	Boundary b2 = cell.bnd.boundary2d();
	for (int z=cell.bnd.zmin; z<=cell.bnd.zmax; z++) {
		Raster8 msk = getPlane(z);
		msk.rescanParticleFill(b2, cell.fills[z], c);
	}
	
	return cell.update_from_fill();
}

long long Raster3D::findParticleFillsZ(std::vector<std::pair<int, std::vector<HSeg>>>& zfills,
		int x0, int y0, int z0, unsigned char newc)
{
	unsigned char oldc = value(x0, y0, z0);
	Raster8 msk = getPlane(z0);
	
	size_t idx = zfills.size();
	zfills.resize(zfills.size() + 1);
	zfills[zfills.size()-1].first = z0;
	msk.findParticleFill(zfills[zfills.size()-1].second, x0, y0, newc);
	
	while(idx < zfills.size()) {
		int z = zfills[idx].first;
		std::vector<HSeg> fill = zfills[idx].second;
		++idx;
		findFillZSlices(zfills, z, fill, oldc, newc);
	}
	
	long long vol = 0;
	
	for (auto& el : zfills) {
		for (HSeg hs : el.second) {
			vol += (hs.xr - hs.xl + 1);
		}
	}
	
	return vol;
}

void Raster3D::findFillZSlices(std::vector<std::pair<int, std::vector<HSeg>>>& zfills,
			int z0, std::vector<HSeg>& fill, unsigned char oldc, unsigned char newc)
{
	std::vector<int> zs;
	if (z0 > 0) zs.push_back(z0-1);
	if (z0+1 < d) zs.push_back(z0+1);
	for (int z : zs) {
		Raster8 msk = getPlane(z);
		for (HSeg hs : fill) {
			unsigned char *p = msk.scanLine(hs.y);
			for (int x=hs.xl; x<=hs.xr; x++) {
				if (p[x] != oldc) continue;
				
				zfills.resize(zfills.size()+1);
				zfills[zfills.size()-1].first = z;
				msk.findParticleFill(zfills[zfills.size()-1].second, x, hs.y, newc);
			}
		}
	}
}

Particle3D Raster3D::findBiggestCell(Boundary3D& bnd, unsigned char fg, unsigned char bk, unsigned char tmpc)
{
	Particle3D cell;
	cell.fills.resize(d);
	long long vol = 0;
	
	for (int z=bnd.zmin; z<=bnd.zmax; z++) {
		for (int y=bnd.ymin; y<=bnd.ymax; y++) {
			unsigned char *p = scanLine(y, z);
			for (int x=bnd.xmin; x<=bnd.xmax; x++) {
				if (p[x] != fg) continue;
				
				std::vector<std::pair<int, std::vector<HSeg>>> zfills;
				long long zvol = findParticleFillsZ(zfills, x, y, z, tmpc);
				if (zvol < vol) continue;
				
				cell.clear();
				for (auto& el : zfills) {
					std::vector<HSeg>& fill = cell.fills[el.first];
					for (HSeg hs : el.second)
						fill.push_back(hs);
				}
				vol = zvol;
			}
		}
	}
	
	cell.update_from_fill();
	return cell;
}

//----------------------------- Global Utilities --------------------------------------

static bool cmp_nbr_pt_3d_by_z(const NbrPoint3D& a, const NbrPoint3D& b)
{
	if (a.dz < b.dz) return true;
	return false;
};
std::vector<NbrPoint3D> hood_sorted_by_z(int nbsz)
{
	std::vector<NbrPoint3D> res;
	
	for (int j=0; j<nbsz; j++)
		res.push_back(hood3d_pts[j]);
	
	std::sort(res.begin(), res.end(), cmp_nbr_pt_3d_by_z);
	
	return res;
}

