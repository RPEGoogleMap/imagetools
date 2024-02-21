#define __geom_main__
#include "geom.h"

//----------------------------- Contour --------------------------------

Contour::Contour(const std::vector<int> & flat_cont, int iStart, int len)
{
	if (iStart < 0) {
		len = flat_cont[0];
		iStart = 1;
	}
	is_closed = false;
	path = path_from_flat(flat_cont, iStart, len);
	if (path.size() > 3) {
		path.push_back(path[0]);
		is_closed = true;
	}
	_update_from_path();
}

void Contour::set_path(std::vector<Point> &_path)
{
	is_closed = false;
	path.resize(_path.size());
	if (path.size() > 0) {
		for (size_t j=0; j<path.size(); j++)
			path[j] = _path[j];
		if (path.size() > 3 && path[0].equals(path[path.size()-1]))
			is_closed = true;
	}
	_update_from_path();
}

void Contour::_update_from_path()
{
	plen = 0.;
	if (path.size() > 0) {
		for(size_t j=1; j<path.size(); j++) {
			plen += path[j-1].dist(path[j]);
		}
		bnd.xmin = bnd.xmax = path[0].x;
		bnd.ymin = bnd.ymax = path[0].y;
		for (Point &p : path) {
			if (bnd.xmin > p.x) bnd.xmin = p.x;
			if (bnd.xmax < p.x) bnd.xmax = p.x;
			if (bnd.ymin > p.y) bnd.ymin = p.y;
			if (bnd.ymax < p.y) bnd.ymax = p.y;
		}
	} else {
		bnd.xmin = bnd.xmax = bnd.ymin = bnd.ymax = 0;
	}
}

void Contour::optimize(double mindist, double maxdist)
{
	if (path.size() < 4) return;
	std::vector<Point> npath;
	npath.reserve(path.size());
	int iLow=0, iHigh=1;
	while (size_t(iHigh) < path.size()) {
		Point &p0 = path[iLow];
		npath.push_back(p0);
		for (; size_t(iHigh) < path.size(); iHigh++) {
			if (p0.dist(path[iHigh]) >= mindist) break;
		}
		while (p0.dist(path[iHigh]) > maxdist) {
			if (iHigh-1 <= iLow) break;
			--iHigh;
		}
		if (iLow == iHigh) break;
		iLow = iHigh;
		++iHigh;
	}
	if (iLow < iHigh && size_t(iHigh) < path.size()) {
		npath.push_back(path[iHigh]);
	}
	if (is_closed && !npath[0].equals(npath[npath.size()-1])) {
		npath.push_back(npath[0]);
	}
	set_path(npath);
}

bool Contour::low_trend(Point *pp1, Point *pp2, double mindist)
{
	if (path.size() == 0) return false;
	Point &p0 = path[0];
	Point p1;
	for (size_t j=1; j<path.size(); j++) {
		p1 = path[j];
		if (p0.dist(p1) >= mindist) break;
	}
	*pp1 = p0;
	*pp2 = p1;
	return true;
}

bool Contour::hi_trend(Point *pp3, Point *pp4, double mindist)
{
	if (path.size() == 0) return false;
	Point &p0 = path[path.size()-1];
	Point p1;
	for (size_t j=path.size()-1; j>0; j--) {
		p1 = path[j-1];
		if (p0.dist(p1) >= mindist) break;
	}
	*pp3 = p1;
	*pp4 = p0;
	return true;
}

/*
         P5  P6
        :      :
      P2         P3
     /             \
   P1               P4
    +---------------+
*/

bool Contour::can_close(double maxdist, double mindist)
{
	if (path.size() < 4) return false;
	double dist = path[0].dist(path[path.size()-1]);
	if (dist <= mindist) return true;
	if (dist > maxdist) return false;
	Point p1, p2, p3, p4;
	if (!low_trend(&p1, &p2)) return false;
	if (!hi_trend(&p3, &p4)) return false;
	double gap = p2.dist(p3);
	if (p1.dist(p4) <= gap) return false;

	Point p5(2*p2.x-p1.x, 2*p2.y-p1.y);
	Point p6(2*p3.x-p4.x, 2*p3.y-p4.y);
	if (p5.dist(p6) >= gap) return false;
	
	p3.x = p2.x + p4.x - p3.x;
	p3.y = p2.y + p4.y - p3.y;
	double a = p1.dist(p2);
	double b = p2.dist(p3);
	double c = p1.dist(p3);
	return c*c*1.2 > a*a + b*b;
}

//-- "Winding Number" algorithm for the inclusion of a point in polygon, adapted from:
// http://geomalgorithms.com/a03-inclusion.html

// tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
int Contour::isLeft(Point P0, Point P1, Point P2)
{
	return ((P1.x - P0.x) * (P2.y - P0.y)
		- (P2.x - P0.x) * (P1.y - P0.y));
}
bool Contour::IsInside(Point P)
{
	if (!is_closed || path.size() < 4) return false;
	int n = int(path.size() - 1);
	std::vector<Point> &V = path;

	int wn = 0;    // the  winding number counter (=0 only when P is outside)

	// loop through all edges of the polygon
	for (int i = 0; i < n; i++) {		// edge from V[i] to V[i+1]
		if (V[i].y <= P.y) {			// start y <= P.y
			if (V[i + 1].y > P.y)		// an upward crossing
				if (isLeft(V[i], V[i + 1], P) > 0)	// P left of  edge
					++wn;				// have  a valid up intersect
		}
		else {							// start y > P.y (no test needed)
			if (V[i + 1].y <= P.y)		// a downward crossing
				if (isLeft(V[i], V[i + 1], P) < 0)	// P right of  edge
					--wn;				// have  a valid down intersect
		}
	}
	return wn != 0;
}
//-- End of "Winding Number" algorithm


//----------------------------- Particle --------------------------------

void Particle::fromContour(Contour &cont)
{
	fill.resize(0);
	bnd = cont.bnd;
	HSeg hs;
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		bool is_in = false;
		hs.y = y;
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (cont.IsInside(Point(x, y))) {
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
	if (fill.size() > 0) {
		x0 = fill[0].xl;
		y0 = fill[0].y;
	} else {
		x0 = y0 = 0;
	}
}

long long Particle::update_from_fill()
{
	if (fill.size() > 0) {
		x0 = fill[0].xl;
		y0 = fill[0].y;
	} else {
		x0 = y0 = 0;
	}
	bnd = fill_boundary(fill);
	return fill_area(fill);
}

bool Particle::IsInside(int x, int y)
{
	if (!bnd.IsInside(x, y)) return false;
	for (HSeg &hs : fill) {
		if (hs.y != y) continue;
		if (x>=hs.xl && x<=hs.xr) return true;
	}
	return false;
}

Point Particle::center_mass()
{
	double xc = 0., yc = 0.;
	int npt = 0;
	for (HSeg &hs : fill) {
		int y = hs.y;
		for (int x=hs.xl; x<=hs.xr; x++) {
			xc += x;
			yc += y;
			++npt;
		}
	}
	if (npt == 0) return Point(0, 0);
	return Point(int(xc/npt), int(yc/npt));
}

long long Particle::overlay_area(std::vector<HSeg> &other_fill)
{
	return fill_overlay_area(fill, other_fill);
}
long long Particle::overlay_area(Particle &other)
{
	if (!bnd.intersects(other.bnd)) return 0;
	return fill_overlay_area(fill, other.fill);
}

long long ParticlePerim::sqdist(ParticlePerim& other)
{
	long long res = -1ll;
	for (Point& pt : perim) {
		for (Point& opt : other.perim) {
			int dx = pt.x - opt.x;
			int dy = pt.y - opt.y;
			long long qsd = (long long)dx * dx + (long long)dy * dy;
			if (res < 0 || qsd < res) res = qsd;
		}
	}
	return res;
}

//----------------------------- Particle3D --------------------------------

long long Particle3D::update_from_fill()
{
	long long vol = 0;
	x0 = y0 = z0 = 0;
	bnd = Boundary3D(0,0,0,0,0,0);
	bool is_first = true;
	for (int z=0; size_t(z)<fills.size(); z++) {
		std::vector<HSeg> &fill = fills[z];
		if (fill.empty()) continue;
		bnd.zmax = z;
		for (HSeg &hs : fill) {
			if (is_first) {
				is_first = false;
				x0 = hs.xl;
				y0 = hs.y;
				z0 = z;
				bnd.xmin = hs.xl;
				bnd.xmax = hs.xr;
				bnd.ymin = bnd.ymax = hs.y;
				bnd.zmin = z;
			} else {
				if (bnd.xmin > hs.xl) bnd.xmin = hs.xl;
				if (bnd.xmax < hs.xr) bnd.xmax = hs.xr;
				if (bnd.ymin > hs.y) bnd.ymin = hs.y;
				if (bnd.ymax < hs.y) bnd.ymax = hs.y;
			}
			vol += (hs.xr - hs.xl + 1);
		}
	}
	return vol;
}

long long Particle3D::overlay_volume(Particle3D &other)
{
	long long vol = 0;
	
	for (int z=0; size_t(z)<fills.size(); z++) {
		if (size_t(z) >= other.fills.size()) break;
		vol += fill_overlay_area(fills[z], other.fills[z]);
	}
	
	return vol;
}

double Particle3D::iou_score(int max_gap)
{
	double sc = 0.;
	int d = int(fills.size());
	for (int z1=0; z1<d-1; z1++) {
		std::vector<HSeg> &fill1 = fills[z1];
		if (fill1.empty()) continue;
		long long a1 = fill_area(fill1);
		for (int z2=z1+1; z2<d; z2++) {
			if (z2 - z1 > max_gap) break;
			std::vector<HSeg> &fill2 = fills[z2];
			if (fill2.empty()) continue;
			long long a2 = fill_area(fill2);
			long long ovl = fill_overlay_area(fill1, fill2);
			if (ovl > 0)
				sc += double(ovl) / (a1 + a2 - ovl);
		}
	}
	return sc;
}

Point3D Particle3D::center_mass()
{
	double xc=0., yc=0., zc=0.;
	long long cnt = 0;
	for (int z=0; size_t(z)<fills.size(); z++) {
		for (HSeg& hs : fills[z]) {
			for (int x=hs.xl; x<=hs.xr; x++) {
				zc += z;
				yc += hs.y;
				xc += x;
				++cnt;
			}
		}
	}
	if (cnt > 0) {
		xc /= cnt;
		yc /= cnt;
		zc /= cnt;
	}
	return Point3D(int(xc), int(yc), int(zc));
}

//----------------------------- Cell --------------------------------

long long Cell::update_from_fill()
{
	vol = 0;
	x0 = y0 = z0 = 0;
	bnd = Boundary3D(0,0,0,0,0,0);
	bool is_first = true;
	
	for (int z=0; size_t(z)<fills.size(); z++) {
		std::vector<HSeg> &fill = fills[z];
		if (fill.empty()) continue;
		bnd.zmax = z;
		vol += fill_area(fill);
		Boundary b = fill_boundary(fill);
		if (is_first) {
			is_first = false;
			bnd.zmin = z;
			bnd.set2d(b);
			HSeg &hs = fill[0];
			x0 = hs.xl;
			y0 = hs.y;
			z0 = z;
		} else {
			bnd.combo2d(b);
		}
	}
	return vol;
}

//----------------------------- Global utility functions --------------------------------

//--- Line/Circle/Ellipse drawings based on: http://members.chello.at/~easyfilter/Bresenham.pdf
std::vector<Point> straight_line_path(int x0, int y0, int x1, int y1)
{
	std::vector<Point> res;
	
	int dx = x1 - x0;
	int sx = 1;
	if (dx < 0) {
		dx = -dx;
		sx = -1;
	}
	int dy = y1 - y0;
	int sy = 1;
	if (dy < 0) sy = -1;
	else dy = -dy;
	int err = dx + dy;
	int err2;
	
	while (true) {
		res.push_back(Point(x0, y0));
		err2 = err + err;
		if (err2 >= dy) {
			if (x0 == x1) break;
			err += dy;
			x0 += sx;
		}
		if (err2 <= dx) {
			if (y0 == y1) break;
			err += dx;
			y0 += sy;
		}
	}
	
	return res;
}

std::vector<Point> circle_path(int xm, int ym, int r)
{
	std::vector<Point> res;
	std::vector<Point> q1;
	std::vector<Point> q2;
	std::vector<Point> q3;
	std::vector<Point> q4;
	
	int x = -r;
	int y = 0;
	int err = 2 - (r+r);
	do {
		q1.push_back(Point(xm-x, ym+y));
		q2.push_back(Point(xm-y, ym-x));
		q3.push_back(Point(xm+x, ym-y));
		q4.push_back(Point(xm+y, ym+x));
		r = err;
		if (r <= y) err += ++y*2+1;
		if (r > x || err > y)
			err += ++x*2+1;
	} while (x < 0);
	
	res.insert(res.end(), q1.begin(), q1.end());
	res.insert(res.end(), q2.begin(), q2.end());
	res.insert(res.end(), q3.begin(), q3.end());
	res.insert(res.end(), q4.begin(), q4.end());
	return res;
}

std::vector<HSeg> circle_fill(int xm, int ym, int r)
{
	std::vector<HSeg> res;
	std::vector<Point> pts;

	int x = -r;
	int y = 0;
	int err = 2 - (r+r);
	do {
		pts.push_back(Point(xm-x, ym+y));
		pts.push_back(Point(xm-y, ym-x));
		pts.push_back(Point(xm+x, ym-y));
		pts.push_back(Point(xm+y, ym+x));
		r = err;
		if (r <= y) err += ++y*2+1;
		if (r > x || err > y)
			err += ++x*2+1;
	} while (x < 0);
	
	std::sort(pts.begin(), pts.end(), [](Point &a, Point &b) {
		return a.less_than(b);   
	});
	
	if (!pts.empty()) {
		Point& p0 = pts[0];
		HSeg hs(p0.y, p0.x, p0.x);
		for (Point& p : pts) {
			if (p.y != hs.y) {
				res.push_back(hs);
				hs = HSeg(p.y, p.x, p.x);
			} else {
				if (p.x < hs.xl) hs.xl = p.x;
				if (p.x > hs.xr) hs.xr = p.x;
			}
		}
		res.push_back(hs);
	}
	
	return res;
}

std::vector<Point> ellipse_path(int xm, int ym, int a, int b)
{
	std::vector<Point> res;
	std::vector<Point> q1;
	std::vector<Point> q2;
	std::vector<Point> q3;
	std::vector<Point> q4;
	
	long long x = -a;
	long long y = 0;
	long long e2 = b;
	long long dx = (1+2*x)*e2*e2;
	long long dy = x*x;
	long long err = dx+dy;
	
	do {
		q1.push_back(Point(int(xm-x), int(ym+y)));
		q2.push_back(Point(int(xm+x), int(ym+y)));
		q3.push_back(Point(int(xm+x), int(ym-y)));
		q4.push_back(Point(int(xm-x), int(ym-y)));
		e2 = 2*err;
		if (e2 >= dx) {
			x++;
			err += dx += 2*(long)b*b;
		}
		if (e2 <= dy) {
			y++;
			err += dy += 2*(long)a*a;
		}
	} while (x <= 0);
	
	if (q1.back().equals(q2.back()))
		q2.resize(q2.size()-1);
	if (q3.back().equals(q4.back()))
		q4.resize(q4.size()-1);
	
	std::reverse(q2.begin(), q2.end());
	std::reverse(q4.begin(), q4.end());
	
	if (q2.back().equals(q3[0]))
		q2.resize(q2.size()-1);
	if(!q4.empty() && q4.back().equals(q1[0]))
		q4.resize(q4.size()-1);
	
	while (y++ < b) {
		q1.push_back(Point(xm, int(ym+y)));
		q3.push_back(Point(xm, int(ym-y)));
	}
	
	res.insert(res.end(), q1.begin(), q1.end());
	res.insert(res.end(), q2.begin(), q2.end());
	res.insert(res.end(), q3.begin(), q3.end());
	res.insert(res.end(), q4.begin(), q4.end());

//	std::cout << "Result" << std::endl;
//	for (Point& p : res)
//		std::cout << p.x << "," << p.y << std::endl;

	return res;
}

std::vector<HSeg> ellipse_fill(int xm, int ym, int a, int b)
{
	std::vector<HSeg> res;
	std::vector<Point> pts;
	
	long long x = -a;
	long long y = 0;
	long long e2 = b;
	long long dx = (1+2*x)*e2*e2;
	long long dy = x*x;
	long long err = dx+dy;
	
	do {
		pts.push_back(Point(int(xm-x), int(ym+y)));
		pts.push_back(Point(int(xm+x), int(ym+y)));
		pts.push_back(Point(int(xm+x), int(ym-y)));
		pts.push_back(Point(int(xm-x), int(ym-y)));
		e2 = 2*err;
		if (e2 >= dx) {
			x++;
			err += dx += 2*((long long)b)*b;
		}
		if (e2 <= dy) {
			y++;
			err += dy += 2*((long long)a)*a;
		}
	} while (x <= 0);
	
	while (y++ < b) {
		pts.push_back(Point(xm, int(ym+y)));
		pts.push_back(Point(xm, int(ym-y)));
	}

	std::sort(pts.begin(), pts.end(), [](Point &a, Point &b) {
		return a.less_than(b);   
	});
	
	if (!pts.empty()) {
		Point& p0 = pts[0];
		HSeg hs(p0.y, p0.x, p0.x);
		for (Point& p : pts) {
			if (p.y != hs.y) {
				res.push_back(hs);
				hs = HSeg(p.y, p.x, p.x);
			} else {
				if (p.x < hs.xl) hs.xl = p.x;
				if (p.x > hs.xr) hs.xr = p.x;
			}
		}
		res.push_back(hs);
	}
	
	return res;
}

//---

void reverse_path(std::vector<Point> & path)
{
	int iLow=0, iHigh=int(path.size()-1);
	while (iLow < iHigh) {
		Point tmp = path[iLow];
		path[iLow] = path[iHigh];
		path[iHigh] = tmp;
		++iLow;
		--iHigh;
	} 
}

void remove_lead_from_path(std::vector<Point> & path, size_t k)
{
	if (k == 0) return;
	size_t i = 0;
	for (; k<path.size(); k++) {
		path[i] = path[k];
		++i;
	}
	path.resize(i);
}

double path_length(std::vector<Point> & path)
{
	double plen = 0.;
	for(size_t j=1; j<path.size(); j++) {
		plen += path[j-1].dist(path[j]);
	}
	return plen;
}

std::vector<Point> path_from_flat(const std::vector<int>& flat_cont, int iStart, int len)
{
	std::vector<Point> path;
	for (int j=0; j<len; j++) {
		path.push_back(Point(flat_cont[iStart], flat_cont[iStart+1]));
		iStart += 2;
	}
	return path;
}

Boundary boundary_around(int x, int y, int dist, int w, int h)
{
	int sz = dist + dist;
	int xtop = w - 1;
	int xmin = x - dist;
	int xmax = x + dist;
	if (xmin < 0) {
		xmin = 0;
		xmax = sz;
		if (xmax > xtop) xmax = xtop;
	}
	if (xmax > xtop) {
		xmax = xtop;
		xmin = xmax - sz;
		if (xmin < 0) xmin = 0;
	}
	int ytop = h - 1;
	int ymin = y - dist;
	int ymax = y + dist;
	if (ymin < 0) {
		ymin = 0;
		ymax = sz;
		if (ymax > ytop) ymax = ytop;
	}
	if (ymax > ytop) {
		ymax = ytop;
		ymin = ymax - sz;
		if (ymin < 0) ymin = 0;
	}
	return Boundary(xmin, ymin, xmax, ymax);
}

Boundary fill_boundary(std::vector<HSeg> & fill)
{
	Boundary b(0,0,0,0);
	if (fill.size() == 0) return b;
	HSeg hs0 = fill[0];
	b.ymin = b.ymax = hs0.y;
	b.xmin = hs0.xl;
	b.xmax = hs0.xr;
	for (HSeg &s : fill) {
		if (b.ymin > s.y) b.ymin = s.y;
		if (b.ymax < s.y) b.ymax = s.y;
		if (b.xmin > s.xl) b.xmin = s.xl;
		if (b.xmax < s.xr) b.xmax = s.xr;
	}
	return b;
}

double fill_centroid(std::vector<HSeg> & fill, double *px, double *py)
{
	double xc = 0., yc = 0.;
	int npt = 0;
	for (HSeg &hs : fill) {
		int y = hs.y;
		for (int x=hs.xl; x<=hs.xr; x++) {
			xc += x;
			yc += y;
			++npt;
		}
	}
	if (npt == 0) {
		*px = *py = 0.;
		return 0.;
	}
	xc /= npt;
	yc /= npt;
	*px = xc;
	*py = yc;

	double r = 0.;
	for (HSeg &hs : fill) {
		double dy2 = yc - hs.y;
		dy2 *= dy2;
		for (int x=hs.xl; x<=hs.xr; x++) {
			double dx = xc - x;
			r += (dy2 + dx*dx);
		}
	}
	return sqrt(2.*r/npt);
}

double fill_circularity(std::vector<HSeg> & fill)
{
	double xc, yc;
	double rad = fill_centroid(fill, &xc, &yc);
	if (rad < 0.0000001) return 0.;
	double radsq = rad * rad;
	int inside = 0;
	for (HSeg &hs : fill) {
		double dy2 = yc - hs.y;
		dy2 *= dy2;
		for (int x=hs.xl; x<=hs.xr; x++) {
			double dx = xc - x;
			if (dx*dx + dy2 <= radsq)
				++inside;
		}
	}

	return double(inside) / (M_PI * radsq);
}

long long fill_overlay_area(std::vector<HSeg> &fill, std::vector<HSeg> &other_fill)
{
	long long a = 0;
	for (HSeg &hs0 : fill) {
		for (HSeg &hs1 : other_fill) {
			if (hs0.y != hs1.y) continue;
			int xmin = hs0.xl < hs1.xl ? hs1.xl : hs0.xl;
			int xmax = hs0.xr > hs1.xr ? hs1.xr : hs0.xr;
			if (xmax < xmin) continue;
			a += (xmax-xmin+1);
		}
	}
	return a;
}

void fill_from_contour(std::vector<HSeg> &fill, Contour &cont)
{
	Boundary bnd = cont.bnd;
	HSeg hs;
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		bool is_in = false;
		hs.y = y;
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (cont.IsInside(Point(x, y))) {
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

void write_particle_data(std::vector<Slice> particles, const char *outfn)
{
	std::filebuf fb;
	fb.open(outfn, std::ios::out);
	if (fb.is_open())
	{
		std::ostream fos(&fb);
		fos << "ID,y,xL,xR\n";
		int id = 1;
		
		for (Slice & ptc : particles) {
			if (ptc.fill.empty()) continue;
			for (HSeg &hs : ptc.fill) {
				fos << id << "," << hs.y << "," << hs.xl << "," << hs.xr << "\n";
			}
			++id;
		}

		fb.close();
	}
}

void write_cell_data(std::vector<Particle3D> &cells, const char *outfn)
{
	std::filebuf fb;
	fb.open(outfn, std::ios::out);
	if (fb.is_open())
	{
		std::ostream fos(&fb);
		fos << "ID,Frame,y,xL,xR\n";
		int id = 1;
		
		for (Particle3D &cell : cells) {
			if (cell.empty()) continue;
			for (int z=0; size_t(z)<cell.fills.size(); z++) {
				for (HSeg &hs : cell.fills[z]) {
					fos << id << "," << z << "," << hs.y << "," << hs.xl << "," << hs.xr << "\n";
				}
			}
			++id;
		}

		fb.close();
	}
}

void write_cell_data(std::vector<Cell> &cells, const char *outfn)
{
	std::filebuf fb;
	fb.open(outfn, std::ios::out);
	if (fb.is_open())
	{
		std::ostream fos(&fb);
		fos << "ID,Frame,y,xL,xR\n";
		int id = 1;
		
		for (Cell &cell : cells) {
			if (cell.empty()) continue;
			for (int z=0; size_t(z)<cell.fills.size(); z++) {
				// if (cell.fills[z].empty()) continue;
				for (HSeg &hs : cell.fills[z]) {
					fos << id << "," << z << "," << hs.y << "," << hs.xl << "," << hs.xr << "\n";
				}
			}
			++id;
		}

		fb.close();
	}
}

int read_cell_data(const char *infn, std::vector<Particle3D> &cells, int w, int h, int d)
{
	CsvCellDataReader rdr(infn);
	if (rdr.eof || !rdr.is_3d) return -1;

// std::cout << "Reading: " << infn << std::endl;

	size_t j0 = cells.size();
	
	int cur_id = -1;
	Particle3D *pcell = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (z<0 || z>=d || y<0 || y>=h || xl<0 || xr>=w || xl>xr) continue;
		if (id != cur_id) {
			cells.resize(cells.size()+1);
			pcell = & cells[cells.size()-1];
			pcell->fills.resize(d);
			cur_id = id;
		}
		if (!pcell) continue;
		pcell->fills[z].push_back(HSeg(y, xl, xr));
	}
	
	for (size_t j=j0; j<cells.size(); j++) {
		Particle3D& cell = cells[j];
		for (int z=0; z<d; z++) cell.fills[z].shrink_to_fit();
		cell.fills.shrink_to_fit();
		cell.update_from_fill();
	}
	return cur_id;
}

int read_particle_data(const char *infn, std::vector<Particle> &particles, int w, int h)
{
	CsvCellDataReader rdr(infn);
	if (rdr.eof || !rdr.is_2d) return -1;

// std::cout << "Reading: " << infn << std::endl;

	size_t j0 = particles.size();
	
	int cur_id = -1;
	Particle *pptc = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (y<0 || y>=h || xl<0 || xr>=w || xl>xr) continue;
		if (id != cur_id) {
			particles.resize(particles.size()+1);
			pptc = & particles[particles.size()-1];
			cur_id = id;
		}
		if (!pptc) continue;
		pptc->fill.push_back(HSeg(y, xl, xr));
	}
	
	for (size_t j=j0; j<particles.size(); j++) {
		Particle& ptc = particles[j];
		ptc.fill.shrink_to_fit();
		ptc.update_from_fill();
	}
	
	return cur_id;
}

//--------------------- Global data --------------------------------

NbrPoint hood_pts[37] = {
	{0, 0},

	{0, -1},
	{1, 0},
	{0, 1},
	{-1, 0},

	{-1, -1},
	{1, -1},
	{1, 1},
	{-1, 1},

	{-1, -2},
	{0, -2},
	{1, -2},
	{2, -1},
	{2, 0},
	{2, 1},
	{1, 2},
	{0, 2},
	{-1, 2},
	{-2, 1},
	{-2, 0},
	{-2, -1},

	{-2, -2},
	{2, -2},
	{2, 2},
	{-2, 2},
	{-1, -3},
	{0, -3},
	{1, -3},
	{3, -1},
	{3, 0},
	{3, 1},
	{1, 3},
	{0, 3},
	{-1, 3},
	{-3, 1},
	{-3, 0},
	{-3, -1}
};

NbrPoint3D hood3d_pts[27] = {
	{0, 0, 0},

	{1, 0, 0},
	{0, -1, 0},
	{-1, 0, 0},
	{0, 1, 0},
	{0, 0, 1},
	{0, 0, -1},

	{-1, -1, 0},
	{1, -1, 0},
	{1, 1, 0},
	{-1, 1, 0},
	{1, 0, 1},
	{0, -1, 1},
	{-1, 0, 1},
	{0, 1, 1},
	{1, 0, -1},
	{0, -1, -1},
	{-1, 0, -1},
	{0, 1, -1},

	{-1, -1, 1},
	{1, -1, 1},
	{1, 1, 1},
	{-1, 1, 1},
	{-1, -1, -1},
	{1, -1, -1},
	{1, 1, -1},
	{-1, 1, -1}
};

//--------------------- Neighborhood Utilities --------------------------------

std::vector<NbrPoint3D> clipped_zhood(int z, int d, int j0, int nbsz)
{
	std::vector<NbrPoint3D> res;
	for (int j=j0; j<nbsz; j++) {
		NbrPoint3D pt = hood3d_pts[j];
		int zz = z + pt.dz;
		if (zz >= 0 && zz < d)
			res.push_back(pt);
	}
	return res;
}
