#ifndef fstack_h
#define fstack_h

#include "raster.h"

struct NbrDist3D
{
	int dx, dy, dz;
	double dist;
	NbrDist3D() {}
	NbrDist3D(int _dx, int _dy, int _dz, double _dist) : dx(_dx), dy(_dy), dz(_dz), dist(_dist) {}
	
	static std::vector<NbrDist3D> list_dist_3d(double xsc, double ysc, double zsc, double maxr);
};

bool comp_nbrdist3d(const NbrDist3D& a, const NbrDist3D& b);

template <class T>
class RasterF3D {
public:
	int w, h, d;
	T* buf;
	long long len;
	long long plane_len;
	bool is_temp;
	double scx, scy, scz;
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
		scx = scy = scz = 1.;
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

void distance_map(RasterF3D<float>& dstack, std::vector<Particle3D>& particles, std::vector<Particle3D>& cells);

#endif

