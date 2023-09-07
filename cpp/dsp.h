#ifndef dsp_h
#define dsp_h

#include "raster.h"

struct DSP_Seed
{
	std::vector<Point3D> pts;
	double score;
	DSP_Seed() : score(0.) {}
};

struct DSP_Score
{
	size_t idx;
	double score;
	DSP_Score() : idx(0), score(0.) {}
	DSP_Score(size_t _idx, double _score) : idx(_idx), score(_score) {}
};

class DSP_Detector
{
public:
	int w, h, d;
	Raster16_3D srcstack;
	Raster3D mstack;
	
	//--- Adjustable parameters ---
	double sensitivity = 3.;		// 1..5; higher values = bigger detected objects
	long long min_volume = 21;		// 0..200; filter out anything smaller than min_volume voxels
	double min_snr = 3.;			// 1.5..5; filter out particles having local SNR lower than this
	
	DSP_Detector(unsigned short *data3d, int zd3d, int hd3d, int wd3d, unsigned char *mask3d) :
		w(wd3d), h(hd3d), d(zd3d),
		srcstack(wd3d, hd3d, zd3d, data3d),
		mstack(wd3d, hd3d, zd3d, mask3d) {
	}
	virtual ~DSP_Detector() {}
	std::vector<Particle3D> detect_particles();
	
	void detect_upd();
	
protected:
	Histogram normalize_data(Raster16_3D& dstack);
	double particle_snr(Particle3D& cell);
};

struct CellBorder
{
	int idx;
	std::vector<std::vector<Point>> bpixels;
	double closest(double xc, double yc, double zc) {
		double bsf = 1000000000.;
		for (int z=0; size_t(z)<bpixels.size(); z++) {
			double dz = zc - z;
			for (Point& p : bpixels[z]) {
				double dx = xc - p.x;
				double dy = yc - p.y;
				double rsq = dx*dx + dy*dy + dz*dz;
				if (rsq < bsf) bsf = rsq;
			}
		}
		return sqrt(bsf);
	}
};

struct DSPtoCell
{
	int dspidx;
	int cellidx;
	int zmin, zmax;
	double xc, yc, zc;
	double value;
	long long volume;
	long long overlay;
	double border_dist;
};

class DSP_Analyzer
{
public:
	int w, h, d;
	Raster16_3D srcstack;
	Raster3D mstack;
	
	std::vector<Particle3D> dsps;
	std::vector<Particle3D> cells;
	std::vector<CellBorder> borders;
	std::vector<Point3D> centroids;

	DSP_Analyzer(unsigned short *data3d, int zd3d, int hd3d, int wd3d, const char *dspcsv, const char *actincsv) :
		w(wd3d), h(hd3d), d(zd3d),
		srcstack(wd3d, hd3d, zd3d, data3d),
		mstack(wd3d, hd3d, zd3d, NULL)
	{
		read_cell_data(dspcsv, dsps, w, h, d);
		read_cell_data(actincsv, cells, w, h, d);
		detect_cell_borders();
	}
	virtual ~DSP_Analyzer() {}
	long long pix_distances(std::vector<long long>& pix_dist_inside,
			std::vector<long long>& pix_dist_outside, int maxdist=20);
	std::vector<DSPtoCell> dsp_to_cell_match(int maxdist=20, double close_by=0.);
	void paint_dsp_mask(std::vector<DSPtoCell>& dsptocells);

protected:
	void detect_cell_borders();
	DSPtoCell find_closest_cell(int dspidx, int maxdist);
};

#endif
