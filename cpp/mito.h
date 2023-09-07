#ifndef mito_h
#define mito_h

void extract_z01(unsigned short *data3d, int zd3d, int hd3d, int wd3d);

struct MitoResult {
	long long total_count;
	long long mito_count;
	std::vector<long long> counts_by_dist;
	double median_dist;
	std::vector<long long> spread_by_z;
	std::vector<long long> total_by_z;
	double z_mean;
	double z_stdev;
};

MitoResult detect_mito(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		unsigned char *mask, int hm, int wm,
		const char *csvfile);

// --- New ---
void segment_mito(unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *csvfile);

void analyze_mito(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *mitocsv, const char *dnacsv, const char *actincsv, const char *datacsv,
		double scx=1., double scy=1., double scz=1.);

void detect_dsp(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile,
		int min_volume = 21, double min_snr = 3., double sensitivity = 3.);

void detect_dsp_upd(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d);

struct DSPResult
{
	long long pix_count;
	std::vector<long long> pix_dist_inside;
	std::vector<long long> pix_dist_outside;
};

DSPResult analyze_dsp(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *dspcsv, const char *actincsv,
		const char *outdspcsv, const char *outcellcsv);

// colorize_dsp(src_scaled (h*w), colorized(3*5h*5w), z, actincsv, indspcsv)
void colorize_dsp_x5(unsigned char *mask, int hm, int wm,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		int z, const char *actincsv, const char *indspcsv);

#endif
