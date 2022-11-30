#ifndef laminb_h
#define laminb_h

struct LaminB_Result
{
	std::vector<double> dnarad;
	std::vector<double> lmnrad;
	std::vector<double> volume;
};

LaminB_Result lamin_b_analysis(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *csvfile, double thresh, double xyscale=1., double zscale=1.);

// void border_distances(unsigned char *mask, int hm, int wm);

// void calc_dist_3d(double xyscale, double zscale, double maxr);

struct LaminB_Result_2
{
	std::vector<double> dnarad;
	std::vector<double> lmndev;
	std::vector<double> volume;
};

LaminB_Result_2 lamin_b_analysis_2(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *csvfile, double thresh, double xyscale=1., double zscale=1.);

#endif
