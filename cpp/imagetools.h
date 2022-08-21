#ifndef imagetools_h
#define imagetools_h

// Return module version as str
const char* version();

// imagetools.detect_particles_csv(mask, csvfile)
//
// mask is a numpy array of shape (height,width) and type uin8, 0 = background, 255 (0xFF) = foreground
// csvfile (str) is a file name where to write the results (particle data)
int detect_particles_csv(unsigned char *mask, int hm, int wm, const char *csvfile);

// imagetools.detect_particles(mask)
// *** same as imagetools.detect_particles_csv(), but return particle data to Python rather than write CSV ***
//
// mask is a numpy array of shape (height,width) and type uin8, 0 = background, 255 (0xFF) = foreground
// Returns:
// particles - a tuple of tuples of ints containing particle data:
// 		( (y1,xL1,xR1, y2,xL2,xR2, ...), (y1,xL1,xR1, y2,xL2,xR2, ...), ... )
//		Inner touple's length is a multiple of 3, each group of 3 ints is a segment in the form y,xL,xR.
//		Both xL and xR are included, segment length is xR-xL+1; a 1-pixel segment will have xL==xR.
std::vector<std::vector<int>> detect_particles(unsigned char *mask, int hm, int wm);

// imagetools.postprocess_particle_borders(mask)
//
// mask - numpy array of shape=(height,width), dtype=uint8, 0-background, 255-foreground
// Postprocess a mask containing detected cell borders, such as from REShAPE:
//		fill small gaps in the borders using dilations/erosions, and thin the borders down to 1 pixel.
void postprocess_particle_borders(unsigned char *mask, int hm, int wm);

// particles = imagetools.masks_to_particles(masks, [x_orig, y_orig, [minarea]])
//
// masks - numpy array, shape=(num_masks, height, width) dtype=uint8, background=0, foreground=255
//		(as returned by Mask_RCNN);
// Optional parameters:
// x_orig, y_orig - tile coordinate origins. Output coordinates are translated by adding the origin -
//		good for reconstructing tiled images (just concatenate the outputs).
// minarea - if detected mask consists of several disconnected areas, filter out those smaller than this;
// Returns:
// particles - a tuple of tuples of ints containing particle data:
// 		( (y1,xL1,xR1, y2,xL2,xR2, ...), (y1,xL1,xR1, y2,xL2,xR2, ...), ... )
//		Inner touple's length is a multiple of 3, each group of 3 ints is a segment in the form y,xL,xR.
//		Both xL and xR are included, segment length is xR-xL+1; a 1-pixel segment will have xL==xR.
std::vector<std::vector<int>> masks_to_particles(unsigned char *ptmask, int npts, int hptm, int wptm,
		int x_orig=0, int y_orig=0, int minarea=100);

// Constants for postproc parameter in assemble_ml()
const int POSTPROC_NONE = 0;
const int POSTPROC_DNA = 1;
const int POSTPROC_ACTIN = 2;

void assemble_ml(std::vector<std::vector<std::vector<int>>> particles_3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *csvfile, int postproc=POSTPROC_NONE,
		double good_iou=0.6, double ok_iou=0.5);

struct Compare3dResult {
	int base_cells;
	int cmp_cells;
	int base_slices;
	int cmp_slices;

	int f_pos, f_neg, fragm, fused;
	std::vector<int> pct_match;

	int f_pos_3d, f_neg_3d, fragm_3d, fused_3d;
	std::vector<int> pct_match_3d;
};

Compare3dResult compare_3d_annotations(int w, int h, int d, const char *base_csv, const char *cmp_csv);

#endif
