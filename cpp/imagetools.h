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
// expand - (optional) expand borders by 1 pixel if true (1-pix lines -> 3-pix lines).
void postprocess_particle_borders(unsigned char *mask, int hm, int wm, bool expand=false);

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
const int POSTPROC_SANDPAPER = 0x10;

// postproc - postprocessing flag.
//		The POSTPROC_SANDPAPER value can be used alone or in combination with DNA/ACTIN (using '|'):
// 			imagetools.assemble_ml(particles_3d, mask3d, csvfile, POSTPROC_SANDPAPER)
// 			imagetools.assemble_ml(particles_3d, mask3d, csvfile, POSTPROC_SANDPAPER|POSTPROC_DNA)
// 			imagetools.assemble_ml(particles_3d, mask3d, csvfile, POSTPROC_SANDPAPER|POSTPROC_ACTIN)
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

// result = imagetools.border_pixels(w, h, d, csvfile)
//
// w,h,d (int) - stack dimensions (width, height, number of frames)
// csvfile (str) - path to .csv containing particle data (ID,Frame,y,xL,xR)
//
// result - tuple of tuples of tuples, formatted like this:
//	(
//		((x0,y0,z0), (x1,y1,z1), ...),	# object 1
//		((x0,y0,z0), (x1,y1,z1), ...),	# object 2
//		((x0,y0,z0), (x1,y1,z1), ...),	# object 3
//		...
//	)
//
// * (xN,yN,zN) -- coordinates of a border pixel
std::vector<std::vector<std::vector<int>>> border_pixels(int w, int h, int d, const char *csvfile);

// imagetools.rs_tops_bottoms(mask3d)
//		Complement RESHAPE_3D GT images with cell tops and bottoms, wherever there are detectable cells, i.e.
//		cells recognized by imagetools.assemble_ml().
//
// mask3d - numpy array, shape=(num_frames, height, width), dtype=np.uint8; background=0, foreground=255
//		as in Ground Truth for RERSHAPE_3D (foreground = cell borders).
// xcv (int) - optional (default 127) pixel value for extrapolated cell parts.
// 		On return, mask3d updated with extrapolated cell tops and bottoms painrted with color specified by xcv.
void rs_tops_bottoms(unsigned char *mask3d, int zm3d, int hm3d, int wm3d, int xcv=0x7F);

// imagetools.sandpaper_cells(w, h, d, in_csv, out_csv)
//
// w,h,d (int) - stack dimensions (width, height, number of frames)
// in_csv (str) - path to .csv containing input cell data (ID,Frame,y,xL,xR)
// out_csv (str) - path to .csv where to write output cell data (ID,Frame,y,xL,xR)
void sandpaper_cells(int w, int h, int d, const char *in_csv, const char *out_csv);

// imagetools.paint_cells(mask3d, in_csv)
//		Paint cells into an empty numpy array, RPE Map style (gray=background, black=cell interiors, white=cell borders)
//
// mask3d - empty numpy array, shape=(num_frames, height, width), dtype=np.uint8.
// in_csv (str) - path to .csv containing input cell data (ID,Frame,y,xL,xR)
void paint_cells(unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *in_csv);

#endif
