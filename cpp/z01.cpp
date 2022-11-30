#include "raster.h"
#include "z01.h"

// Used by extract()
static void smooth_z_mask(Raster8 &xmsk)
{
	Raster8 tmsk(xmsk.w, xmsk.h, NULL);
	memcpy(tmsk.buf, xmsk.buf, xmsk.len);
	for (int y=2; y<xmsk.h-2; y++) {
		unsigned char *p = tmsk.scanLine(y);
		for (int x=2; x<xmsk.w-2; x++) {
			int csum = 0;
			for (int j=0; j<HOOD_SIZE_FATCROSS; j++) {
				csum += int(xmsk.value(x+hood_pts[j].dx, y+hood_pts[j].dy));
			}
			p[x] = (unsigned char)(csum / HOOD_SIZE_FATCROSS);
		}
	}
	memcpy(xmsk.buf, tmsk.buf, xmsk.len);
}

// Input: 3D array of 16-bit (n_frames, height, width) containing Z01 channel raw data;
// Output:	Frame 0 replaced with aggregated 2D Z01 data
//			Frame 1 - Lower Z (ranging 0 through n_frames-1)
//			Frame 2 - Upper Z (ranging 0 through n_frames-1)
void extract(unsigned short *data3d, int zd3d, int hd3d, int wd3d)
{
	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);

	// Detect frame range loz..hiz suitable for computing Otsu
	// Try to exclude first 3 frames often containing heavy static
	double *fr_dev = new double[dstack.d];
	double hidev = 0.;
	int loz = 0, hiz = 0;
	for (int z=dstack.d-1; z>=0; z--) {
		Raster16 dat = dstack.getPlane(z);
		double dev;
		dat.mean_std(&dev);
		fr_dev[z] = dev;
		if (hidev < dev && z > 2) {
			hidev = dev;
			hiz = z;
		}
	}
	double middev = hidev * 0.75;
	double upperdev = hidev * 1.02;
	for (loz=hiz-1; loz>=0; loz--) {
		if (fr_dev[loz] < middev || fr_dev[loz] >= upperdev) break;
	}
	++loz;
	for (; hiz<dstack.d; hiz++) {
		if (fr_dev[hiz] < middev) break;
	}
	--hiz;
	delete [] fr_dev;
	
	// std::cout << "loz=" << loz << " hiz=" << hiz << std::endl;
	if (loz < 0) loz = 0;
	if (hiz >= dstack.d) hiz = dstack.d - 1;
	if (loz > hiz) loz=hiz;
	
	Histogram hist(0x1000);
	for (int z=loz; z<=hiz; z++) {
		Raster16 dat = dstack.getPlane(z);
		for (int y=0; y<dat.h; y++) {
			hist.add_row16(dat.scanLine(y), dat.w);
		}
	}
	
	unsigned short otsu = hist.otsu16();
	// Half Otsu seems to work best for any week W1-W6
	otsu /= 2;

	Raster3D mstack(wd3d, hd3d, zd3d, NULL);
	for (long long i=0; i<mstack.len; i++) {
		mstack.buf[i] = (dstack.buf[i] > otsu) ? 0xFF : 0;
	}
	
	Raster8 tmsk(wd3d, hd3d, NULL);
	Raster8 xmsk(wd3d, hd3d, NULL);
	xmsk.fill(0xFF);
	
	// Start from the top, follow the pattern.
	// Exclude anything that conflicts with pattern found 8 (or more) frames above the current frame.
	for (int z=mstack.d-1; z>=0; z--) {
		Raster8 msk = mstack.getPlane(z);
		msk.expandBorders(0xFF, 0x80, HOOD_SIZE_MOORE, 0);
		if (z < mstack.d-1) {
			Raster8 top = mstack.getPlane(z+1);
			for (long long i=0; i<msk.len; i++)
				if (top.buf[i] == 0xFF && msk.buf[i] != 0xFF) msk.buf[i] = 0x80;
		}
		if (z > 0) {
			Raster8 bot = mstack.getPlane(z-1);
			for (long long i=0; i<msk.len; i++)
				if (bot.buf[i] == 0xFF && msk.buf[i] != 0xFF) msk.buf[i] = 0x80;
		}
		msk.expandBorders(0x80, 0x40, HOOD_SIZE_MOORE, 0);
		msk.replaceColor(0x40, 0xFF);
		msk.replaceColor(0x80, 0xFF);
		msk.fillBorder(0x10, 1);
		msk.filterParticles(0xFF, 0x80, 1000, 0);
		
		if (z+6 < mstack.d-1)
		{
			Raster8 top = mstack.getPlane(z+1);
			Raster8 top2 = mstack.getPlane(z+2);
			for (int y0=1; y0<msk.h-1; y0++) {
				unsigned char *p = msk.scanLine(y0);
				for (int x0=1; x0<msk.w-1; x0++) {
					if (p[x0] != 0x80) continue;
					Slice ptc;
					ptc.area = msk.detectParticle(ptc, x0, y0, 0xA0);
					
					int ncgood = 0;
					for (HSeg &hs : ptc.fill) {
						unsigned char *p1 = top.scanLine(hs.y);
						unsigned char *p2 = top2.scanLine(hs.y);
						for (int x=hs.xl; x<=hs.xr; x++) {
							if (p1[x] >= 0x80 || p2[x] >= 0x80)
								++ncgood;
						}
					}
					
					int ncbad = ptc.area - xmsk.countColor(ptc, 0xFF);
					
					if (ncbad > int(ptc.area * 0.4) ||
							(ncbad > ncgood && ncbad > int(ptc.area * 0.05))) {
						msk.paintParticle(ptc, 0x40);
					}
				}
			}
			msk.replaceColor(0xA0, 0x80);
			msk.replaceColor(0x40, 0);
		}
		
		
		for (int pass=0; pass<3; pass++) {
			msk.expandBorders(0x80, 0x40, HOOD_SIZE_RAD3, 0);
			msk.replaceColor(0x40, 0x80);
		}
		msk.filterParticles(0, 0x20, 10000, 0x40);
		msk.replaceColor(0x40, 0x80);
		msk.replaceColor(0x20, 0);
		
		if (z+8 < mstack.d-1) {
			Raster8 top = mstack.getPlane(z+8);
			for (long long i=0; i<xmsk.len; i++)
				if (top.buf[i] >= 0x80 && xmsk.buf[i] == 0xFF)
					xmsk.buf[i] = (unsigned char)(z);
		}
		
		Raster16 dat = dstack.getPlane(z);
		for (long long i=0; i<msk.len; i++) {
			if (msk.buf[i] == 0x80 && dat.buf[i] > otsu)
				msk.buf[i] = 0xFF;
		}
	}
	
	// Find upper Z for every pixel in the detected pattern
	xmsk.fill(0xFF);
	for (int z=mstack.d-1; z>=0; z--) {
		int zz = z + 1;
		if (zz >= mstack.d) zz = mstack.d - 1;
		else if (zz < 3) zz = 3;
		Raster8 msk = mstack.getPlane(z);
		for (long long i=0; i<xmsk.len; i++)
			if (msk.buf[i] >= 0x80 && xmsk.buf[i] == 0xFF)
				xmsk.buf[i] = (unsigned char)(zz);
	}
	
	// Try to fill the gaps (averaging Z-values at the perimeter)
	xmsk.fillBorder(0xFE, 1);
	for (int y0=1; y0<xmsk.h-1; y0++) {
		unsigned char *p = xmsk.scanLine(y0);
		for (int x0=1; x0<xmsk.w-1; x0++) {
			if (p[x0] != 0xFF) continue;
			Slice ptc;
			ptc.area = xmsk.detectParticle(ptc, x0, y0, 0xF0);
			
			Boundary bnd = ptc.bnd;
			bnd.expand(3);
			xmsk.clip(bnd);
			bnd.expand(-2);
			
			int npts = 1;
			int csum = dstack.d - 1;

			for (int y=bnd.ymin; y<=bnd.ymax; y++) {
				for (int x=bnd.xmin; x<=bnd.xmax; x++) {
					unsigned char c = xmsk.value(x,y);
					if (c >= 0xF0) continue;
					for (int j=1; j<HOOD_SIZE_FATCROSS; j++) {
						if (xmsk.value(x+hood_pts[j].dx, y+hood_pts[j].dy) >= 0xF0) {
							++npts;
							csum += int(c);
							break;
						}
					}
				}
			}
			if (npts > 0) {
				csum = (csum / npts) + 3;
				if (csum >= dstack.d) csum = dstack.d - 1;
			}

			xmsk.paintParticle(ptc, (unsigned char)csum);
		}
	}
	
	// Fill border values by duplicating the values next to border
	for (int y0=1; y0<xmsk.h-1; y0++) {
		unsigned char *p = xmsk.scanLine(y0);
		p[0] = p[1];
		int x = xmsk.w-1;
		p[x] = p[x-1];
	}
	memcpy(xmsk.scanLine(0), xmsk.scanLine(1), xmsk.w * sizeof(unsigned char));
	memcpy(xmsk.scanLine(xmsk.h-1), xmsk.scanLine(xmsk.h-2), xmsk.w * sizeof(unsigned char));
	
	// Smooth it a little
	smooth_z_mask(xmsk);
	
	// Now find lower Z 
	tmsk.fill(0xFF);
	for (int z=0; z<mstack.d-1; z++) {
		int zz = z - 1;
		if (zz < 0) zz = 0;
		else if (zz+3 >= mstack.d-1) zz = mstack.d-4;
		Raster8 msk = mstack.getPlane(z);
		for (long long i=0; i<tmsk.len; i++)
			if (msk.buf[i] >= 0x80 && tmsk.buf[i] == 0xFF)
				tmsk.buf[i] = (unsigned char)(zz);
	}
	// Fill gaps with values UpperZ - 8
	for (long long i=0; i<xmsk.len; i++) {
		if (tmsk.buf[i] == 0xFF) {
			int c = int(xmsk.buf[i]) - 8;
			if (c < 0) c = 0;
			tmsk.buf[i] = (unsigned char)(c);
		}
	}
	
	// Smooth lower Z
	smooth_z_mask(tmsk);

	// Aggregate pixel values by taking an average of 3 brightest pixels in the range LowerZ..UpperZ
	// Store the result in frame 0
	for (int y=0; y<dstack.h; y++) {
		unsigned short* tgt = dstack.scanLine(y, 0);
		unsigned char *plo = tmsk.scanLine(y);
		unsigned char *phi = xmsk.scanLine(y);
		for (int x=0; x<dstack.w; x++) {
			int z0 = int(plo[x]);
			int z1 = int(phi[x]);
			if (z1-z0 < 4) {
				z1 = z0 + 4;
				if (z1 >= dstack.d) z1 = dstack.d - 1;
			}
			if (z1-z0 < 2) {
				z0 = z1 - 2;
				if (z0 < 0) z0 = 0;
			}
			unsigned short v1=0, v2=0, v3=0;
			for (int z=z0; z<=z1; z++) {
				unsigned short v = dstack.value(x, y, z);
				if (v > v1) {
					v3 = v2;
					v2 = v1;
					v1 = v;
				} else if (v > v2) {
					v3 = v2;
					v2 = v;
				} else if (v > v3) {
					v3 = v;
				}
			}
			tgt[x] = (unsigned short)((int(v1) + int(v2) + int(v3)) / 3);
		}
	}

	// Store LowerZ values in frame 1, and UpperZ, in frame 2
	unsigned short *fr1 = dstack.scanLine(0, 1);
	unsigned short *fr2 = dstack.scanLine(0, 2);
	for (long long i=0; i<xmsk.len; i++) {
		fr1[i] = (unsigned short)(tmsk.buf[i]); // << 3;
		fr2[i] = (unsigned short)(xmsk.buf[i]); // << 3;
	}
}
