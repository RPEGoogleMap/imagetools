#include "raster.h"
#include "fstack.h"
#include "csv.h"
#include "dsp.h"

#include "mito.h"

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
void extract_z01(unsigned short *data3d, int zd3d, int hd3d, int wd3d)
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


struct ZParticle : public Particle
{
	int area;
	int zmin, zmax;
	ZParticle() : Particle(), area(0), zmin(-1), zmax(-1) {}
};

MitoResult detect_mito(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		unsigned char *mask, int hm, int wm,
		const char *csvfile)
{
	Raster3D mstack(wm3d, hm3d, zm3d, mask3d);
	Raster8 zmsk(wm, hm, mask);
	
	std::vector<Particle3D> nuclei;
	read_cell_data(csvfile, nuclei, mstack.w, mstack.h, mstack.d);
	std::cout << "Read " << nuclei.size() << " nuclei from " << csvfile << std::endl;

	int z0=-1, z1=-1;
	std::vector<HSeg> fill;
	
	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.expandBorders(0, 0x40, HOOD_SIZE_MOORE, 0x80);
		msk.replaceColor(0x40, 0);
		
		msk.expandBorders(0xFF, 0xE0, HOOD_SIZE_RAD3, 0x80);
		msk.replaceColor(0xE0, 0xFF);
		//
		for (int y=0; y<mstack.h; y++) {
			unsigned char *p = msk.scanLine(y);
			unsigned char *zp = zmsk.scanLine(y);
			for (int x=0; x<mstack.w; x++) {
				if (zp[x] != 0 && p[x] == 0x80) p[x] = 0;
			}
		}
		//
		msk.expandBorders(0, 0x40, HOOD_SIZE_MOORE, 0x80);
		//
		msk.expandBorders(0xFF, 0xE0, HOOD_SIZE_RAD3, 0x80);
		//
		
		bool found = false;
		
		msk.fillBorder(0, 1);
		for (int y=1; y<mstack.h-1; y++) {
			unsigned char *p = msk.scanLine(y);
			for (int x=1; x<mstack.w-1; x++) {
				if (p[x] != 0x80) continue;
				found = true;
				
				for (int j=1; j<HOOD_SIZE_NEUMANN; j++) {
					if (msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy) == 0xE0) {
						fill.clear();
						msk.findParticleFill(fill, x, y, 0xD0);
						msk.paintParticleFill(fill, 0xE0);
						break;
					}
				}
			}
		}
		msk.replaceColor(0xE0, 0xFF);
		if (found) {
			//msk.replaceColor(0x80, 0x40);
			if (z0 < 0) z0 = z;
			z1 = z;
		}
	}
	std::cout << "z0=" << z0 << " z1=" << z1 << std::endl;
	
	for (int z=z0; z<=z1; z++) {
		std::vector<NbrPoint3D> nbpts = clipped_zhood(z, mstack.d, 1, HOOD3D_6);
		
		for (int y=1; y<mstack.h-1; y++) {
			unsigned char *p = mstack.scanLine(y, z);
			for (int x=1; x<mstack.w-1; x++) {
				if (p[x] != 0x80) continue;
				
				for (NbrPoint3D& pt : nbpts) {
					if (mstack.value(x+pt.dx, y+pt.dy, z+pt.dz) == 0xFF) {
						p[x] = 0xE0;
						break;
					}
				}
			}
		}
	}
	
	mstack.replaceColor(0xE0, 0xFF);
	mstack.replaceColorInv(0xFF, 0);
	
	for (Particle3D& nuc : nuclei) {
		mstack.paintParticleInto(nuc, 0x40, 0xFF);
		mstack.paintParticleInto(nuc, 0x20, 0);
	}
	
	double zarea = 0.;
	std::vector<ZParticle> zcells;
	zmsk.fillBorder(0x10, 1);
	for (int y=1; y<zmsk.h-1; y++) {
		unsigned char *p = zmsk.scanLine(y);
		for (int x=1; x<zmsk.w-1; x++) {
			if (p[x] != 0) continue;
			
			zcells.resize(zcells.size()+1);
			ZParticle& zpt = zcells[zcells.size()-1];
			zpt.area = zmsk.detectParticle(zpt, x, y, 0x50);
			zarea += double(zpt.area);
			zmsk.paintParticle(zpt, 0x20);
		}
	}
	
	int minar = 0;
	int valid = 0;
	if (zcells.size() > 0) {
		zarea /= double(zcells.size());
		minar = int(zarea / 3.);
		std::cout << "zarea=" << zarea << " minar=" << minar << std::endl;
		for (ZParticle& zpt : zcells) {
			if (zpt.area < minar) {
				zmsk.paintParticle(zpt, 0xFF);
				zpt.clear();
			} else {
				++valid;
			}
		}
	}
	
	std::cout << "Detected " << valid << " (" << zcells.size() << ") cells on Z01 mask." << std::endl;
	
	z0 = z1 = -1;
	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		for (ZParticle& zpt : zcells) {
			if (zpt.fill.empty()) continue;
			if (msk.countColor(zpt, 0) != zpt.area) {
				if (zpt.zmin < 0) zpt.zmin = z;
				zpt.zmax = z;
				if (z0 < 0) z0 = z;
				z1 = z;
			}
		}
	}
	std::cout << "Now z0=" << z0 << " z1=" << z1 << std::endl;
	
	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		for (ZParticle& zpt : zcells) {
			if (zpt.fill.empty()) continue;
			if (z < zpt.zmin || z > zpt.zmax)
				msk.paintParticle(zpt, 0x40);
		}
		//for (long long i=0; i<msk.len; i++) {
		//	if (zmsk.buf[i] == 0xFF) msk.buf[i] = 0x80;
		//}
	}
	
	if (z0 < 1) z0 = 1;
	if (z1+1 >= mstack.d) z1 = mstack.d - 2;
	
	unsigned char val0 = 0xFF;
	unsigned char val1 = 0xC0;
	unsigned char cval = val0;
	unsigned char nval = val1;
	for (int pass=0; pass<50; pass++) {
		//std::cout << "Pass " << pass << std::endl;
		int nbsz = (pass & 1) ? HOOD3D_18 : HOOD3D_6;
		
		for (int z=z0; z<=z1; z++) {
			for (int y=1; y<mstack.h-1; y++) {
				unsigned char *p = mstack.scanLine(y, z);
				for (int x=1; x<mstack.w-1; x++) {
					if (p[x] != 0) continue;
					
					for (int j=1; j<nbsz; j++) {
						if (mstack.value(x+hood3d_pts[j].dx, y+hood3d_pts[j].dy, z+hood3d_pts[j].dz) == cval) {
							p[x] = nval;
							break;
						}
					}
				}
			}
		}
		
		cval = nval;
		nval -= 3;
	}
	
	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		for (long long i=0; i<msk.len; i++) {
			if (zmsk.buf[i] == 0xFF) msk.buf[i] = 0x40;
		}
	}
	mstack.replaceColor(0, nval);
	
	std::cout << "Constructing result" << std::endl;
	
	MitoResult res;
	
	long long tot = 0;
	long long pos = 0;
	std::vector<long long>& cnt_by_dist = res.counts_by_dist;
	cnt_by_dist.resize(200);
	std::vector<long long>& pos_by_z = res.spread_by_z;
	pos_by_z.resize(200);
	std::vector<long long>& tot_by_z = res.total_by_z;
	tot_by_z.resize(200);
	for (size_t i=0; i<200; i++) {
		cnt_by_dist[i] = 0;
		pos_by_z[i] = 0;
		tot_by_z[i] = 0;
	}
	
	std::cout << "Counting..." << std::endl;
	
	for (ZParticle& zpt : zcells) {
		if (zpt.fill.empty() || zpt.zmin < 0) continue;
		for (int z=zpt.zmin; z<=zpt.zmax; z++) {
			Raster8 msk = mstack.getPlane(z);
			int dz = z - zpt.zmin;
			
			for (HSeg& hs : zpt.fill) {
				unsigned char *p = msk.scanLine(hs.y);
				for (int x=hs.xl; x<=hs.xr; x++) {
					unsigned char c = p[x];
					if (c < nval) continue;
					int d = 0;
					if (c <= val1) {
						d = int(val1 - c) / 3 + 1;
					} else {
						++pos_by_z[dz];
						++pos;
					}
					++tot_by_z[dz];
					++cnt_by_dist[d];
					++tot;
				}
			}
		}
	}
	
	res.total_count = tot;
	res.mito_count = pos;
	
	std::cout << "Pos: " << pos << " Tot: " << tot << std::endl;

	for (int i=1; size_t(i+1)<cnt_by_dist.size(); i+=2) {
		long long dv = (cnt_by_dist[i+1] - cnt_by_dist[i]) / 2;
		cnt_by_dist[i] += dv;
		cnt_by_dist[i+1] -= dv;
	}
	
	for (int i=199; i>0; i--) {
		if (cnt_by_dist[i] != 0) {
			cnt_by_dist.resize(i+1);
			break;
		}
	}
	for (int i=199; i>0; i--) {
		if (tot_by_z[i] != 0) {
			pos_by_z.resize(i+1);
			tot_by_z.resize(i+1);
			break;
		}
	}

	double median = 0.;
	double mean_z = 0.;
	double stdev_z = 0.;
	std::cout << "Pos: " << pos << " Tot: " << tot << std::endl;
	if (pos > 0)
	{
		long long half_neg = (tot - pos) / 2;
		long long acc = 0;
		for (int i=1; size_t(i)<cnt_by_dist.size(); i++) {
			acc += cnt_by_dist[i];
			if (acc >= half_neg) {
				median = double(i) - double(acc - half_neg)/cnt_by_dist[i];
				std::cout << "Z: " << i << " Median: " << median << " acc=" << acc << " half_neg=" << half_neg << std::endl;
				break;
			}
		}
		
		for (int z=0; size_t(z)<pos_by_z.size(); z++) {
			mean_z += double(z) * pos_by_z[z];
		}
		mean_z /= double(pos);
		for (int z=0; size_t(z)<pos_by_z.size(); z++) {
			double dz = double(z) - mean_z;
			stdev_z += dz*dz * pos_by_z[z];
		}
		stdev_z = sqrt(stdev_z / double(pos));
	}
	std::cout << "Mean: " << mean_z << " StdDev: " << stdev_z << std::endl;
	res.median_dist = median;
	res.z_mean = mean_z;
	res.z_stdev = stdev_z;
	
	return res;
}

void segment_mito(unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *csvfile)
{
	Raster3D mstack(wm3d, hm3d, zm3d, mask3d);
	Raster8 msk0 = mstack.getPlane(0);
	
	for (int pass=0; pass<2; pass++)
	{
		int nbsz = pass&1 ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE;
		for (int z=0; z<mstack.d; z++) {
			Raster8 msk = mstack.getPlane(z);
			msk.fillBorder(0, 1);
			msk.expandBorders(0, 0x40, nbsz, 0x80);
			msk.expandBorders(0x80, 0xC0, nbsz, 0xFF);
			msk.replaceColor(0x40, 0);
			msk.replaceColor(0xC0, 0x80);
		}
	}
	
	std::vector<Particle3D> particles;
	//std::vector<long long> volumes;
	for (int z=0; z<mstack.d; z++) {
		for (int y=1; y<mstack.h-1; y++) {
			unsigned char *p = mstack.scanLine(y, z);
			for (int x=1; x<mstack.w-1; x++) {
				if (p[x] != 0xFF) continue;
				
				particles.resize(particles.size() + 1);
				Particle3D& pt = particles[particles.size()-1];
				long long vol = mstack.detectParticle(pt, x, y, z, 0xF0);
				if (vol < 25) {
					mstack.paintParticle(pt, 0x40);
					particles.resize(particles.size() - 1);
				} else {
					mstack.paintParticle(pt, 0xC0);
					//volumes.push_back(vol);
				}
			}
		}
	}
	
	for (Particle3D& pt : particles) {
		Boundary bnd = pt.bnd.boundary2d();
		bnd.expand(1);
		msk0.clip(bnd, 1);
		for (int z=pt.bnd.zmin; z<=pt.bnd.zmax; z++) {
			std::vector<HSeg>& fill = pt.fills[z];
			Raster8 msk = mstack.getPlane(z);
			
			msk.paintParticleFill(fill, 0xFF);
			for (int y=bnd.ymin; y<=bnd.ymax; y++) {
				unsigned char *p = msk.scanLine(y);
				for (int x=bnd.xmin; x<=bnd.xmax; x++) {
					if (p[x] != 0x80) continue;
					bool mine=false;
					bool their=false;
					
					for (int j=1; j<HOOD_SIZE_MOORE; j++) {
						unsigned char c = msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy);
						if (c == 0xFF) mine = true;
						else if (c == 0xC0) their = true;
					}
					if (mine) {
						p[x] = their ? 0 : 0xFF;
					}
				}
			}
			msk.rescanParticleFill(bnd, fill, 0xFF);
			msk.paintParticleFill(fill, 0xC0);
		}
	}
	
	//std::cout << "Detected " << particles.size() << " particles." << std::endl;
	
	mstack.fill(0x80);
	for (Particle3D& pt : particles)
		mstack.paintParticle(pt, 0);
	mstack.expandBorders(0x80, 0, 0xFF, HOOD3D_26);

	std::cout << "Write " << particles.size() << " objects to " << csvfile << std::endl;
	write_cell_data(particles, csvfile);
}

static double mito_dist_to_center(Raster16_3D& dstack, Particle3D& mito,
		std::vector<Particle3D>& cells, std::vector<Point3D>& cell_centers,
		double scx, double scy, double scz)
{
	double res = 1e12;
	Boundary mbnd = mito.bnd.boundary2d();
	mbnd.expand(5);
	
	for (size_t j=0; j<cells.size(); j++) {
		Particle3D& cell = cells[j];
		if (!cell.bnd.intersects2d(mbnd)) continue;
		Point3D& cc = cell_centers[j];
		double dist = 0.;
		double mass = 0.;
		for (int z=0; z<dstack.d; z++) {
			if (cell.fills[z].empty()) continue;
			double dz = (z - cc.z) * scz;
			dz *= dz;
			for (HSeg& hs : cell.fills[z]) {
				double dy = (hs.y - cc.y) * scy;
				dy *= dy;
				unsigned short *buf = dstack.scanLine(hs.y, z);
				for (int x=hs.xl; x<=hs.xr; x++) {
					double dx = (x - cc.x) * scx;
					dx *= dx;
					mass += double(buf[x]);
					dist += (dx + dy + dz) * double(buf[x]);
				}
			}
		}
		dist = sqrt(dist / mass);
		if (dist < res) res = dist;
	}
	
	return res;
}

void analyze_mito(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *mitocsv, const char *dnacsv, const char *actincsv, const char *datacsv,
		double scx, double scy, double scz)
{
	CsvWriter wr(datacsv);
	if (!wr.is_open()) {
		std::cout << "Cannot write to " << datacsv << std::endl;
		return;
	}
	wr.SetDoubleFormat("%0.3lf");

	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);
	
	std::vector<Particle3D> mitos;
	read_cell_data(mitocsv, mitos, dstack.w, dstack.h, dstack.d);
	std::cout << "Read " << mitos.size() << " Mito particles from " << mitocsv << std::endl;

	std::vector<Particle3D> cells;
	read_cell_data(actincsv, cells, dstack.w, dstack.h, dstack.d);
	std::cout << "Read " << cells.size() << " Cells from " << actincsv << std::endl;
	
	std::vector<Point3D> cell_centers(cells.size());
	for (size_t j=0; j<cells.size(); j++) {
		cell_centers[j] = cells[j].center_mass();
	}

	std::vector<Particle3D> nuclei;
	read_cell_data(dnacsv, nuclei, dstack.w, dstack.h, dstack.d);
	std::cout << "Read " << nuclei.size() << " Nuclei from " << dnacsv << std::endl;
	
	std::vector<Point3D> nuc_centers(nuclei.size());
	for (size_t j=0; j<nuclei.size(); j++) {
		nuc_centers[j] = nuclei[j].center_mass();
	}

	std::cout << "Write Mito data to " << datacsv << std::endl;

	wr.append("ID");
	wr.append("XC");
	wr.append("YC");
	wr.append("ZC");
	wr.append("Volume");
	wr.append("NucleusCenterDist");
	wr.append("CellCenterDist");
	wr.next();
	
	int id = 1;
	for (Particle3D& mito : mitos) {
		double xc, yc, zc;
		double mass = dstack.centerMass(mito, &xc, &yc, &zc);
		double vol = mito.volume() * scx * scy * scz;
		
		wr.append(id);
		wr.append(xc);
		wr.append(yc);
		wr.append(zc);
		wr.append(vol);
		
		double nuc_dist = mito_dist_to_center(dstack, mito, nuclei, nuc_centers, scx, scy, scz);
		if (nuc_dist >= 1e6)
			wr.append("NaN");
		else
			wr.append(nuc_dist);
		double cell_dist = mito_dist_to_center(dstack, mito, cells, cell_centers, scx, scy, scz);
		if (cell_dist >= 1e6)
			wr.append("NaN");
		else
			wr.append(cell_dist);

		wr.next();
		++id;
	}
	
	wr.close();
}

//-------- DSP (and other similar Organelles) -------

void detect_dsp(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile,
		int min_volume, double min_snr, double sensitivity)
{
	DSP_Detector dsp(data3d, zd3d, hd3d, wd3d, mask3d);
	dsp.min_volume = (long long)min_volume;
	dsp.min_snr = min_snr;
	dsp.sensitivity = sensitivity;
	
	std::vector<Particle3D> particles = dsp.detect_particles();
	
	unsigned short otsu = dsp.srcstack.otsu();
	double sc = 16. / otsu;
	unsigned short *p = dsp.srcstack.buf;
	for (long long i=0; i<dsp.srcstack.len; i++) {
		double v = (*p) * sc;
		if (v > 255.) v = 255.;
		*p++ = (unsigned short)(v);
	}
	
	size_t valid = 0;
	for (Particle3D& pt : particles)
		if (!pt.empty()) ++valid;
	std::cout << "Write " << valid << " objects to " << csvfile << std::endl;
	write_cell_data(particles, csvfile);
}

void detect_dsp_upd(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d)
{
	DSP_Detector dsp(data3d, zd3d, hd3d, wd3d, mask3d);
	dsp.detect_upd();
}


DSPResult analyze_dsp(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *dspcsv, const char *actincsv,
		const char *outdspcsv, const char *outcellcsv)
{
	DSPResult res;
	
	
	DSP_Analyzer da(data3d, zd3d, hd3d, wd3d, dspcsv, actincsv);
	std::cout << "Loaded " << da.cells.size() << " cells, " << da.dsps.size() << " ORGs." << std::endl;
	
	res.pix_count = da.pix_distances(res.pix_dist_inside, res.pix_dist_outside, 20);
	
	std::vector<DSPtoCell> dsptocells = da.dsp_to_cell_match(20, -2.);
	
	CsvWriter dspwr(outdspcsv);
	if (dspwr.is_open()) {
		std::cout << "Write: " << outdspcsv << std::endl;
		dspwr.SetDoubleFormat("%0.3lf");
		dspwr.append("ORG_ID");
		dspwr.append("Cell_ID");
		dspwr.append("CenterX");
		dspwr.append("CenterY");
		dspwr.append("CenterZ");
		dspwr.append("MinZ");
		dspwr.append("MaxZ");
		dspwr.append("Value");
		dspwr.append("Volume");
		dspwr.append("Overlay");
		dspwr.append("PixToBorder");
		dspwr.next();
		for (DSPtoCell& dc : dsptocells) {
			dspwr.append(dc.dspidx+1);
			dspwr.append((dc.cellidx >=0) ? dc.cellidx+1 : -1);
			dspwr.append(dc.xc);
			dspwr.append(dc.yc);
			dspwr.append(dc.zc);
			dspwr.append(dc.zmin);
			dspwr.append(dc.zmax);
			dspwr.append(dc.value);
			dspwr.append(dc.volume);
			dspwr.append(dc.overlay);
			dspwr.append(dc.border_dist);
			dspwr.next();
		}
		dspwr.close();
	}
	CsvWriter cellwr(outcellcsv);
	if (cellwr.is_open()) {
		std::cout << "Write: " << outcellcsv << std::endl;
		cellwr.append("Cell_ID");
		cellwr.append("CenterX");
		cellwr.append("CenterY");
		cellwr.append("CenterZ");
		cellwr.append("Volume");
		cellwr.next();
		for (int cellidx=0; size_t(cellidx)<da.cells.size(); cellidx++) {
			Particle3D& cell = da.cells[cellidx];
			Point3D& ptc = da.centroids[cellidx];
			cellwr.append(cellidx+1);
			cellwr.append(ptc.x);
			cellwr.append(ptc.y);
			cellwr.append(ptc.z);
			cellwr.append(cell.volume());
			cellwr.next();
		}
		cellwr.close();
	}
	
	da.paint_dsp_mask(dsptocells);
	for (long long i=0; i<da.srcstack.len; i++) {
		da.srcstack.buf[i] = (unsigned short)(da.mstack.buf[i]);
	}
	
	return res;
}

static void rescale_x5(Raster8& msk, Raster8& cmsk)
{
	for (int y=0; y<msk.h; y++) {
		unsigned char *sp = msk.scanLine(y);
		unsigned char *tp0 = cmsk.scanLine(y*5);
		unsigned char *tp = tp0;
		for (int x=0; x<msk.w; x++) {
			unsigned char c = *sp++;
			*tp++ = c;
			*tp++ = c;
			*tp++ = c;
			*tp++ = c;
			*tp++ = c;
		}
		for (int ty=y*5+1; ty<y*5+5; ty++) {
			memcpy(cmsk.scanLine(ty), tp0, cmsk.w * sizeof (unsigned char));
		}
	}
}

void colorize_dsp_x5(unsigned char *mask, int hm, int wm,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		int z, const char *actincsv, const char *indspcsv)
{
	Raster8 msk(wm, hm, mask);
	Raster3D rgb(wm3d, hm3d, zm3d, mask3d);
	Raster8 red = rgb.getPlane(0);
	Raster8 green = rgb.getPlane(1);
	Raster8 blue = rgb.getPlane(2);
	
	for (long long i=0; i<msk.len; i++) {
		if (msk.buf[i] > 0xBF) msk.buf[i] = 0xBF;
	}

	Raster8 amsk(wm, hm, NULL);
	amsk.fill(0);
	
	std::vector<Particle3D> cells;
	read_cell_data(actincsv, cells, msk.w, msk.h, 50);
	std::cout << "Read " << cells.size() << " cells from " << actincsv << std::endl;
	
	std::vector<Point3D> cell_cms(cells.size());
	for (size_t idx=0; idx<cells.size(); idx++) {
		Particle3D& cell = cells[idx];
		cell_cms[idx] = cell.center_mass();
		
		if (z >= 50) continue;
		amsk.paintParticleFill(cell.fills[z], 0x20);
	}
	amsk.expandBorders(0, 0x40, HOOD_SIZE_NEUMANN, 0x20);
	
	rescale_x5(msk, red);
	rescale_x5(msk, green);
	for (long long i=0; i<msk.len; i++) {
		msk.buf[i] += amsk.buf[i];
	}
	rescale_x5(msk, blue);
	
	std::vector<DSPtoCell> dsps;
	CsvReader drdr(indspcsv);
	if (drdr.next()) {
		CsvRow hdrs = drdr.GetRow();
		drdr.setHeaders(hdrs);
		
		// ORG_ID,Cell_ID,CenterX,CenterY,CenterZ,MinZ,MaxZ,Value,Volume,Overlay,PixToBorder
		int i_dsp_id = drdr.header_index("ORG_ID");
		int i_cell_id = drdr.header_index("Cell_ID");
		int i_cx = drdr.header_index("CenterX");
		int i_cy = drdr.header_index("CenterY");
		int i_cz = drdr.header_index("CenterZ");
		int i_minz = drdr.header_index("MinZ");
		int i_maxz = drdr.header_index("MaxZ");
		int i_dist = drdr.header_index("PixToBorder");
		
		while(drdr.next()) {
			int dspid = drdr.GetInt(i_dsp_id, -1);
			if (dspid <= 0) continue;
			DSPtoCell dc;
			dc.dspidx = dspid - 1;
			dc.cellidx = drdr.GetInt(i_cell_id, -1) - 1;
			dc.zmin = drdr.GetInt(i_minz);
			dc.zmax = drdr.GetInt(i_maxz);
			dc.xc = drdr.GetDouble(i_cx);
			dc.yc = drdr.GetDouble(i_cy);
			dc.zc = drdr.GetDouble(i_cz);
			dc.border_dist = drdr.GetDouble(i_dist);
			dsps.push_back(dc);
		}
	}
	std::cout << "Read " << dsps.size() << " ORGs from " << indspcsv << std::endl;
	
	Boundary bnd = red.getBoundary();
	for (DSPtoCell& dc : dsps) {
		if (z < dc.zmin || z > dc.zmax) continue;

		int xc = int(dc.xc * 5. + 2.5);
		int yc = int(dc.yc * 5. + 2.5);
		int zc = int(dc.zc);
		Raster8* pmsk = &green;
		if (dc.border_dist < 0) pmsk = &red;
		std::vector<Point> outer = circle_path(xc, yc, 20);
		for (Point& p : outer) {
			if (bnd.IsInside(p))
				pmsk->setValue(p.x, p.y, 0xFF);
		}
		if (z == zc) {
			std::vector<HSeg> inner = circle_fill(xc, yc, 3);
			pmsk->paintParticleFill(inner, 0xFF);
		}
		
		if (dc.cellidx < 0 || size_t(dc.cellidx) >= cells.size()) continue;
		Particle3D& cell = cells[dc.cellidx];
		Point3D& cm3d = cell_cms[dc.cellidx];
		Point cm = Point(int(cm3d.x*5. + 2.5), int(cm3d.y*5. + 2.5));
		if (!cell.fills[z].empty()) {
			double xm, ym;
			fill_centroid(cell.fills[z], &xm, &ym);
			cm = Point(int(xm*5.), int(ym*5.));
		}
		std::vector<Point> path = straight_line_path(xc, yc, cm.x, cm.y);
		for (Point& p : path) {
			if (!bnd.IsInside(p)) break;
			double dist = p.dist(xc, yc);
			if (dist < 20.) continue;
			if (dist > 50.) break;
			pmsk->setValue(p.x, p.y, 0xFF);
		}
	}
	
}





