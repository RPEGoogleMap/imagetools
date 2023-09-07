#include "fstack.h"

bool comp_nbrdist3d(const NbrDist3D& a, const NbrDist3D& b)
{
	if (a.dz < b.dz) return true;
	if (a.dz > b.dz) return false;
	if (a.dist < b.dist) return true;
	if (a.dist > b.dist) return false;
	if (a.dy < b.dy) return true;
	if (a.dy > b.dy) return false;
	if (a.dx < b.dx) return true;
	return false;
};

std::vector<NbrDist3D> NbrDist3D::list_dist_3d(double xsc, double ysc, double zsc, double maxr)
{
	std::vector<NbrDist3D> res;
	
	for (int z=-5; z<=5; z++) {
		double dz = z * zsc;
		for (int y=-5; y<=5; y++) {
			double dy = y * ysc;
			for (int x=-5; x<=5; x++) {
				double dx = x * xsc;
				if (dx == 0 && dy == 0 && dz == 0) continue;
				double dist = sqrt(dz*dz + dy*dy + dx*dx);
				if (dist <= maxr)
					res.push_back(NbrDist3D(x, y, z, float(dist)));
			}
		}
	}
	std::sort(res.begin(), res.end(), comp_nbrdist3d);
	
	return res;
}

void distance_map(RasterF3D<float>& dstack, std::vector<Particle3D>& particles, std::vector<Particle3D>& cells)
{
	dstack.fill(-1.);
	
	int zmin = dstack.d - 1;
	int zmax = 0;
	for (Particle3D& cell : cells) {
		dstack.paintParticle(cell, 1000000.);
		if (cell.bnd.zmin < zmin) zmin = cell.bnd.zmin;
		if (cell.bnd.zmax > zmax) zmax = cell.bnd.zmax;
	}
	for (Particle3D& pt : particles)
		dstack.paintParticle(pt, 0.);
	std::cout << "zmin=" << zmin << " zmax=" << zmax << std::endl;
	Boundary clip(0, 0, dstack.w-1, dstack.h-1);

	double rad = sqrt(dstack.scx*dstack.scx + dstack.scy*dstack.scy + dstack.scz*dstack.scz) * 1.2;
	std::vector<NbrDist3D> nbd = NbrDist3D::list_dist_3d(dstack.scx, dstack.scy, dstack.scz, rad);
	std::cout << "Nbr size: " << nbd.size() << std::endl;
	
	RasterF3D<float>* gstack = new RasterF3D<float>(dstack.w, dstack.h, dstack.d, NULL);
	
	int iter;
	for (iter=0; iter<100; iter++) {
		memcpy(gstack->buf, dstack.buf, gstack->len * sizeof(float));
		
		long long nmod = 0ll;
		
		for (int z=zmin; z<=zmax; z++) {
			int j0=-1, j1=-1;
			for (int j=0; size_t(j)<nbd.size(); j++) {
				int z1 = z + nbd[j].dz;
				if (j0<0 && z1>=0) j0 = j;
				if (z1 >= dstack.d) break;
				j1 = j;
			}
			
			for (int y=0; y<dstack.h; y++) {
				for (int x=0; x<dstack.w; x++) {
					float sdist = dstack.value(x, y, z);
					if (sdist <= 0.) continue;
					float cdist = sdist;
					for (int j=j0; j<=j1; j++) {
						NbrDist3D& nb = nbd[j];
						int xx = x + nb.dx;
						int yy = y + nb.dy;
						if (!clip.IsInside(xx, yy)) continue;
						float tv = dstack.value(xx, yy, z+nb.dz);
						if (tv < 0.) continue;
						float _cdist = tv + nb.dist;
						if (_cdist < cdist) cdist = _cdist;
					}
					if (cdist < sdist) {
						gstack->setValue(x, y, z, cdist);
						++nmod;
					}
				}
			}
			
		}
		std::cout << "iter " << iter << " nmod=" << nmod << std::endl;
		if (nmod == 0) break;
		memcpy(dstack.buf, gstack->buf, gstack->len * sizeof(float));
	}
	
	std::cout << "iter=" << iter << std::endl;
	
	delete gstack;
}

