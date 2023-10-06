class FitBasis
{
public:
	virtual Int size() = 0;
	virtual VecDoub basis(Doub x) = 0;
};

class FitBasisMD
{
public:
	virtual Int size() = 0;
	virtual VecDoub basis(VecDoub_I &xx) = 0;
};

struct Fitsvd {
	Int ndat, ma;
	Doub tol;
	VecDoub_I *x,&y,&sig;
	// VecDoub (*funcs)(const Doub);
	FitBasis *fb;
	VecDoub a;
	MatDoub covar;
	Doub chisq;

	Fitsvd(VecDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig,
		// VecDoub funks(const Doub),
		FitBasis *_fb,
		const Doub TOL=1.e-12)
	: ndat(yy.size()), x(&xx), xmd(NULL), y(yy), sig(ssig),
		// funcs(funks),
		fb(_fb),
		tol(TOL) {}

	void fit() {
		Int i,j,k;
		Doub tmp,thresh,sum;
		if (x) ma = fb->size(); // funcs((*x)[0]).size();
		else ma = fbmd->size(); // funcsmd(row(*xmd,0)).size();
		a.resize(ma);
		covar.resize(ma,ma);
		MatDoub aa(ndat,ma);
		VecDoub b(ndat),afunc(ma);
		for (i=0;i<ndat;i++) {
			if (x) afunc = fb->basis((*x)[i]); // funcs((*x)[i]);
			else afunc = fbmd->basis(row(*xmd, i)); // funcsmd(row(*xmd,i));
			tmp=1.0/sig[i];
			for (j=0;j<ma;j++) aa[i][j]=afunc[j]*tmp;
			b[i]=y[i]*tmp;
		}
		SVD svd(aa);
		thresh = (tol > 0. ? tol*svd.w[0] : -1.);
		svd.solve(b,a,thresh);
		chisq=0.0;
		for (i=0;i<ndat;i++) {
			sum=0.;
			for (j=0;j<ma;j++) sum += aa[i][j]*a[j];
			chisq += SQR(sum-b[i]);
		}
		for (i=0;i<ma;i++) {
			for (j=0;j<i+1;j++) {
				sum=0.0;
				for (k=0;k<ma;k++) if (svd.w[k] > svd.tsh)
					sum += svd.v[i][k]*svd.v[j][k]/SQR(svd.w[k]);
				covar[j][i]=covar[i][j]=sum;
			}
		}

	}

	MatDoub_I *xmd;
	// VecDoub (*funcsmd)(VecDoub_I &);
	FitBasisMD *fbmd;

	Fitsvd(MatDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig,
		//VecDoub funks(VecDoub_I &),
		FitBasisMD *_fbmd,
		const Doub TOL=1.e-12)
	: ndat(yy.size()), x(NULL), xmd(&xx), y(yy), sig(ssig),
		//funcsmd(funks),
		fbmd(_fbmd),
		tol(TOL) {}

	VecDoub row(MatDoub_I &a, const Int i) {
		Int j,n=a.ncols();
		VecDoub ans(n);
		for (j=0;j<n;j++) ans[j] = a[i][j];
		return ans;
	}
};
