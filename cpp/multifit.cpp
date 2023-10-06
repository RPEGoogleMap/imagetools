#include "nrc/nr3.h"
#include "nrc/svd.h"
#include "nrc/fitsvd.h"
#include "multifit.h"

#define TINY 1.0e-20

static double int_power(double x, int pow)
{
	double r = 1.;
	if (pow <= 0) return r;
	r = x;
	while (--pow > 0) r *= x;
	return r;
}

int lu_decomposition(double *a, int n, int *indx, double *d, double *_vv)
{
	int i, imax, j, k;
	double big, dum, sum, temp;
	double *vv = _vv;	// vv stores the implicit scaling of each row.
	int rc = 0;

	if (!_vv) vv = new double[n];  // (double *)my_malloc(n * sizeof(double));
	*d = 1.0;	// No row interchanges yet.

	for (i = 0; i < n; i++)
	{
		// Loop over rows to get the implicit scaling information.
		big = 0.0;
		for (j = 0; j < n; j++)
		{
			temp = fabs(a[i*n + j]);
			if (temp > big) big = temp;
		}
		if (big == 0.0)
		{
			// No nonzero largest element.
			rc = -1;
			goto zerr;
		}
		vv[i] = 1.0 / big;	// Save the scaling.
	}

	// This is the loop over columns of Crout's method.
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < j; i++)
		{
			sum = a[i*n + j];
			for (k = 0; k < i; k++) sum -= a[i*n + k] * a[k*n + j];
			a[i*n + j] = sum;
		}
		big = 0.0;	// Initialize for the search for largest pivot element.
		for (i = j; i < n; i++)
		{
			sum = a[i*n + j];
			for (k = 0; k < j; k++)
				sum -= a[i*n + k] * a[k*n + j];
			a[i*n + j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big)
			{
				// Is the figure of merit for the pivot better than the best so far?
				big = dum;
				imax = i;
			}
		}
		if (j != imax)
		{
			// Do we need to interchange rows?
			for (k = 0; k < n; k++)
			{
				// Yes, do so...
				dum = a[imax*n + k];
				a[imax*n + k] = a[j*n + k];
				a[j*n + k] = dum;
			}
			*d = -(*d);		// ...and change the parity of d.
			dum = vv[imax];
			vv[imax] = vv[j];	// Also interchange the scale factor.
			vv[j] = dum;
		}
		indx[j] = imax;

		// If the pivot element is zero the matrix is singular
		// (at least to the precision of the algorithm).
		// For some applications on singular matrices,
		// it is desirable to substitute TINY for zero.
		if (a[j*n + j] == 0.0) a[j*n + j] = TINY;
		if (j != n)
		{
			// Now, finally, divide by the pivot element.
			dum = 1.0 / (a[j*n + j]);
			for (i = j + 1; i < n; i++) a[i*n + j] *= dum;
		}

	}	//Go back for the next column in the reduction.

zerr:
	if (!_vv) delete[] vv;	// my_free(vv);
	return rc;
}

void lu_backsubstitution(double *a, int n, int *indx, double *b)
{
	int i, ii = -1, ip, j;
	double sum;

	// When ii is set to a non-negative value, it will become the
	// index of the first nonvanishing element of b. We now
	// do the forward substitution. The only new wrinkle is
	// to unscramble the permutation as we go.
	for (i = 0; i < n; i++)
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0)
		{
			for (j = ii; j <= i - 1; j++) sum -= a[i*n + j] * b[j];
		}
		else if (sum)
		{
			// A nonzero element was encountered, so from now on we
			// will have to do the sums in the loop above.
			ii = i;
		}
		b[i] = sum;
	}

	for (i = n - 1; i >= 0; i--)
	{
		// Now we do the backsubstitution.
		sum = b[i];
		for (j = i + 1; j < n; j++) sum -= a[i*n + j] * b[j];
		b[i] = sum / a[i*n + i];	// Store a component of the solution vector X.
	}	// All done!
}

int savgol(double *c, int np, int nl, int nr, int ld, int m)
{
	int imj, ipj, j, k, kk, mm, *indx, dim = m + 1;
	double d, fac, sum, *a, *b;

	if (np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m)
		return -1;	// Bad arguments

	indx = new int[dim];
	a = new double[dim*dim];
	b = new double[dim];

	// Set up the normal equations of the desired least-squares fit.
	for (ipj = 0; ipj <= (m << 1); ipj++)
	{
		// std::cout << "ipj=" << ipj << std::endl;
		sum = (ipj ? 0.0 : 1.0);
		for (k = 1; k <= nr; k++) sum += int_power((double)k, ipj);
		for (k = 1; k <= nl; k++) sum += int_power((double)-k, ipj);
		mm = IMIN(ipj, 2 * m - ipj);
		for (imj = -mm; imj <= mm; imj += 2)
			a[((ipj + imj) / 2)*dim + (ipj - imj) / 2] = sum;
	}

/*
	std::cout << "-- a --" << std::endl;
	for (int j = 0; j < dim; j++) {
		for (int i = 0; i < dim; i++)
			std::cout << a[j*dim + i] << " ";
		std::cout << std::endl;
	}
*/

	// Solve them: LU decomposition.
	lu_decomposition(a, dim, indx, &d, NULL);

	// Right-hand side vector is unit vector,
	// depending on which derivative we want.
	for (j = 0; j < dim; j++) b[j] = 0.0;
	b[ld] = 1.0;

	// Get one row of the inverse matrix.
	lu_backsubstitution(a, dim, indx, b);

	// Zero the output array (it may be bigger than number of coefficients).
	for (kk = 0; kk < np; kk++) c[kk] = 0.0;

	for (k = -nl; k <= nr; k++)
	{
		// Each Savitzky-Golay coefficient is the dot product of powers
		// of an integer with the inverse matrix row.
		sum = b[0];
		fac = 1.0;
		for (mm = 0; mm < m; mm++) sum += b[mm + 1] * (fac *= k);
		kk = k + nl; // normal order: -nl, -(nl-1), ... -1, 0, 1, ... (nr-1), nr
					 // To store in wrap-around order: kk = ((np-k) % np);
		c[kk] = sum;
	}

	delete[] b;
	delete[] a;
	delete[] indx;
	return 0;
}


class FitBasisPoly : public FitBasis
{
private:
	int dim;
public:
	FitBasisPoly(int dim) { this->dim = dim; }
	virtual Int size() override { return dim; }
	virtual VecDoub basis(double x) override
	{
		VecDoub p(dim);
		p[0] = 1.;
		for (Int j = 1; j < dim; j++)
			p[j] = p[j - 1] * x;
		return p;
	}
};

void PolynomialFit::fit(int nsamp, double *x, double *y)
{
	VecDoub xx(nsamp, x);
	VecDoub yy(nsamp, y);
	double sigma = 1.;
	VecDoub ssig(nsamp, sigma);
	FitBasisPoly fp(dim);
	Fitsvd svd(xx, yy, ssig, &fp);
	svd.fit();
	for (int i = 0; i < dim; i++)
		coeffs[i] = svd.a[i];
}

double PolynomialFit::value(double x)
{
	double res = coeffs[0];
	double xx = 1.;
	for (int i = 1; i < dim; i++) {
		xx *= x;
		res += xx * coeffs[i];
	}
	return res;
}

class FitBasisPoly2D : public FitBasisMD
{
private:
	int ord;
	int dim;
	double *pwx, *pwy;
public:
	FitBasisPoly2D(int ord, double *pwx, double *pwy) {
		this->ord = ord;
		dim = ord2dim2D(ord);
		this->pwx = pwx;
		this->pwy = pwy;
	}
	virtual Int size() override { return dim; }
	virtual VecDoub basis(VecDoub_I &xx)
	{
		VecDoub res(dim);
		double x = xx[0];
		double y = xx[1];
		int o;
		pwx[0] = pwy[0] = 1;
		for (o = 1; o <= ord; o++) {
			pwx[o] = pwx[o - 1] * x;
			pwy[o] = pwy[o - 1] * y;
		}
		int i = 0;
		res[i++] = 1.;
		for (o = 1; o <= ord; o++) {
			for (int j = 0; j <= o; j++) {
				res[i++] = pwx[o - j] * pwy[j];
			}
		}
// std::cout << " dim=" << std::to_string(dim) << " i=" << std::to_string(i) << std::endl;
		return res;
	}
};

void PolynomialFit2D::fit(int nsamp, double *x, double *y, double *z)
{
	MatDoub xx(nsamp, 2);
	VecDoub zz(nsamp);
	double sigma = 1.;
	VecDoub ssig(nsamp, sigma);
	int i;
	for (i = 0; i < nsamp; i++) {
		xx[i][0] = x[i];
		xx[i][1] = y[i];
		zz[i] = z[i];
	}
	FitBasisPoly2D fbmd(ord, pwx, pwy);
	Fitsvd svd(xx, zz, ssig, &fbmd);
	svd.fit();
	for (i = 0; i < dim; i++)
		coeffs[i] = svd.a[i];
}

void PolynomialFit2D::fit(SampleIterator2D& iter)
{
	int nsamp = iter.GetNumSamples();
	MatDoub xx(nsamp, 2);
	VecDoub zz(nsamp);
	double sigma = 1.;
	VecDoub ssig(nsamp, sigma);
	int i;
	iter.reset();
	for (i = 0; i < nsamp; i++) {
		iter.next(&(xx[i][0]), &(xx[i][1]), &(zz[i]));
	}
	FitBasisPoly2D fbmd(ord, pwx, pwy);
	Fitsvd svd(xx, zz, ssig, &fbmd);
	svd.fit();
	for (i = 0; i < dim; i++)
		coeffs[i] = svd.a[i];
}

double PolynomialFit2D::value(double x, double y)
{
	int o;
	pwx[0] = pwy[0] = 1;
	for (o = 1; o <= ord; o++) {
		pwx[o] = pwx[o - 1] * x;
		pwy[o] = pwy[o - 1] * y;
	}
	int i = 0;
	double res = coeffs[i++];
	for (o = 1; o <= ord; o++) {
		for (int j = 0; j <= o; j++) {
			res += coeffs[i++] * pwx[o - j] * pwy[j];
		}
	}
	return res;
}

SavitzkyGolay2D::SavitzkyGolay2D(int nbs, int ord)
{
	this->nbs = nbs;
	this->ord = ord;
	wsz = nbs + nbs + 1;
	ndat = wsz * wsz;
	ma = ord2dim2D(ord);

	aa.resize(ndat, ma);
	int jdat = 0;
	for (int y = -nbs; y <= nbs; y++) {
		for (int x = -nbs; x <= nbs; x++) {
			int idat = 0;
			for (int mm = 0; mm <= ord; mm++) {
				for (int xord = 0; xord <= mm; xord++) {
					aa[jdat][idat++] = int_power(double(x), xord) * int_power(double(y), mm - xord);
				}
			}
			++jdat;
		}
	}
	SVD svd(aa);
	svd.pinv(aa);
}
