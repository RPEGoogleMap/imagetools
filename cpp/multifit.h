#ifndef multifit_h
#define multifit_h

#include "nrc/nr3.h"

#define IMIN(a,b) ((a)<(b)?(a):(b))

/*
Given a matrix a[0..n-1][0..n-1], this routine replaces it by
the LU decomposition of a rowwise permutation of itself.

a and n are input.

a is output, arranged as LU;

indx[0..n-1] is an output vector that records the row permutation
effected by the partial pivoting;

d is output as +/-1 depending on whether the number of row interchanges
was even or odd, respectively.

_vv is an array of size n for temporary storage; It can be NULL, in this
case it is allocated and deallocated internally.

This routine is used in combination with lu_backsubstitution() to solve
linear equations or invert a matrix.
*/
int lu_decomposition(double *a, int n, int *indx, double *d, double *_vv=NULL);

/*
Solves the set of n linear equations A*X = B.

Here a[0..n-1][0..n-1] is input, not as the matrix A but rather
as its LU decomposition, determined by the routine lu_decomposition().

indx[0..n-1] is input as the permutation vector returned by lu_decomposition().

b[0..n-1] is input as the right-hand side vector B, and returns
with the solution vector X.

a, n, and indx are not modified by this routine and can be left in place
for successive calls with different right-hand sides b.

This routine takes into account the possibility that b will begin with
many zero elements, so it is efficient for use in matrix inversion.
*/
void lu_backsubstitution(double *a, int n, int *indx, double *b);

/*
Returns in c[1..np] a set of Savitzky-Golay filter coeffcients. nl is the number of leftward (past)
data points used, while nr is the number of rightward (future) data points, making the total
number of data points used nl+nr+1. ld is the order of the derivative desired (e.g., ld = 0
for smoothed function). m is the order of the smoothing polynomial, also equal to the highest
conserved moment; usual values are m = 2 or m = 4.

Note that the derivatives (ld > 0) obtained after applying S.-G. filter with coefficients
returned by this function, must be multiplied by the factor {-1^ld * ld!}
("!" here is "factorial"), e.g.:
ld=1 : -1
ld=2 : 2
ld=3 : -6
ld=4 : 24
ld=5 : -120
ld=6 : 720
*/
int savgol(double *c, int np, int nl, int nr, int ld, int m);


class PolynomialFit
{
private:
	int dim;
	double *coeffs;
public:
	PolynomialFit(int ord) { dim = ord + 1; coeffs = new double[dim]; }
	virtual ~PolynomialFit() { delete[] coeffs; }
	void fit(int nsamp, double *x, double *y);
	double value(double x);
	double size() { return dim; }
	double(operator [])(int idx) { return coeffs[idx]; }
};

class SampleIterator2D
{
public:
	virtual int GetNumSamples() = 0;
	virtual void reset() = 0;
	virtual bool next(double *px, double *py, double *pz) = 0;
};

inline int ord2dim2D(int ord) { return (ord + 1)*(ord + 2) / 2; }

class PolynomialFit2D
{
private:
	int ord;
	int dim;
	double *coeffs;
	double *pwx, *pwy;
public:
	PolynomialFit2D(int ord) {
		this->ord = ord;
		dim = ord2dim2D(ord);
		coeffs = new double[dim];
		pwx = new double[ord + 1];
		pwy = new double[ord + 1];
	}
	virtual ~PolynomialFit2D() {
		delete[] coeffs;
		delete[] pwx;
		delete[] pwy;
	}
	void fit(int nsamp, double *x, double *y, double *z);
	void fit(SampleIterator2D& iter);
	double value(double x, double y);
	double size() { return dim; }
	double(operator [])(int idx) { return coeffs[idx]; }
};

/*
Compute Savitzky-Golay filter coefficients for a 2D neighborhood,
given the order of the approximating polynom.

nbs = half neighborhood size; defines a square neighborhood from (x-nbs,y-nbs) to (x+nbs,y+nbs).
	E.g. if nbs=4, the neighborhood is +/-4 in all directions, total size is 9*9.
ord = polynomial order (for both X and Y). The number of polynomial terms is (ord+1)*(ord+2)/2,
	e.g. if ord=3, the polynomial has 10 terms, (1, y, x, y*y, y*x, x*x, y*y*y, y*y*x, y*x*x, x*x*x).

Method GetCoeffs() copies Savitzky-Golay coefficients for the requested derivative ld
into a user-provided array c. The array must be at least NbSize() size.
The ld parameter can range from 0 to MaxDerOrd().
	ld=0 - (smoothed) function itself;
	ld=1 - partial y first derivative d/dy;
	ld=2 - partial x first derivative d/dx;
	ld=3 - partial y second derivative d^2/dy^2;
	ld=4 - partial x-y cross-derivative d^2/dx*dy
	ld=5 - partial x second derivative d^2/dx^2
	ld=6 - partial y third derivative d^3/dy^3
	... etc.
The coefficients in the array c are arranged as
	[(x-nbs,y-nbs)...(x,y-nbs)...(x+nbs,y-nbs)]
	[(x-nbs,y)...(x,y)...(x+nbs,y)]
	[(x-nbs,y+nbs)...(x,y+nbs)...(x+nbs,y+nbs)]
*/

class SavitzkyGolay2D
{
private:
	int nbs, ord;
	int wsz, ndat, ma;
	MatDoub aa;
public:
	SavitzkyGolay2D(int nbs, int ord);
	virtual ~SavitzkyGolay2D() {}
	int GetNeighb() { return nbs; }
	int GetOrder() { return ord; }
	int NbOneSide() { return wsz; }
	int NbSize() { return ndat; }
	int MaxDerOrd() { return ma - 1; }
	bool GetCoeffs(double *c, int ld)
	{
		if (ld < 0 || ld >= ma) return false;
		for (int i = 0; i < ndat; i++) c[i] = aa[ld][i];
		return true;
	}
	double *CoeffPtr(int ld) { return aa[ld]; }
};

#endif
