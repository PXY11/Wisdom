#include"dataStructure.h"

//---------------------dvec�ڲ���Ա��������------------------------------
dvec::dvec(const Array1D& a) : tvec<double>(static_cast<int>(a.size())) {
	for (int i = 0; i < n; i++) (*this)[i] = a[i];
}

dvec::dvec(const Array2D& a) : tvec<double>(static_cast<int>(a.size()* a[0].size())) {
	int		Nr = static_cast<int>(a.size());
	int		Nc = static_cast<int>(a[0].size());
	int		i, j;
	for (i = 0; i < Nr; i++) {
		for (j = 0; j < Nc; j++) (*this)[Nc * i + j] = a[i][j];
	}
}

double dvec::sum()
{
	double sum0 = 0;
	for (int i = 0; i < n; i++) sum0 = sum0+(*this)[i];
	return sum0;
}

//---------------------ivec�ڲ���Ա��������------------------------------
ivec::ivec(const Array1D& a) : tvec<int>(static_cast<int>(a.size())) {
	for (int i = 0; i < n; i++) (*this)[i] = (int)a[i];
}

ivec::ivec(const Array2D& a) : tvec<int>(static_cast<int>(a.size()* a[0].size())) {
	int		Nr = static_cast<int>(a.size());
	int		Nc = static_cast<int>(a[0].size());
	int		i, j;
	for (i = 0; i < Nr; i++) {
		for (j = 0; j < Nc; j++) (*this)[Nc * i + j] = a[i][j];
	}
}

//---------------------cvec�ڲ���Ա��������------------------------------
void cvec::append(const char* s) {
	int		l0 = p ? static_cast<int>(strlen(p)) : 0; //ԭ�ַ�������
	int		l1 = static_cast<int>(strlen(s));//��Ҫ��ӵ��ַ�������
	l0 = MAX(0, MIN(n, l0)); //l0,nһ���������ͬ��
	l0 = MAX(0, MIN(MAX_STR / 2 - 1, l0));
	l1 = MAX(0, MIN(MAX_STR / 2 - 1, l1));
	resize(l0 + l1 + 1); //����չcvec�����е��ַ����ĳ��ȣ���ȱλ���Ȳ�0
	strncpy(p + l0, s, l1); //����Դ�ַ���s�е�l1�����ȵ��ַ������Ƶ�p�ַ�����λ��Ϊl0��ʼ�ĵط�������ĩβ��\0;
	p[l0 + l1] = '\0';

}

//-------------------imat�ڲ���Ա��������-------------------------------
void imat::resize(int r0, int c0) {
	int		i, j;
	int* p0 = new int[static_cast<size_t>(r0*c0)];
	for (i = 0; i < MIN(r, r0); i++) {
		for (j = 0; j < MIN(c, c0); j++) p0[c0 * i + j] = (*this)[i][j];//��ͬ��λ�ñ���ԭֵ
		for (j = MIN(c, c0); j < c0; j++) p0[c0 * i + j] = 0;//�����У�����ֵ0�������У���ִ��
	}
	for (i = MIN(r, r0); i < r0; i++) for (j = 0; j < c0; j++) p0[c0 * i + j] = 0;//�����У�����ֵ0����С�У���ִ��
	delete[] p;
	p = p0;
	r = r0;
	c = c0;
}

imat& imat::operator=(const imat& mat)
{
	//int r;
	//int c;
	resize(mat.r, mat.c);
	for (int i = 0; i < mat.r; i++)
		for (int j = 0; j < mat.c; j++)(*this)[i][j] = mat[i][j];// p[i * mat.r + j] = mat.p[i * mat.r + j];
	return *this;
}

imat::imat(const imat& mat) : r(mat.r), c(mat.c), p(new int[static_cast<size_t>(mat.r*mat.c)])
{
	for (int i = 0; i < mat.r; i++)
		for (int j = 0; j < mat.c; j++)
			//p[i * mat.c + j] = mat.p[i * mat.c + j];
			(*this)[i][j] = mat[i][j];//�������ַ�ʽ������
}


//-------------------dmat�ڲ���Ա��������--------------------------------
dmat::dmat(const Array2D& a) : r(static_cast<int>(a.size())), c(r <= 0 ? 0 : static_cast<int>(a[0].size())), p(new double[r * c]) {
	for (int i = 0; i < r; i++) for (int j = 0; j < c; j++) (*this)[i][j] = a[i][j]; //this����ǰclass�ľ����ĳ��������ܹ�ʹ��[i][j]
};

void dmat::resize(int r0, int c0) {
	int		i, j;
	double* p0 = new double[r0 * c0];
	for (i = 0; i < MIN(r, r0); i++) {
		for (j = 0; j < MIN(c, c0); j++) p0[c0 * i + j] = (*this)[i][j];//��ͬ��λ�ñ���ԭֵ
		for (j = MIN(c, c0); j < c0; j++) p0[c0 * i + j] = 0;//�����У�����ֵ0�������У���ִ��
	}
	for (i = MIN(r, r0); i < r0; i++) for (j = 0; j < c0; j++) p0[c0 * i + j] = 0;//�����У�����ֵ0����С�У���ִ��
	delete[] p;
	p = p0;
	r = r0;
	c = c0;
}

// i -> i*r mod r*c-1.
void dmat::transpose(void) {
	double	temp;
	if (r > 1 && c > 1) {
		int		n = r * c - 1;
		int		i;
		dvec	q(r * c);
		for (i = 0; i < r * c; i++) q[i] = p[i];
		for (i = 1; i < n; i++)
			p[i] = q[(i * c) % n];
	}
	temp = r; r = c; c = temp;
}

dmat& dmat::operator=(const dmat& mat)
{
	//int r;
	//int c;
	resize(mat.r, mat.c);
	for (int i = 0; i < mat.r; i++)
		for (int j = 0; j < mat.c; j++)(*this)[i][j] = mat[i][j];// p[i * mat.r + j] = mat.p[i * mat.r + j];
	return *this;
}

dmat::dmat(const dmat& mat) : r(mat.r), c(mat.c), p(new double[static_cast<size_t>(mat.r * mat.c)])
{
	for (int i = 0; i < mat.r; i++)
		for (int j = 0; j < mat.c; j++)
			//p[i * mat.c + j] = mat.p[i * mat.c + j];
			(*this)[i][j] = mat[i][j];//�������ַ�ʽ������
}


//====================================================================================
/* allocate a double matrix */
double** New2DMatrix(long nrow, long ncol)
{
	double** m;

	/* allocate pointers to rows */
	m = new double* [nrow];

	/* allocate rows and set pointers to them */
	m[0] = new double[nrow * ncol];
	for (long i = 1; i < nrow; i++)
	{
		m[i] = m[i - 1] + ncol;
	}
	return m;
}

void Delete2DMatrix(double** m)
{
	delete[] m[0];
	delete[] m;
}





//-----------------dmatv�ڲ���Ա��������--------------------------------
dmatV::dmatV(int n_, const Array2D** mat) : n(n_), p(new const dmat* [n]) {
	//����n�����ȵĿռ䣬ÿ�������dmat��ָ�룬Ҳ���Ƿ�dmat�����ַ
	//p[i]�����ŵ��ǵ�ַ��˭�ĵ�ַ�أ�һ��dmat����ĵ�ַ��p[i] = new dmat(*mat[i])//*mat[i]��ʾ����Array2D����
	int		i, j;
	for (i = 0; i < n; i++) {
		if (!mat || !mat[i]) { p[i] = 0; continue; }
		int	done = 0;
		for (j = 0; !done && j < i; j++) {
			if (mat[j] == mat[i]) {// ��ʾ�������mat������������ַָ����ͬһ��Array2D����Ļ����Ͳ��ÿ����¿ռ��ˣ������ú����p[i]ָ��ǰ��p[j]ָ�����ͬ�Ķ���
				p[i] = p[j];
				done = 1;
			}
		}
		if (done) continue;
		//		dmat	Mi(*mat[i]);	// For UNIX compiler 13-jul-04.
		//		p[i]= new dmat(Mi);
		p[i] = new dmat(*mat[i]);//����dmat�Ĺ��캯����
	}
}

dmatV::~dmatV(void) {
	int		i, j;
	if (p) {
		for (i = 0; i < n; i++) {
			if (!p[i]) continue;
			int	later = 0;
			for (j = i + 1; !later && j < n; j++) if (p[j] == p[i]) later = 1;
			if (later) continue;
			delete p[i];//��Ϊ�ڹ��캯����ʹ����������ε�new,���������������У�ͬ����Ҫ��������ε�delete
		}
	}
	delete[] p;
}