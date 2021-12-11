#include"dataStructure.h"

//---------------------dvec内部成员函数定义------------------------------
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

//---------------------ivec内部成员函数定义------------------------------
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

//---------------------cvec内部成员函数定义------------------------------
void cvec::append(const char* s) {
	int		l0 = p ? static_cast<int>(strlen(p)) : 0; //原字符串长度
	int		l1 = static_cast<int>(strlen(s));//需要添加的字符串长度
	l0 = MAX(0, MIN(n, l0)); //l0,n一般情况是相同的
	l0 = MAX(0, MIN(MAX_STR / 2 - 1, l0));
	l1 = MAX(0, MIN(MAX_STR / 2 - 1, l1));
	resize(l0 + l1 + 1); //先扩展cvec对象中的字符串的长度，空缺位置先补0
	strncpy(p + l0, s, l1); //将来源字符串s中的l1个长度的字符，复制到p字符串中位置为l0开始的地方；并在末尾补\0;
	p[l0 + l1] = '\0';

}

//-------------------imat内部成员函数定义-------------------------------
void imat::resize(int r0, int c0) {
	int		i, j;
	int* p0 = new int[static_cast<size_t>(r0*c0)];
	for (i = 0; i < MIN(r, r0); i++) {
		for (j = 0; j < MIN(c, c0); j++) p0[c0 * i + j] = (*this)[i][j];//相同的位置保持原值
		for (j = MIN(c, c0); j < c0; j++) p0[c0 * i + j] = 0;//增加列，均赋值0；缩减列，不执行
	}
	for (i = MIN(r, r0); i < r0; i++) for (j = 0; j < c0; j++) p0[c0 * i + j] = 0;//增加行，均赋值0；减小行，不执行
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
			(*this)[i][j] = mat[i][j];//以上两种方式均可以
}


//-------------------dmat内部成员函数定义--------------------------------
dmat::dmat(const Array2D& a) : r(static_cast<int>(a.size())), c(r <= 0 ? 0 : static_cast<int>(a[0].size())), p(new double[r * c]) {
	for (int i = 0; i < r; i++) for (int j = 0; j < c; j++) (*this)[i][j] = a[i][j]; //this代表当前class的具体的某个类对象。能够使用[i][j]
};

void dmat::resize(int r0, int c0) {
	int		i, j;
	double* p0 = new double[r0 * c0];
	for (i = 0; i < MIN(r, r0); i++) {
		for (j = 0; j < MIN(c, c0); j++) p0[c0 * i + j] = (*this)[i][j];//相同的位置保持原值
		for (j = MIN(c, c0); j < c0; j++) p0[c0 * i + j] = 0;//增加列，均赋值0；缩减列，不执行
	}
	for (i = MIN(r, r0); i < r0; i++) for (j = 0; j < c0; j++) p0[c0 * i + j] = 0;//增加行，均赋值0；减小行，不执行
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
			(*this)[i][j] = mat[i][j];//以上两种方式均可以
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





//-----------------dmatv内部成员函数定义--------------------------------
dmatV::dmatV(int n_, const Array2D** mat) : n(n_), p(new const dmat* [n]) {
	//开辟n个长度的空间，每个里面放dmat的指针，也就是放dmat对象地址
	//p[i]里面存放的是地址，谁的地址呢？一个dmat对象的地址。p[i] = new dmat(*mat[i])//*mat[i]表示的是Array2D对象
	int		i, j;
	for (i = 0; i < n; i++) {
		if (!mat || !mat[i]) { p[i] = 0; continue; }
		int	done = 0;
		for (j = 0; !done && j < i; j++) {
			if (mat[j] == mat[i]) {// 表示的是如果mat的连续两个地址指向了同一个Array2D对象的话，就不用开辟新空间了，而是用后面的p[i]指向前面p[j]指向的相同的对象
				p[i] = p[j];
				done = 1;
			}
		}
		if (done) continue;
		//		dmat	Mi(*mat[i]);	// For UNIX compiler 13-jul-04.
		//		p[i]= new dmat(Mi);
		p[i] = new dmat(*mat[i]);//调用dmat的构造函数，
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
			delete p[i];//因为在构造函数中使用了两个层次的new,所以在析构函数中，同样需要有两个层次的delete
		}
	}
	delete[] p;
}