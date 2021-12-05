#pragma once//ͬһ���ļ����ᱻ�������

#pragma warning(disable:4018)// �з���/�޷��Ų�ƥ�� ����. ��unsigned int ���� �� int ���ͱ���֮��ıȽϡ�
#pragma warning(disable:4244)// ��������ת�Ϳ������ݶ�ʧ�� �署doubleת��Ϊint. ��ʹ��static_cast<int>����ǿ������ת����

//��base�ļ��ο��������ļ�
//SharedDefines.h---������Array1D��


#define NOMINMAX
#include <Windows.h>
//#define NOMINMAX���߱�����(��ʵ������Ԥ������)����min��max�Ķ���,��ֻ����#include "windows.h"֮ǰִ�д˲���ʱ�Ż�Ӧ������

#include <algorithm>//C++��׼�⣺����������vector��list��set��map�������뵽��һЩ����,�����,�滻�����򡢼����ȳ��õĹ���ȫ��������
#include <fstream>//C++�е�IO������
#include <strstream>//C++��C���Է����������
#include <valarray>
#include <float.h> //C���Ա�׼��ͷ�ļ��������˸����͵ķ�Χ�����ȡ�������ֵ������
using namespace std;//valarray��C++ STL��׼�⣻������Ӵ���䣬Ҳ����ʹ��std::valarray

//1. ʹ��#define���峣���ͳ�����ѧ����

//#define  ��ʶ��  ����   //ע��, ���û�зֺ�
//�����ԡ�#����ͷ�ľ�ΪԤ����ָ��,Ԥ���벻�Ǳ��룬���Ǳ���ǰ�Ĵ����������������ʽ����֮ǰ��ϵͳ�Զ���ɵġ�Ԥ������ִ�еĲ������Ǽ򵥵ġ��ı����滻��
// 
//#define��typedef������
//#define�ñ�ʶ����ʾһ���ַ�������Ϊ�ꣻ typedef���ڸ�һ����������µ����֣������Ķ��������ڱ���׶Ρ�

//#define��const���峣��������
//#define��Ԥ����׶��滻��const�ڱ���ʱȷ��ֵ
//#define�����ͣ���������Ͱ�ȫ�ԣ�const���������ͣ�����ʱ���
//#define����õ���ο�����const�ھ�̬�洢�������ڴ棬�������й������ڴ���ֻ��һ������


#define MAX(a,b)	( (a)<(b) ? (b) : (a) )
#define MIN(a,b)	( (a)<(b) ? (a) : (b) )
#define MIN4(a,b,c,d)	( MIN( MIN((a),(b)) , MIN((c),(d)) ) ) //rdv.cpp��ʹ��
#define ABS(a)		( (a)< 0  ?-(a) : (a) )
#define SGN(a)		( (a)< 0  ? -1  :  1  )
#define SWAP(a,b)   {temp=(a); (a)=(b); (b)=temp;} //ԭ�ȶ�����NumUtilities.cpp��
#define EPS         (1.e-8 )
#define EPS1		( 1e-8 )
#define MAX_STR	    (1<<12)		// Maximum string length. 2**12=4096. ��ʾ��������1�����ƶ���λ��
#define JMAX        200//ԭ�ȶ�����NumUtilities.cpp��ͷ��
#define	ROOT_2PI	( 2.506628274631 )//ԭ�ȶ�����util.h; ��CumNorm.h��ͬ�������ֵ��#define DBA_SQRT2PI
#define MAXN		(1000)  //����������lu.cpp�ж����
#define	MM(i,j)		(M[(i)*n+(j)])
#define	STATUS_MAX	(250)	// Maximum length of Excel status string.
#define EXITS(s) {sprintf(status,"%.*s",STATUS_MAX,(s)); return;}
#define nrerror		throw


typedef enum
{
	Yes = 1,
	No,

	Put = 0,
	Call = 1,

	UpOut = 1,
	UpIn,
	DnOut,
	DnIn,

	NonQuanto = 1,
	Quanto,
	Compo,

	// Barrier
	NotTouched = 0,
	Touched,

	// Rebate
	PayAtTouch = 1,
	PayAtExpiry,

	LookBack = 1,
	HindSight,
	CallMinPutMax,

	L_FactorTree = 1,
	PureImpliedTree,

	FwdSkew = 1,
	FwdATMVol,

	Term_L = 1,
	Fwd_Fwd_L,

	// Option settlement.
	CashSettle = 1,
	PhysicalSettle,

	// Option type.
	European = 1,
	American,
	Bermudan,
	CallPut_able,
	Euro_JumpD,        // jump diffusion.

	// For callable bonds.
	Callable = 1,
	Putable,

	// Volatility dynamics.
	StaticVol = 1,      // Static Vol Surface.
	StickyDelta,        // Vol Surface Moves With Spot -- Convolution.
	StickySkew,         // ATM-Vol From The Static Vol Surface, Skew Moves With Spot -- Sophis.

	// Model methods used by the PDE solver.
	DiscreteDiv = 1,    // Handle discrete dividends rigorously.
	ProportionalDiv,    // Turn discrete dividends into proportional dividends.
	ContinuousYield,    // Merge discrete dividends with repo and treat as continuous yield.

	// Boundary searching for PDE.
	UpLimit = 1,
	LowLimit,

	// PDE attached tail model.
	None = 1,		//	No tail.
	AsianSimple = 2,		//	Use simple Asain - continuous moment matching.
	AsianCG = 3,		//	Use conditioning geometric Asian model.
	LookbackFlatVol = 4,		//	Using flat vol
	LookbackCB1Vol = 5,		//	Using cascading barrier and one vol  for barrier
	LookbackCB2Vol = 6,		//	Using cascading barrier and two vols for barrier

	//For Barrier Vol model
	VolAtStrike = 1,
	VolAtBarrier = 2,
	AsRainbowBarrier = 3,	//	Use rainbow barrier option
	FromRBProb = 4,	//	Using rainbow barrier option
	FromRBProbCorr = 5,	//	Using rainbow barrier option

	//For Lookback Vol model
	FlatVol = 1,
	WithSkew1 = 2,	//	Cascading Barrier with simple barrier option
	WithSkew2 = 3,	//	Cascading Barrier with rainbow barrier option

	// Which vol functional form: Tanh or PolyBlend.
	VolFuncTanh = 55,
	VolFuncPoly = 33,
}Enums;

// Accrual basis.
typedef enum {
	ANFPED_basis_30E_360 = 0,
	ANFPED_basis_30_360 = 1,
	ANFPED_basis_actual_360 = 2,
	ANFPED_basis_actual_365 = 3,
	ANFPED_basis_actual_actual = 4
} ANFPED_basis;

const double pi = 3.141592653589790;
const double days_pa = 365.;

//inline�������������Ƶ�����ú�����ջ�ռ���������⡣
//�������Դ������͸���Ϊ���ۣ�����ʡȥ�˺������õĿ������Ӷ���ߺ�����ִ��Ч�ʡ����ִ�к������ڴ����ʱ�䣬����ں������õĿ����ϴ���ôЧ�ʵ��ջ����١�
//�������ڲ��ܳ���ѭ��

//ģ�庯�����巽ʽ��template <class T> T f(T a,T b){}
template <class T>T inline maxi(T a, T b) { return a > b ? a : b; }
template <class T>T inline mini(T a, T b) { return a < b ? a : b; }
template<class T>inline int signof(T& a) { if (a < 0) return(-1); return(1); }

//C++������洢�ļ��ַ�ʽ
//�������飺int a[40];��С�̶����ٶȿ�
//����ָ�룺 int *a; a = new int[40]; �����a���Եȼ����������飻����a����ʾ�����׸�Ԫ�صĵ�ַ��
//vector: #include<vector> vector<int> a(N); a.push_back(888); ʹ��new/delete�����ڴ���� ���ȿɱ䣬��Ч�ʵ�
//array: #include<array> array<int,2> a = {1,2}; ���ȹ̶�������ȫ�ӿڣ�Ч��ͬ�������顣
//valarray: #include<valarray> valarray<int> a(10); ������ֵ��������飻 ��ͨ��resize�ı��С��ע�������С�����ݽ���ʧ��shift/cshift���������ݵ�ƽ��
typedef valarray< double > Array1D;
typedef valarray< Array1D > Array2D;
typedef valarray< Array2D > Array3D;
typedef valarray< long > Array1L;
typedef valarray< Array1L > Array2L;



//2. ���峣�õ����ݽṹ

// ************************************************************************************
// Memory allocation & other common functions.
// ************************************************************************************


//2.1 �ɱ䳤��һά���ݽṹ
//tvec---�ɱ䳤�ȵ�����;ģ��
//dvec---�̳���tvec<double>
//ivec---�̳���tvec<int>
//cvec---�̳���tvec<char> ��append���ַ���              //������util.h��δ���壬��utils.h�ж��壻��excel.h���õ�cvec
//cvecv---�̳���tvec<cvec>  �ɱ䳤��һϵ���ַ������ɵ�����

//2.2 �ɱ䳤�Ķ�ά���ݽṹ
//imat---�ɱ�ߴ�Ķ�ά���飬���int����ֵ
//dmat---�ɱ�ߴ�Ķ�ά���飬���double����ֵ
//dmatv---Ӧ����һϵ��dmat���ɵ����飬����Ϊʲô�ö�άָ��

//�Զ���ɱ䳤����ĳ��ýӿ�
//1. ���캯������һ������n1��ʼ���������ٳ���Ϊn1�Ŀռ�
//2. ���캯������ͨ��Array1D,Array2D�ȹ���һά�ɱ䳤����
//3. resize(n0): �����ɱ䳤���鳤��Ϊn0��n0<n1���������(n1-n0)��ֵ;n0>n1����
//4. ��ʽת��operator: tvec������������ڲ�����p,ͨ��[],*��ȡp�ڲ���ֵ
//5. ���ƹ��캯����dvec a(b); ����b��һ���Ѿ����ڵ�dvec����
//6. operator=���أ�dvec a = b;



// Variable size vector.  Behaves as Type*.
template<class Type> class tvec {
	public:
		//ֻ��Ψһ�Ĵ��ι��캯��
		tvec(int n1 = 0) : n(n1), p(new Type[n1]) {};
		//���캯����ʼ���б���ʼ���б�ĳ�Ա��ʼ��˳��:C++ ��ʼ�����Աʱ���ǰ���������˳���ʼ���ģ������ǰ��ճ����ڳ�ʼ���б��е�˳��
		~tvec(void) { delete[] p; }
		void resize(int n0) {
			//������valarray�е�resize: 
			//valarray��resize���賤�ȣ��������³�ʼ����myarray.resize(3,1);����Ϊ3��1.û��1�Ļ���Ĭ��3��0. 
			//�����resize����ԭ�����е�ֵ������ȥ�����������㡣
			int		i;
			Type* p0 = new Type[n0];
			for (i = 0; i < (n < n0 ? n : n0); i++) p0[i] = (*this)[i];//������λ�õ�Ԫ�ػ������ڱ���λ��
			for (i = (n < n0 ? n : n0); i < n0; i++) p0[i] = 0;//������ȱ䳤����0��
			delete[] p;
			p = p0;
			n = n0;
		}
		//operator��ʹ�÷���֮һ����ʽ����ת����operator casting; ��һ�ֳ����ô�Ϊ�����������ؼ�operator overloading
		//����Ҫ�ĵط�tvec������Ե�������ı���pʹ�ã���tvec a; a[0]����ȡ��p��ָ��������е�ֵ��
		inline operator Type* (void) const { return p; }	// Violate const.
		int		n;//����ĳ���
		Type*   p;//���鳤�ȵ��׵�ַ

		//���ƹ��캯��
		tvec(const tvec& v) : n(v.n), p(v.n > 0 ? new Type[v.n] : 0)
		{
			int	i;
			for (i = 0; i < v.n; i++) p[i] = v.p[i];
		}

		//������=������
		tvec& operator=(const tvec& v)
		{
			int		i;
			resize(v.n);
			for (i = 0; i < v.n; i++) p[i] = v.p[i];
			return *this;
		}

	private:
};

// Variable size vector.  Behaves as double*.
class dvec : public tvec<double> {
public:
	//1.���಻�̳и���Ĺ��캯����
	//2.������๹�캯�����Σ����������ʽ���ø���Ĺ��캯�������÷�����ʾΪ":father(n1)". �����������캯����ʵ���˸��๹�캯���ĵ��á�
	dvec(int n1 = 0) : tvec<double>(n1) {};//A(b)ʵ��������������b��ʼ��A�ĳ�Աx;��x=b; �ǳ�ʼ���б�ʽ��
	dvec(const Array1D& a);//��Array1D����dvec
	dvec(const Array2D& a);//��Array2D����dvec
	double sum();//���
private:
	dvec(double);//���޶��壡��
};

// Variable size vector.  Behaves as int*.
class ivec : public tvec<int> {
public:
	ivec(int n1 = 0) : tvec<int>(n1) {};
	ivec(const Array1D& a);
	ivec(const Array2D& a);//��Array2D����ivec
private:
	ivec(double);//���޶���
};

//Variable size Vector. ������ char*
//C++�洢�ַ����ļ��ַ�ʽ
// char str0[2] = {'a','b','\0'};
// char str1[] = "ab\0";//��������Ϊ���飬��ü���\0
// char *str2 = "ab"; //char * �� char []������char*�Ǳ�����ֵ���Ըı䣬 char[]�ǳ�����ֵ���ܸı�

//strlen, strcpy, strncpy������string.h�С�
//string.h, cstring, string �����������
//string.h�ļ���C���ԵĿ⣻cstring��C++��string.h�ļ����������װ;
//string����ǰ�������б��ʲ���ˡ�����C++�Լ�������װ���࣬ͬ�������ַ��������������õ��˺ܶ�Ĳ��������صȷ�����ʵ�ַ�����C��׼���е�string.h�кܴ���


//Variable size Vector. Behaves as char*. 
//���Ե���append����չԭ�ַ���
class cvec : public tvec<char> {
public:
	cvec(int n1 = 0) : tvec<char>(n1) {};
	//strlen,strcpy,strncpy������������string.h�У�����δ����string.hҲ����ʹ�ã��������������!?
	//strlen:�ó��Ľ������������\0;����Ҫ����1��
	//strcpy:�Ѻ���'\0'���������ַ������Ƶ���һ����ַ�ռ䣬����ֵ������Ϊchar*;ǰ����ΪĿ���ַ�������Ϊ������Դ��ַ��
	//strncpy: ����Դ�ַ���ǰn���ַ����Ƶ���һ����ַ�ռ䣬����ֵ������Ϊchar*;ǰ����ΪĿ���ַ�������Ϊ������Դ��ַ��
	cvec(const char* s) : tvec<char>(s ? static_cast<int>(strlen(s) + 1) : 0) { if (s) strcpy(p, s); };
	cvec(const char* s, int n) : tvec<char>(n + 1) { strncpy(p, s, n); p[n] = 0; };
	void append(const char* s); //�ú����Ǵ�utils.cpp�и��Ƶ�util.cpp��
};


// Vector of cvec:
class cvecv : public tvec<cvec> {
public:
	cvecv(int n2 = 0) : tvec<cvec>(n2) {};//���������ø���ĺ��ι��캯��
};
// ************************************************************************************

// Variable size matrix of integers
class imat {
public:
	//imat(int r1, int c1) : r(r1), c(c1), p(new int[r1 * c1]) {}
	//����������־���using operator"*" on a 4 byte value and then casting the value to a 8 byte value.
	//ԭ��int[]�з���ı�����size_t���ͣ�ռ��8 �ֽڣ� ������int�������int���ͣ�ռ��4���ֽڣ���4�ֽڱ����ŵ���Ӧ�÷�8�ֽڵ�λ�þͳ����˾��档

	imat(int r1, int c1) : r(r1), c(c1), p(new int[static_cast<long>(r1) * static_cast<long>(c1)]) {}//4�ֽ�ǿ��ת��Ϊ8�ֽ�
	//32λ��������int/float/unsigned int/long/unsigned long: 4�ֽڣ� long long/double��8�ֽ�
	//64λ������: int/float/unsigned int:4�ֽڣ� long/long long/unsigned long/double:8�ֽ�
	~imat(void) { delete[] p; }
	void resize(int r, int c);
	//const�޶��� ֻ�ܶ�ȡ�������޸ĳ�Ա����
	inline int* operator[](int r_) const { return p+ static_cast<long>(r_*c); }
	inline operator int* (void) const { return p; }	

	//��Ա�������У��У���ʼָ��
	int	r;
	int	c;	
	int* p;
private:
	//��������������Ϊ�Լ����壬֮ǰֻ��������δ���塣
	imat(const imat& d);//���ƹ���
	imat& operator=(const imat&);//operator=����������
};





// ************************************************************************************

// Variable size matrix
//����һ���ɱ�size�Ķ�άdouble����
class dmat {
public:

	dmat(int r1, int c1) : r(r1), c(c1), p(new double[static_cast<long>(r1*c1)]) {} //�ó��ȣ���ȹ���һ�����󡣳�Ա���������ð�ű�ʾ��ֵ
	dmat(const Array2D& a); //���Զ����valarray<valarray<double>>�����͹���dmat
	~dmat(void) { delete[] p; }//�ͷſռ�
	void resize(int r, int c); //�ɷŴ󣺷Ŵ�ʱ���㣻����С����Сʱά��ԭλ�ò���
	void transpose(void);//ת��
	inline double* operator[](int r_) const { return p + static_cast<long>(r_*c); }
	//����һ��ָ��double��ָ�룬ָ��Ӽ�����ʾ��������ƽ����ͬ������λ�ã�
	//(*this)[i][j] �״�ʹ��ƽ�Ƶ���i�п�ͷ�ĵ�һ��Ԫ��-->��*this��[i] ���ص��Ǹ�dmat�����е�����ָ��-->�ټ�[j]��ʾ����ָ�����ƽ��j��λ��
	//this[i][j]���Ӧ���Ƿ���ָ��ĳ���ض�Ԫ�ص�ָ�룬ȡֵ��Ӧ����*(this[i][j]),Ӧ����(*this)[i][j]�ȼ۰ɣ�������һ����ʽת������

	inline operator double* (void) const { return p; }	// Violate const.��ʽת�����������ת��Ϊ�ڲ���ָ�����p,�������ϵ�(*this)���ص���p
	int			r;
	int			c;	// Can't change this once data entered.
	double* p;
private:
	dmat(const dmat& d);
	dmat& operator=(const dmat&);
};


/* ��ED_Analytics.h�У�������MemoryAllocator.cpp�еĿ��ٶ�ά����ĳ��� */
//һ��û��Ҫ�ã�ֱ����dmat����󼴿ɡ�
double** New2DMatrix(long nrow, long ncol);
void Delete2DMatrix(double** m);



//***********************************************************
// Vector of matrices.
//��ʾ��һ��dmat�����������
class dmatV {
public:
	dmatV(int n, const Array2D** mat = 0);
	~dmatV(void);
	inline operator const dmat** (void) const { return p; }
	int	n;//
private:
	const dmat** p;//dmat���ɵ����飬����ֻҪһ��*�Ϳ�������
	//����dmat��������һ����ά����
	//dmat** p��ʾ����һ��ָ��һϵ��dmat�����ַ����ָ�룻pָ���ָ����ָ���dmat�������ڴ��в�һ������
	//dmat*p �����p��ʾdmat������׵�ַ����ô��Щdmat�������ڴ����������ġ�
private:
	//dmatV(const dmatV&);
	//dmatV operator=(const dmatV&);
};

//***********************************************************

// N<=3 dimensional rectangle of numbers.
class dcube {
public:
	dcube(int a_ = 1, int b_ = 1, int c_ = 1) : a(a_), b(b_), c(c_),
		p(new double[a * b * c]) {
		if (!p) throw " Out of memory";
	}
	~dcube(void) { delete[] p; }
	inline double& operator()(int i = 0, int j = 0, int k = 0) {
		return p[i * b * c + j * c + k];
	}
	inline double operator()(int i = 0, int j = 0, int k = 0) const {
		return p[i * b * c + j * c + k];
	}
	inline double* ptr(int i = 0, int j = 0, int k = 0) const {	// Violate const.
		return p + i * b * c + j * c + k;
	}
	inline double* ptr(int i = 0, int j = 0, int k = 0) {
		return p + i * b * c + j * c + k;
	}
	int		a, b, c;
	double* p;
	//inline friend void swap(dcube& a, dcube& b) { double* p = a.p; a.p = b.p; b.p = p; }
private:
	//dcube(const dcube&);
	//dcube& operator=(const dcube&);
};