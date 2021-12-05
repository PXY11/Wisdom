#pragma once//同一个文件不会被包含多次

#pragma warning(disable:4018)// 有符号/无符号不匹配 警告. 如unsigned int 类型 与 int 类型变量之间的比较。
#pragma warning(disable:4244)// 数据类型转型可能数据丢失。 如讲double转化为int. 可使用static_cast<int>进行强制类型转换。

//该base文件参考了以下文件
//SharedDefines.h---定义了Array1D等


#define NOMINMAX
#include <Windows.h>
//#define NOMINMAX告诉编译器(或实际上是预处理器)跳过min和max的定义,但只有在#include "windows.h"之前执行此操作时才会应用它。

#include <algorithm>//C++标准库：包含了所有vector、list、set、map操作能想到的一些函数,如查找,替换、排序、计数等常用的功能全部在里面
#include <fstream>//C++中的IO库类型
#include <strstream>//C++中C语言风格的输入输出
#include <valarray>
#include <float.h> //C语言标准库头文件，定义了浮点型的范围，精度。用于数值分析。
using namespace std;//valarray是C++ STL标准库；可以添加此语句，也可以使用std::valarray

//1. 使用#define定义常量和常用数学函数

//#define  标识符  常量   //注意, 最后没有分号
//凡是以“#”开头的均为预处理指令,预编译不是编译，而是编译前的处理。这个操作是在正式编译之前由系统自动完成的。预编译所执行的操作就是简单的“文本”替换。
// 
//#define与typedef的区别
//#define用标识符表示一个字符串，成为宏； typedef用于给一个变量起个新的名字，便于阅读，发生在编译阶段。

//#define与const定义常量的区别
//#define在预处理阶段替换；const在编译时确定值
//#define无类型，不检查类型安全性；const有数据类型，编译时检查
//#define多次用到多次拷贝；const在静态存储区分配内存，程序运行过程中内存中只有一个拷贝


#define MAX(a,b)	( (a)<(b) ? (b) : (a) )
#define MIN(a,b)	( (a)<(b) ? (a) : (b) )
#define MIN4(a,b,c,d)	( MIN( MIN((a),(b)) , MIN((c),(d)) ) ) //rdv.cpp需使用
#define ABS(a)		( (a)< 0  ?-(a) : (a) )
#define SGN(a)		( (a)< 0  ? -1  :  1  )
#define SWAP(a,b)   {temp=(a); (a)=(b); (b)=temp;} //原先定义在NumUtilities.cpp中
#define EPS         (1.e-8 )
#define EPS1		( 1e-8 )
#define MAX_STR	    (1<<12)		// Maximum string length. 2**12=4096. 表示将二进制1向左移动的位数
#define JMAX        200//原先定义在NumUtilities.cpp的头部
#define	ROOT_2PI	( 2.506628274631 )//原先定义在util.h; 在CumNorm.h中同样定义该值，#define DBA_SQRT2PI
#define MAXN		(1000)  //以下两个在lu.cpp中定义的
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

//inline内联函数，解决频繁调用函数对栈空间的消耗问题。
//内联是以代码膨胀复制为代价，仅仅省去了函数调用的开销，从而提高函数的执行效率。如果执行函数体内代码的时间，相比于函数调用的开销较大，那么效率的收获会很少。
//函数体内不能出现循环

//模板函数定义方式：template <class T> T f(T a,T b){}
template <class T>T inline maxi(T a, T b) { return a > b ? a : b; }
template <class T>T inline mini(T a, T b) { return a < b ? a : b; }
template<class T>inline int signof(T& a) { if (a < 0) return(-1); return(1); }

//C++中数组存储的几种方式
//内置数组：int a[40];大小固定，速度快
//数组指针： int *a; a = new int[40]; 这里的a可以等价于内置数组；两个a均表示数组首个元素的地址。
//vector: #include<vector> vector<int> a(N); a.push_back(888); 使用new/delete进行内存管理； 长度可变，灵活，效率低
//array: #include<array> array<int,2> a = {1,2}; 长度固定，更安全接口，效率同内置数组。
//valarray: #include<valarray> valarray<int> a(10); 面向数值计算的数组； 可通过resize改变大小，注意调整大小后，数据将丢失。shift/cshift可用于数据的平移
typedef valarray< double > Array1D;
typedef valarray< Array1D > Array2D;
typedef valarray< Array2D > Array3D;
typedef valarray< long > Array1L;
typedef valarray< Array1L > Array2L;



//2. 定义常用的数据结构

// ************************************************************************************
// Memory allocation & other common functions.
// ************************************************************************************


//2.1 可变长的一维数据结构
//tvec---可变长度的数组;模板
//dvec---继承于tvec<double>
//ivec---继承于tvec<int>
//cvec---继承于tvec<char> 可append的字符串              //本身在util.h中未定义，在utils.h有定义；在excel.h中用到cvec
//cvecv---继承于tvec<cvec>  可变长的一系列字符串构成的数组

//2.2 可变长的二维数据结构
//imat---可变尺寸的二维数组，存放int型数值
//dmat---可变尺寸的二维数组，存放double型数值
//dmatv---应该书一系列dmat构成的数组，但是为什么用二维指针

//自定义可变长数组的常用接口
//1. 构造函数：由一个整数n1初始化，并开辟长度为n1的空间
//2. 构造函数：可通过Array1D,Array2D等构造一维可变长数组
//3. resize(n0): 调整可变长数组长度为n0，n0<n1则丢弃后面的(n1-n0)个值;n0>n1则补零
//4. 隐式转换operator: tvec对象可类似于内部对象p,通过[],*来取p内部的值
//5. 复制构造函数：dvec a(b); 其中b是一个已经存在的dvec对象
//6. operator=重载：dvec a = b;



// Variable size vector.  Behaves as Type*.
template<class Type> class tvec {
	public:
		//只有唯一的带参构造函数
		tvec(int n1 = 0) : n(n1), p(new Type[n1]) {};
		//构造函数初始化列表；初始化列表的成员初始化顺序:C++ 初始化类成员时，是按照声明的顺序初始化的，而不是按照出现在初始化列表中的顺序。
		~tvec(void) { delete[] p; }
		void resize(int n0) {
			//区别于valarray中的resize: 
			//valarray的resize重设长度，并且重新初始化。myarray.resize(3,1);重设为3个1.没有1的话，默认3个0. 
			//这里的resize保留原数组中的值，多了去掉；不够补零。
			int		i;
			Type* p0 = new Type[n0];
			for (i = 0; i < (n < n0 ? n : n0); i++) p0[i] = (*this)[i];//将本来位置的元素还放置于本来位置
			for (i = (n < n0 ? n : n0); i < n0; i++) p0[i] = 0;//如果长度变长，补0；
			delete[] p;
			p = p0;
			n = n0;
		}
		//operator的使用方法之一：隐式类型转换符operator casting; 另一种常用用处为：操作符重载即operator overloading
		//在需要的地方tvec对象可以当作里面的变量p使用，即tvec a; a[0]就是取得p所指向的数组中的值。
		inline operator Type* (void) const { return p; }	// Violate const.
		int		n;//数组的长度
		Type*   p;//数组长度的首地址

		//复制构造函数
		tvec(const tvec& v) : n(v.n), p(v.n > 0 ? new Type[v.n] : 0)
		{
			int	i;
			for (i = 0; i < v.n; i++) p[i] = v.p[i];
		}

		//操作符=的重载
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
	//1.子类不继承父类的构造函数；
	//2.如果父类构造函数含参，子类必须显式调用父类的构造函数；调用方法表示为":father(n1)". 以下三个构造函数均实现了父类构造函数的调用。
	dvec(int n1 = 0) : tvec<double>(n1) {};//A(b)实际上是做的是用b初始化A的成员x;即x=b; 是初始化列表方式。
	dvec(const Array1D& a);//用Array1D构造dvec
	dvec(const Array2D& a);//用Array2D构造dvec
	double sum();//求和
private:
	dvec(double);//暂无定义！？
};

// Variable size vector.  Behaves as int*.
class ivec : public tvec<int> {
public:
	ivec(int n1 = 0) : tvec<int>(n1) {};
	ivec(const Array1D& a);
	ivec(const Array2D& a);//用Array2D构造ivec
private:
	ivec(double);//暂无定义
};

//Variable size Vector. 类似于 char*
//C++存储字符串的几种方式
// char str0[2] = {'a','b','\0'};
// char str1[] = "ab\0";//以上两种为数组，最好加上\0
// char *str2 = "ab"; //char * 与 char []的区别：char*是变量，值可以改变， char[]是常量，值不能改变

//strlen, strcpy, strncpy定义在string.h中。
//string.h, cstring, string 三个库的区别：
//string.h文件是C语言的库；cstring是C++对string.h的简略升级与包装;
//string就与前面两个有本质差别了。它是C++自己开发封装的类，同样用于字符串操作，其中用到了很多的操作符重载等方法，实现方法和C标准库中的string.h有很大差别。


//Variable size Vector. Behaves as char*. 
//可以调用append来扩展原字符串
class cvec : public tvec<char> {
public:
	cvec(int n1 = 0) : tvec<char>(n1) {};
	//strlen,strcpy,strncpy三个均定义在string.h中，这里未引入string.h也可以使用，是在哪里引入的!?
	//strlen:得出的结果不包括最后的\0;所以要加上1；
	//strcpy:把含有'\0'结束符的字符串复制到另一个地址空间，返回值的类型为char*;前参数为目标地址，后参数为数据来源地址。
	//strncpy: 把来源字符串前n个字符复制到另一个地址空间，返回值的类型为char*;前参数为目标地址，后参数为数据来源地址。
	cvec(const char* s) : tvec<char>(s ? static_cast<int>(strlen(s) + 1) : 0) { if (s) strcpy(p, s); };
	cvec(const char* s, int n) : tvec<char>(n + 1) { strncpy(p, s, n); p[n] = 0; };
	void append(const char* s); //该函数是从utils.cpp中复制到util.cpp中
};


// Vector of cvec:
class cvecv : public tvec<cvec> {
public:
	cvecv(int n2 = 0) : tvec<cvec>(n2) {};//子类必须调用父类的含参构造函数
};
// ************************************************************************************

// Variable size matrix of integers
class imat {
public:
	//imat(int r1, int c1) : r(r1), c(c1), p(new int[r1 * c1]) {}
	//上面的语句出现警告using operator"*" on a 4 byte value and then casting the value to a 8 byte value.
	//原因int[]中放入的变量是size_t类型，占用8 字节； 而两个int相乘仍是int类型，占用4个字节；将4字节变量放到本应该放8字节的位置就出现了警告。

	imat(int r1, int c1) : r(r1), c(c1), p(new int[static_cast<long>(r1) * static_cast<long>(c1)]) {}//4字节强制转换为8字节
	//32位编译器：int/float/unsigned int/long/unsigned long: 4字节； long long/double：8字节
	//64位编译器: int/float/unsigned int:4字节； long/long long/unsigned long/double:8字节
	~imat(void) { delete[] p; }
	void resize(int r, int c);
	//const限定了 只能读取，不能修改成员变量
	inline int* operator[](int r_) const { return p+ static_cast<long>(r_*c); }
	inline operator int* (void) const { return p; }	

	//成员变量：行，列，起始指针
	int	r;
	int	c;	
	int* p;
private:
	//以下两个函数均为自己定义，之前只是声明，未定义。
	imat(const imat& d);//复制构造
	imat& operator=(const imat&);//operator=操作符重载
};





// ************************************************************************************

// Variable size matrix
//创建一个可变size的二维double矩阵
class dmat {
public:

	dmat(int r1, int c1) : r(r1), c(c1), p(new double[static_cast<long>(r1*c1)]) {} //用长度，宽度构造一个矩阵。成员函数后面加冒号表示赋值
	dmat(const Array2D& a); //用自定义个valarray<valarray<double>>的类型构造dmat
	~dmat(void) { delete[] p; }//释放空间
	void resize(int r, int c); //可放大：放大时补零；可缩小：缩小时维持原位置不变
	void transpose(void);//转置
	inline double* operator[](int r_) const { return p + static_cast<long>(r_*c); }
	//返回一个指向double的指针，指针加减法表示在数组中平移相同整数的位置；
	//(*this)[i][j] 首次使用平移到第i行开头的第一个元素-->（*this）[i] 返回的是该dmat对象中的数组指针-->再加[j]表示数组指针向后平移j个位置
	//this[i][j]这个应该是返回指向某个特定元素的指针，取值不应该是*(this[i][j]),应该与(*this)[i][j]等价吧！？见下一个隐式转化函数

	inline operator double* (void) const { return p; }	// Violate const.隐式转换，将类对象转化为内部的指针变量p,所以以上的(*this)返回的是p
	int			r;
	int			c;	// Can't change this once data entered.
	double* p;
private:
	dmat(const dmat& d);
	dmat& operator=(const dmat&);
};


/* 在ED_Analytics.h中，定义在MemoryAllocator.cpp中的开辟二维数组的程序 */
//一般没必要用，直接用dmat类对象即可。
double** New2DMatrix(long nrow, long ncol);
void Delete2DMatrix(double** m);



//***********************************************************
// Vector of matrices.
//表示由一列dmat对象构造的数组
class dmatV {
public:
	dmatV(int n, const Array2D** mat = 0);
	~dmatV(void);
	inline operator const dmat** (void) const { return p; }
	int	n;//
private:
	const dmat** p;//dmat构成的数组，不是只要一个*就可以了吗？
	//首先dmat表面上是一个二维矩阵
	//dmat** p表示的是一个指向“一系列dmat对象地址”的指针；p指向的指针所指向的dmat对象在内存中不一定连续
	//dmat*p ；如果p表示dmat数组的首地址，那么这些dmat对象在内存中是连续的。
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