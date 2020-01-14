#ifndef HAVE_AD_BLAS_HPP
#define HAVE_AD_BLAS_HPP
// Autogenerated - do not edit by hand !
#define GLOBAL_HASH_TYPE unsigned int
#define GLOBAL_COMPRESS_TOL 16
#define GLOBAL_UNION_OR_STRUCT union
#define stringify(s) #s
#define xstringify(s) stringify(s)
#define THREAD_NUM 0
#define GLOBAL_INDEX_VECTOR std::vector<GLOBAL_INDEX_TYPE>
#define GLOBAL_INDEX_TYPE unsigned int
#define ASSERT2(x, msg)                          \
  if (!(x)) {                                    \
    Rcerr << "ASSERTION FAILED: " << #x << "\n"; \
    Rcerr << "POSSIBLE REASON: " << msg << "\n"; \
    abort();                                     \
  }
#define GLOBAL_MAX_NUM_THREADS 48
#define INDEX_OVERFLOW(x) \
  ((size_t)(x) >= (size_t)std::numeric_limits<GLOBAL_INDEX_TYPE>::max())
#define ASSERT(x)                                \
  if (!(x)) {                                    \
    Rcerr << "ASSERTION FAILED: " << #x << "\n"; \
    abort();                                     \
  }
#define GLOBAL_REPLAY_TYPE ad_aug
#define GLOBAL_MIN_PERIOD_REP 10
#define INHERIT_CTOR(A, B)                                       \
  A() {}                                                         \
  template <class T1>                                            \
  A(const T1 &x1) : B(x1) {}                                     \
  template <class T1, class T2>                                  \
  A(const T1 &x1, const T2 &x2) : B(x1, x2) {}                   \
  template <class T1, class T2, class T3>                        \
  A(const T1 &x1, const T2 &x2, const T3 &x3) : B(x1, x2, x3) {} \
  template <class T1, class T2, class T3, class T4>              \
  A(const T1 &x1, const T2 &x2, const T3 &x3, const T4 &x4)      \
      : B(x1, x2, x3, x4) {}
#define GLOBAL_SCALAR_TYPE double
#include <Eigen/Dense>
#include "global.hpp"

namespace TMBad {

/** \brief Request a contiguous block on the tape.

    1. Check if `x` already is on the tape and satisfies the storage
   requirement.
    2. If **no** invoke a deep copy of `x` to the tape **and** *update* `x` with
   the new tape addresses.

    \return A reference to `x` as a contiguous block on the tape.

    \note The update step is critical as it ensures that a given
    matrix can be used several times without invoking a deep copy more
    than once.
*/
template <class Matrix>
global::ad_range contiguousBlock(const Matrix &x) {
  bool yes = true;
  Index j_previous = -1;
  for (size_t i = 0; i < (size_t)x.size(); i++) {
    if (!x(i).ontape()) {
      yes = false;
      break;
    }
    Index j = ad_plain(x(i)).index;
    if (i > 0) {
      if (j != j_previous + 1) {
        yes = false;
        break;
      }
    }
    j_previous = j;
  }
  if (yes) {
    return global::ad_range(ad_plain(x(0)), x.rows(), x.cols());
  }

  ad_plain ans;
  for (size_t i = 0; i < (size_t)x.size(); i++) {
    ad_plain xi_cpy = x(i).copy();

    x(i).override_by(xi_cpy);
    if (i == 0) ans = xi_cpy;
  }
  return global::ad_range(ans, x.rows(), x.cols());
}

using Eigen::Dynamic;
using Eigen::Map;
using Eigen::Matrix;
typedef Matrix<double, Dynamic, Dynamic> dmatrix;
typedef Matrix<global::Replay, Dynamic, Dynamic> vmatrix;

template <class Target>
void fill(Target &y, const global::ad_range x) {
  ASSERT((size_t)y.size() == (size_t)x.size());
  ad_plain xx = x;
  for (size_t i = 0; i < (size_t)y.size(); i++) {
    y(i) = xx;
    xx.index++;
  }
}

template <bool XT, bool YT, bool ZT>
struct MatMul;
template <bool XT, bool YT, bool ZT>
void matmul(const vmatrix &x, const vmatrix &y, Map<vmatrix> z) {
  global::ad_range xc = contiguousBlock(x);
  global::ad_range yc = contiguousBlock(y);
  global::ad_range out = get_glob()->add_to_stack<MatMul<XT, YT, ZT> >(xc, yc);
  fill(z, out);
}

/** \brief Multiply two matrices of ad variables */
vmatrix matmul(const vmatrix &x, const vmatrix &y);

/** \brief Multiply two matrices of scalar types */
dmatrix matmul(const dmatrix &x, const dmatrix &y);

/** Expand all 8 combinations */
template <bool XT, bool YT, bool ZT>
void matmul(Map<const dmatrix> x, Map<const dmatrix> y, Map<dmatrix> z) {
  if (XT && YT && ZT) z.transpose() = x.transpose() * y.transpose();
  if (!XT && YT && ZT) z.transpose() = x * y.transpose();
  if (XT && !YT && ZT) z.transpose() = x.transpose() * y;
  if (XT && YT && !ZT) z = x.transpose() * y.transpose();
  if (!XT && !YT && ZT) z.transpose() = x * y;
  if (XT && !YT && !ZT) z = x.transpose() * y;
  if (!XT && YT && !ZT) z = x * y.transpose();
  if (!XT && !YT && !ZT) z = x * y;
}

template <bool XT, bool YT, bool ZT>
struct MatMul : global::Operator<2, -1> {
  static const bool dynamic = true;
  static const int max_fuse_depth = 0;
  int n1, n2, n3;
  static const int ninput = 2;
  MatMul(global::ad_range X, global::ad_range Y) {
    set_dim(X.rows(), X.cols(), Y.rows(), Y.cols());
  }
  MatMul(int n1, int n2, int n3) : n1(n1), n2(n2), n3(n3) {}
  Index input_size() const { return 2; }
  Index output_size() const {
    int Xrows, Xcols, Yrows, Ycols, Zrows, Zcols;
    get_dim(Xrows, Xcols, Yrows, Ycols, Zrows, Zcols);
    return Zrows * Zcols;
  }
  static const bool have_input_size_output_size = true;
  void set_dim(int Xrows, int Xcols, int Yrows, int Ycols) {
    n1 = Xrows;
    n2 = Xcols;
    n3 = (YT ? Yrows : Ycols);
  }
  void get_dim(int &Xrows, int &Xcols, int &Yrows, int &Ycols, int &Zrows,
               int &Zcols) const {
    Xrows = n1;
    Xcols = n2;

    int Xop_rows = Xrows, Xop_cols = Xcols;
    if (XT) std::swap(Xop_rows, Xop_cols);

    int Yop_rows = Xop_cols, Yop_cols = n3;

    Yrows = Yop_rows;
    Ycols = Yop_cols;
    if (YT) std::swap(Yrows, Ycols);

    int Zop_rows = Xop_rows, Zop_cols = Yop_cols;

    Zrows = Zop_rows;
    Zcols = Zop_cols;
    if (ZT) std::swap(Zrows, Zcols);
  }
  template <class Type>
  void forward(ForwardArgs<Type> &args) {
    int Xrows, Xcols, Yrows, Ycols, Zrows, Zcols;
    get_dim(Xrows, Xcols, Yrows, Ycols, Zrows, Zcols);
    typedef Map<Matrix<Type, Dynamic, Dynamic> > MapMatrix;
    typedef Map<const Matrix<Type, Dynamic, Dynamic> > ConstMapMatrix;
    ConstMapMatrix X(args.x_ptr(0), Xrows, Xcols);
    ConstMapMatrix Y(args.x_ptr(1), Yrows, Ycols);
    MapMatrix Z(args.y_ptr(0), Zrows, Zcols);
    matmul<XT, YT, ZT>(X, Y, Z);
  }
  template <class Type>
  void reverse(ReverseArgs<Type> &args) {
    int Xrows, Xcols, Yrows, Ycols, Zrows, Zcols;
    get_dim(Xrows, Xcols, Yrows, Ycols, Zrows, Zcols);
    typedef Map<Matrix<Type, Dynamic, Dynamic> > MapMatrix;
    typedef Map<const Matrix<Type, Dynamic, Dynamic> > ConstMapMatrix;
    ConstMapMatrix X(args.x_ptr(0), Xrows, Xcols);
    ConstMapMatrix Y(args.x_ptr(1), Yrows, Ycols);
    ConstMapMatrix W(args.dy_ptr(0), Zrows, Zcols);
    MapMatrix DX(args.dx_ptr(0), Xrows, Xcols);
    MapMatrix DY(args.dx_ptr(1), Yrows, Ycols);

    Matrix<Type, Dynamic, Dynamic> DX0_tmp(DX.rows(), DX.cols());
    Matrix<Type, Dynamic, Dynamic> DY0_tmp(DY.rows(), DY.cols());
    MapMatrix DX0(&DX0_tmp(0), DX.rows(), DX.cols());
    MapMatrix DY0(&DY0_tmp(0), DY.rows(), DY.cols());

    matmul<ZT, !YT, XT>(W, Y, DX0);
    matmul<!XT, ZT, YT>(X, W, DY0);

    DX += DX0;
    DY += DY0;
  }

  void dependencies(Args<> &args, Dependencies &dep) const {
    int Xrows, Xcols, Yrows, Ycols, Zrows, Zcols;
    get_dim(Xrows, Xcols, Yrows, Ycols, Zrows, Zcols);
    dep.add_segment(args.input(0), Xrows * Xcols);
    dep.add_segment(args.input(1), Yrows * Ycols);
  }
  static const bool have_dependencies = true;
  /** \brief This operator **has** implicit dependencies */
  static const bool implicit_dependencies = true;
  /** \brief It is **not* safe to remap the inputs of this operator */
  static const bool allow_remap = false;

  void forward(ForwardArgs<Writer> &args) { ASSERT(false); }
  void reverse(ReverseArgs<Writer> &args) { ASSERT(false); }
  const char *op_name() { return "MatMul"; }
};

}  // namespace TMBad
#endif  // HAVE_AD_BLAS_HPP
