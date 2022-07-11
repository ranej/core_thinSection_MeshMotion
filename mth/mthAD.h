#ifndef MTH_AD_H
#define MTH_AD_H

#include <cmath>
#include <iostream>
#include "canArray.h"

/** \file mthAD.h
  * \brief A forward automatic differentiation variable */

namespace mth {

/** \brief forward automatic differentiation variable */
template <class T=double, unsigned int N=0>
class AD
{
  public:
    /** \brief the variable value */
    double x_;
    /** \brief the derivative array */
    T dx_[N];
    /** \brief the number of derivatives */
    enum { degree = N };
    /** \brief default constructor */
    AD():x_(0.) {zero();}
    /** \brief default constructor from a double */
    AD(double x):x_(x) {zero();}
    /** \brief copy constructor */
    AD(AD<T, N> const& other) {copy(other);}
    template <class B>
    AD(AD<B, N> const& other) {copy(other);}
    /** \brief get the size of the derivative array */
    unsigned int size() const {return N;}
    /** \brief set as the ith variable of N */
    void diff(unsigned int i, unsigned int n=0)
    {
      (void)n;
      zero();
      dx_[i] = 1.;
    }
    /** \brief get the value of the variable (mutable) */
    double& val() {return x_;}
    /** \brief get the value of the variable (immutable) */
    const double& val() const {return x_;}
    /** \brief get the ith derivative value (mutable) */
    T& dx(unsigned int i) {return dx_[i];}
    /** \brief get the ith derivative value (immutable) */
    const T& dx(unsigned int i) const {return dx_[i];}
    /** \brief resize for static AD (no-op) */
    inline void resize(unsigned int i) {(void)i;}
    /** \brief type conversion from AD to double*/
    operator double() const {return double(x_);}
    /** \brief assignment to a double */
    AD<T, N>& operator=(double other)
    {
      x_ = other;
      zero();
      return *this;
    }
    /** \brief assignment to another AD variable */
    AD<T, N>& operator=(AD<T, N> const& other)
    {
      copy(other);
      return *this;
    }

    template <class B>
    AD<T, N>& operator=(AD<B, N> const& other)
    {
      copy(other);
      return *this;
    }
    /** \brief addition assignment with a double */
    AD<T, N>& operator+=(double other)
    {
      x_ += other;
      return *this;
    }
    /** \brief addition assignment with another AD variable */
    template< class B>
    AD<T, N>& operator+=(AD<B, N> const& other)
    {
      x_ += other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] += other.dx_[i];
      return *this;
    }
    /** \brief subtraction assignment with a double */
    AD<T, N>& operator-=(double other)
    {
      x_ -= other;
      return *this;
    }
    /** \brief subtraction assignment with another AD variable */
    template< class B>
    AD<T, N>& operator-=(AD<B, N> const& other)
    {
      x_ -= other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] -= other.dx_[i];
      return *this;
    }
    /** \brief multiplication assignment with a double */
    AD<T, N>& operator*=(double other)
    {
      x_ *= other;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] *= other;
      return *this;
    }
    /** \brief multiplication assignment with another AD variable */
    template<class B>
    AD<T, N>& operator*=(AD<B, N> const& other)
    {
      x_ *= other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] = dx_[i]*other + x_*other.dx_[i];
    }
    /** \brief division assignment with a double */
    AD<T, N>& operator/=(double other)
    {
      x_ /= other;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] /= other;
      return *this;
    }
    /** \brief division assignment with another AD variable */
    template<class B>
    AD<T, N>& operator/=(AD<B, N> const& other)
    {
      x_ /= other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] = (dx_[i]*other - other*other.dx_[i]) / (other*other);
    }
  private:
    void zero()
    {
      for (unsigned int i=0; i < N; ++i)
        dx_[i] = 0.;
    }

    template <class B>
    void copy(AD<B, N> const& other)
    {
      this->x_ = other.x_;
      for (unsigned int i=0; i < N; ++i)
        this->dx_[i] = (T)other.dx_[i];
    }
};

/** \brief forward automatic differentiation variable with dynamic variable array */
template <class T>
class AD<T,0>
{
  public:
    /** \brief the variable value */
    double x_;
    /** \brief the dynamic derivative array*/
    can::Array<T> dx_;
    const static double _zero_;
    /** \brief default constructor */
    AD():x_(0), dx_() {zero();}
    /** \brief constructs from a double*/
    AD(double val):x_(val), dx_() {zero();}
    /** \brief constructs from other AD*/
    template<class B>
    AD(AD<B, 0> const& other) {zero(); copy(other);}
    unsigned size() { return dx_.size();}
    unsigned size() const { return dx_.size();}
    void diff(unsigned int i, unsigned int n)
    {
      dx_.resize_copy(n);
      zero();
      dx_[i] = 1;
    }
    /** \brief get the value of the variable (mutable) */
    double& val() {return x_;}
    /** \brief get the value of the variable (immutable) */
    const double& val() const {return x_;}
    /** \brief get the ith derivative value (mutable) */
    T& dx(unsigned int i)
    {
      if(i >= size())
      {
        unsigned int temp_size = size();
        dx_.resize_copy(i + 1);
        for(unsigned int x = temp_size; x <= i; x++)
          dx_[x] = 0.;
      }
      return dx_[i];
    }
    /** \brief get the ith deriative value (immutable) */
    const T dx(unsigned int i) const {
      if(i >= size())
        return T(_zero_);
      return dx_[i];
    }
    /** \brief type conversion from AD to double */
    operator double() const {return double(x_);}
    /** \brief assignment to a double */
    AD<T, 0>& operator=(double other)
    {
      x_ = other;
      zero();
      return *this;
    }
    /** \brief assignment to another AD variable */
    template<class B>
    AD<T, 0>& operator=(AD<B, 0> const& other)
    {
      resize(other.size());
      copy(other);
      return *this;
    }
    /** \brief addition assignment with a double */
    AD<T, 0>& operator+=(double other)
    {
      x_ += other;
      return *this;
    }
    /** \brief addition assignment with another AD variable */
    AD<T, 0>& operator+=(AD<T, 0> const& other)
    {
      resize(other.size());
      x_ += other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] += other.dx_[i];
      return *this;
    }
    /** \brief subtraction assignment with a double */
    AD<T, 0>& operator-=(double other)
    {
      x_ -= other;
      return *this;
    }
    /** \brief subtraction assignment with another AD variable */
    AD<T, 0>& operator-=(AD<T, 0> const& other)
    {
      if (other.size() > size())
        resize(other.size());
      x_ -= other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] -= other.dx_[i];
      return *this;
    }
    /** \brief multiplication assignment with a double */
    AD<T, 0>& operator*=(double other)
    {
      x_ *= other;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] *= other;
      return *this;
    }
    /** \brief multiplication assignment with another AD variable */
    AD<T, 0>& operator*=(AD<T, 0> const& other)
    {
      resize(other.size());
      x_ *= other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = dx_[i]*other.x_ + x_*other.dx_[i];
      return *this;
    }
    /** \brief division assignment with a double */
    AD<T, 0>& operator/=(double other)
    {
      x_ /= other;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] /= other;
      return *this;
    }
    /** \brief division assignment with another AD variable */
    AD<T, 0>& operator/=(AD<T, 0> const& other)
    {
      resize(other.size());
      x_ /= other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = (dx_[i]*other.x_ - x_*other.dx_[i]) / (other.x_*other.x_);
      return *this;
    }
    void resize(unsigned int n)
    {
      dx_.resize(n);
    }
  private:
    void zero()
    {
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = 0.;
    }
    template<class B> //templating allows for different derivative sizes
    void copy(AD<B, 0> const& other)
    {
      if(size() != other.size())
        dx_.resize(other.size());
      x_ = other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = (T)other.dx_[i];
    }
};

template<class T>
const double mth::AD<T, 0>::_zero_ = 0.;

/**********************
  * UNARY OPERATIONS *
***********************/

/** \brief unary subtraction */
template <class T, unsigned int N>
AD<T, N> operator-(AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  tmp  = AD<T, N>(-1.) * A;
  return tmp;
}

/***********************
  * BINARY OPERATIONS *
************************/

/** \brief binary addition between a double and an AD variable */
template <class T, unsigned int N>
AD<T, N> operator+(double L, AD<T, N> const& R)
{
  AD<T, N> tmp;
  tmp.resize(R.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(R.dx(i));
  tmp.val() = L + R.val();
  return tmp;
}

/** \brief binary addition between an AD variable and a double */
template <class T, unsigned int N>
AD<T, N> operator+(AD<T, N> const& L, double R)
{
  AD<T, N> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(L.dx(i));
  tmp.val() = L.val() + R;
  return tmp;
}
/** \brief binary addition between two AD variables */
template <class T, class B, unsigned int N>
AD<T, N> operator+(AD<T, N> const& L, AD<B, N> const& R)
{
  AD<T, N> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(L.dx(i) + R.dx(i));
  tmp.val() = L.val() + R.val();
  return tmp;
}

/** \brief binary subtraction between a double and an AD variable */
template <class T, unsigned int N>
AD<T, N> operator-(double L, AD<T, N> const& R)
{
  AD<T, N> tmp;
  tmp.resize(R.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(-R.dx(i));
  tmp.val() = L - R.val();
  return tmp;
}

/** \brief binary subtraction between an AD variable and a double */
template <class T, unsigned int N>
AD<T, N> operator-(AD<T, N> const& L, double R)
{
  AD<T, N> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(L.dx(i));
  tmp.val() = L.val() - R;
  return tmp;
}

/** \brief binary subtraction between two AD variables */
template <class T, class B, unsigned int N>
AD<T, N> operator-(AD<T, N> const& L, AD<B, N> const& R)
{
  AD<T, N> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(L.dx(i) - R.dx(i));
  tmp.val() = L.val() - R.val();
  return tmp;
}

/** \brief binary multiplication between a double and an AD variable */
template <class T, unsigned int N>
AD<T, N> operator*(double L, AD<T, N> const& R)
{
  AD<T, N> tmp;
  tmp.resize(R.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(L * R.dx(i));
  tmp.val() = L * R.val();
  return tmp;
}

/** \brief binary multiplication between an AD variable and a double */
template <class T, unsigned int N>
AD<T, N> operator*(AD<T, N> const& L, double R)
{
  AD<T, N> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(L.dx(i) * R);
  tmp.val() = L.val() * R;
  return tmp;
}

/** \brief binary multiplication between two AD variables */
template <class T, class B, unsigned int N>
AD<T, N> operator*(AD<T, N> const& L, AD<B, N> const& R)
{
  AD<T, N> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(L.dx(i) * R + L * R.dx(i));
  tmp.val() = L.val() * R.val();
  return tmp;
}

/** \brief binary division between a double and an AD variable */
template <class T, unsigned int N>
AD<T, N> operator/(double L, AD<T, N> const& R)
{
  AD<T, N> tmp;
  tmp.resize(R.size());
  T R_tmp = R; //Recursive R, used to prevent infinite recurrsion.
  for (unsigned int i = 0; i < R.size(); i++)
    tmp.dx(i) = T(-L * R.dx(i) * (1. / R_tmp) * (1. / R_tmp));
  tmp.val() = L / R.val();
  return tmp;
}

/** \brief binary division between an AD variable and a double */
template <class T, unsigned int N>
AD<T, N> operator/(AD<T, N> const& L, double R)
{
  AD<T, N> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < L.size(); ++i)
    tmp.dx(i) = T(L.dx(i) / R);
  tmp.val() = L.val() / R;
  return tmp;
}

/** \brief binary division between two AD variables */
template <class T,  class B, unsigned int N>
AD<T, N> operator/(AD<B, N> const& L, AD<T, N> const& R)
{
  AD<T, N> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(((L.dx(i) * R) - (L * R.dx(i))) * (1. / R) * (1. / R));
  tmp.val() = L.val() / R.val();
  return tmp;
}

/********************
  * FANCY FUNCIONS *
*********************/

/** \brief wrapper to standard exp function */
double exp(double x)
{
  return std::exp(x);
}

/** \brief exponent of an AD variable */
template <class T, unsigned int N>
AD<T, N> exp(AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  T A_tmp = A;
  tmp.val() = std::exp(A.val());
  for (unsigned int i=0; i < A.size(); ++i)
    tmp.dx(i) = T(A.dx(i) * exp(A_tmp));
  return tmp;
}

/** \brief wrapper for stander log function */
double log(double A)
{
  return std::log(A);
}

/** \brief logarithm of an AD variable */
template <class T, unsigned int N>
AD<T, N> log(AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  T A_tmp = A;
  tmp.val() = std::log(A.val());
  for (unsigned int i=0; i < A.size(); ++i)
    tmp.dx(i) = T(A.dx(i) / A_tmp);
  return tmp;
}

/** \brief wrapper to standard pow function */
double pow(double A, double e)
{
  return std::pow(A, e);
}

/** \brief AD variable raised to an integer power */
template <class T, unsigned int N>
AD<T, N> pow(AD<T, N> const& A, const int e)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  T A_tmp = A;
  tmp.val() = std::pow(A.val(), e);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(AD<T, N>(e)*A.dx(i)*pow(A_tmp, (double)e-1.));
  return tmp;
}

/** \brief AD variable raised to a double power */
template <class T, unsigned int N>
AD<T, N> pow(AD<T, N> const& A, const double e)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  T A_tmp = A; 
  tmp.val() = pow(A.val(), e);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(e*A.dx(i)*pow(A_tmp, e - 1.));
  return tmp;
}

/** \brief integer raised to an AD power */
template <class T, unsigned int N>
AD<T, N> pow(const int base, AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  T A_tmp = A;
  tmp.val() = std::pow((double)base, A.val());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(std::log((double)base) * pow((double)base, A_tmp)  * A.dx(i));
  return tmp;
}

/** \brief double raised to an AD power */
template <class T, unsigned int N>
AD<T, N> pow(const double base, AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  T A_tmp = A;
  tmp.val() = std::pow((double)base, A.val());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(std::log(base) * pow(base, A_tmp) * A.dx(i));
  return tmp;
}


/** \brief AD variable raised to an AD variable power */
template <class T, unsigned int N>
AD<T, N> pow(AD<T, N> const& A, AD<T, N> const& e)
{
  AD<T, N> tmp;
  unsigned int max = A.size() > e.size() ? A.size() : e.size();
  tmp.resize(max);
  T A_tmp = A;
  T e_tmp = e;
  tmp.val() = std::pow(A.val(), e.val());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(e.dx(i) * log(A_tmp) * pow(A_tmp, e_tmp) +
      e_tmp * A.dx(i) * pow(A_tmp, e_tmp - 1.));
  return tmp;
}

/** \brief wrapper for standard sqrt function */
double sqrt(double A)
{
  return std::sqrt(A);
}

/** \brief square root of an AD variable */
template <class T, unsigned int N>
AD<T, N> sqrt(AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  tmp.val() = std::sqrt(A.val());
  T A_tmp = A;
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = T(.5 * A.dx(i) / sqrt(A_tmp));
  return tmp;
}

/** \brief wrapper for standard sin function */
double sin(double A)
{
  return std::sin(A);
}

/** \brief wrapper for standard cos function */
double cos(double A)
{
  return std::cos(A);
}

/** \brief sin of an AD variable */
template <class T, unsigned int N>
AD<T, N> sin(AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.resize(A.size());
  tmp.val() = std::sin(A.val());
  T A_tmp = A;
  for(unsigned int i = 0; i < tmp.size(); i++)
    tmp.dx(i) = T(cos(A_tmp) * A.dx(i));
  return tmp;
}

/** \brief cos of an AD variable */
template <class T, unsigned int N>
AD<T, N> cos(AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.val() = std::cos(A.val());
  T A_tmp = A;
  tmp.resize(A.size());
  for(unsigned int i = 0; i < tmp.size(); i++)
    tmp.dx(i) = T(-sin(A_tmp) * A.dx(i));
  return tmp;
}

/** \brief tan of an AD variable */
template <class T, unsigned int N>
AD<T, N> tan(AD<T, N> const& A)
{
  AD<T, N> tmp;
  tmp.val() = std::tan(A.val());
  T A_tmp = A;
  tmp.resize(A.size());
  for(unsigned int i = 0; i < tmp.size(); i++)
    tmp.dx(i) = T(A.dx(i) * (1. / (cos(A_tmp) * cos(A_tmp))));
  return tmp;
}

/** \brief absolute value of an AD variable */
template <class T, unsigned int N>
AD<T, N> abs(AD<T, N> const& A)
{
  int sign = A.val() > 0 ? 1 : 0;
  if (sign) return A;
  else return (-A);
}

/**************************
  * COMPARISON OPERATORS *
***************************/

/** \brief double less than an AD variable */
template <class T, unsigned int N>
bool operator<(double L, AD<T, N> const& R)
{
  return L < R.val();
}

/** \brief AD variable less than a double */
template <class T, unsigned int N>
bool operator<(AD<T, N> const& R, double L)
{
  return R.val() < L;
}

/** \brief AD variable less than an AD variable */
template <class T, unsigned int N>
bool operator<(AD<T, N> const& R, AD<T, N> const& L)
{
  return R.val() < L.val();
}

/** \brief double less than or equal to an AD variable */
template <class T, unsigned int N>
bool operator<=(double L, AD<T, N> const& R)
{
  return L <= R.val();
}

/** \brief AD variable less than or equal to a double */
template <class T, unsigned int N>
bool operator<=(AD<T, N> const& R, double L)
{
  return R.val() <= L;
}

/** \brief AD variable less than or equal to an AD variable */
template <class T, unsigned int N>
bool operator<=(AD<T, N> const& R, AD<T, N> const& L)
{
  return R.val() <= L.val();
}

/** \brief double greater than an AD variable */
template <class T, unsigned int N>
bool operator>(double L, AD<T, N> const& R)
{
  return L > R.val();
}

/** \brief AD variable greater than a double */
template <class T, unsigned int N>
bool operator>(AD<T, N> const& R, double L)
{
  return R.val() > L;
}

/** \brief AD variable greater than an AD variable */
template <class T, unsigned int N>
bool operator>(AD<T, N> const& R, AD<T, N> const& L)
{
  return R.val() > L.val();
}

/** \brief double greater than or equal to an AD variable */
template <class T, unsigned int N>
bool operator>=(double L, AD<T, N> const& R)
{
  return L >= R.val();
}

/** \brief AD variable greater than or equal to a double */
template <class T, unsigned int N>
bool operator>=(AD<T, N> const& R, double L)
{
  return R.val() >= L;
}

/** \brief AD variable greater than or equal to an AD variable */
template <class T, unsigned int N>
bool operator>=(AD<T, N> const& R, AD<T, N> const& L)
{
  return R.val() >= L.val();
}
}
#endif
