#pragma once

#include <cassert>
#include <vector>

template<typename T>
class TVecRef {
public:
  TVecRef(int cols, T* ptr) : m_cols(cols), m_ptr(ptr) {
  }

  int cols() const { return m_cols; }
  int size() const { return cols(); }

  bool is_same_shape(const TVecRef<T>& other) const {
    return cols() == other.cols();
  }

  T* data() { return m_ptr; }
  const T* data() const { return m_ptr; }

  T& operator[](int i) {
    assert(i >= 0);
    assert(i < m_cols);
    return m_ptr[i];
  }

private:
  int m_cols = 0;
  T* m_ptr = nullptr;
};

template<typename T>
class TGridRef {
public:
  TGridRef(int rows, int cols, T* ptr) : m_rows(rows), m_cols(cols), m_ptr(ptr) {
  }

  int cols() const { return m_cols; }
  int rows() const { return m_rows; }
  int size() const { return rows() * cols(); }

  bool is_same_shape(const TGridRef<T>& other) const {
    return rows() == other.rows() and cols() == other.cols();
  }

  T* data() { return m_ptr; }
  const T* data() const { return m_ptr; }

  TVecRef<T> operator[](int i) {
    assert(i >= 0);
    assert(i < m_rows);
    return TVecRef<T>(m_cols, &m_ptr[i * m_cols]);
  }

private:
  int m_rows = 0;
  int m_cols = 0;
  T* m_ptr = nullptr;
};

template<typename T>
class TGrid {
public:
  TGrid() = default;
  TGrid(int rows, int cols, T val) : m_rows(rows), m_cols(cols) {
    m_array.resize(rows*cols, val);
  }

  void reshape(int rows, int cols) {
    assert(rows >= 0);
    assert(cols >= 0);
    m_rows = rows;
    m_cols = cols;
    m_array.clear();
    m_array.resize(rows*cols);
  }

  int rows() const { return m_rows; }
  int cols() const { return m_cols; }
  int size() const { return rows() * cols(); }

  bool is_same_shape(const TGrid<T>& other) const {
    return rows() == other.rows() and cols() == other.cols();
  }

  T* data() { return m_array.data(); }
  const T* data() const { return m_array.data(); }

  TVecRef<T> operator[](int i) {
    assert(i >= 0);
    assert(i < m_rows);
    return TVecRef<T>(m_cols, &m_array[i * m_cols]);
  }

private:
  int m_rows = 0;
  int m_cols = 0;
  std::vector<T> m_array;
};

template<typename T>
class TCube {
public:
  TCube() = default;
  TCube(int depths, int rows, int cols, T val) : m_depths(depths), m_rows(rows), m_cols(cols) {
    m_array.resize(depths*rows*cols, val);
  }

  void reshape(int depths, int rows, int cols) {
    assert(depths >= 0);
    assert(rows >= 0);
    assert(cols >= 0);
    m_depths = depths;
    m_rows = rows;
    m_cols = cols;
    m_array.clear();
    m_array.resize(depths*rows*cols);
  }

  int depths() const { return m_depths; }
  int rows() const { return m_rows; }
  int cols() const { return m_cols; }
  int size() const { return depths() * rows() * cols(); }

  bool is_same_shape(const TCube<T>& other) const {
    return depths() == other.depths() and rows() == other.rows() and cols() == other.cols();
  }

  T* data() { return m_array.data(); }
  const T* data() const { return m_array.data(); }

  TGridRef<T> operator[](int depth) {
    assert(depth >= 0);
    assert(depth < m_depths);
    T *p = m_ptr ? m_ptr : m_array.data();
    return TGridRef<T>(m_rows, m_cols, &p[depth*m_rows*m_cols]);
  }

private:
  int m_depths = 0;
  int m_rows = 0;
  int m_cols = 0;
  T* m_ptr = nullptr;
  std::vector<T> m_array;
};

template<typename T>
void elem_copy(TCube<T>& a, const T* src)
{
  int n = a.size();
  T* dst = a.data();
  for(int i = 0; i < n; ++i)
    dst[i] = src[i];
}

template<typename T>
bool elem_add(TCube<T>& a, const TCube<T>& b)
{
  if (!a.is_same_shape(b))
    return false;
  int n = a.size();
  T* dst = a.data();
  const T* src = b.data();
  for(int i = 0; i < n; ++i)
    dst[i] += src[i];
  return true;
}

template<typename T>
bool elem_sub(TCube<T>& a, TCube<T>& b)
{
  if (!a.is_same_shape(b))
    return false;
  int n = a.size();
  T* dst = a.data();
  const T* src = b.data();
  for(int i = 0; i < n; ++i)
    dst[i] -= src[i];
  return true;
}

template<typename T, typename Pred>
bool elem_set(TCube<T>& a, const TCube<T>& b)
{
  if (!a.is_same_shape(b))
    return false;
  int n = a.size();
  T* dst = a.data();
  const T* src = b.data();
  for(int i = 0; i < n; ++i)
    dst[i] = src[i];
  return true;
}

template<typename T, typename Pred>
bool elem_set(TCube<T>& a, const TCube<T>& b, Pred f)
{
  if (!a.is_same_shape(b))
    return false;
  int n = a.size();
  T* dst = a.data();
  const T* src = b.data();
  for(int i = 0; i < n; ++i)
    dst[i] = f(src[i]);
  return true;
}

template<typename T, typename Pred>
bool elem_add(TCube<T>& a, const TCube<T>& b, Pred f)
{
  if (!a.is_same_shape(b))
    return false;
  int n = a.size();
  T* dst = a.data();
  const T* src = b.data();
  for(int i = 0; i < n; ++i)
    dst[i] += f(src[i]);
  return true;
}

template<typename T, typename Pred>
bool elem_sub(TCube<T>& a, const TCube<T>& b, Pred f)
{
  if (!a.is_same_shape(b))
    return false;
  int n = a.size();
  T* dst = a.data();
  const T* src = b.data();
  for(int i = 0; i < n; ++i)
    dst[i] -= f(src[i]);
  return true;
}


typedef TVecRef<double> dVecRef;
typedef TVecRef<float> fVecRef;
typedef TVecRef<int> iVecRef;

typedef TGridRef<double> dGridRef;
typedef TGridRef<float> fGridRef;
typedef TGridRef<int> iGridRef;

typedef TGrid<double> dGrid;
typedef TGrid<float> fGrid;
typedef TGrid<int> iGrid;

typedef TCube<double> dCube;
typedef TCube<float> fCube;
typedef TCube<int> iCube;
