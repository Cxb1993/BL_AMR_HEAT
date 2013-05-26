/*
 * FortranHelper.hpp
 *
 *  Created on: 2013/04/10
 *      Author: homu
 */

#ifndef FORTRANHELPER_HPP_
#define FORTRANHELPER_HPP_

#include "REAL.H"
#include "SPACE.H"

//#include <vector>
#include <cassert>

class FortDataIndexer {
protected:
	enum { MAX_DIMENSION = 4, };

	int dimension;

	int first[MAX_DIMENSION];
	int last[MAX_DIMENSION];
	int stride[MAX_DIMENSION];

	int total;
	int offset;

public:
	FortDataIndexer(int dim) : dimension(dim) {
		assert(dim<=MAX_DIMENSION);

		total = 0;
		offset = 0;
	}

	inline int getDimension() const {
		return dimension;
	}

	inline int size(int dim) const {
		return last[dim]+1 - first[dim];
	}

	inline const int* low() const {
		return first;
	}
	inline int low(int dim) const {
		return first[dim];
	}
	inline const int* high() const {
		return last;
	}
	inline int high(int dim) const {
		return last[dim];
	}

	inline void setBound(int dim, int low, int high) {
		first[dim] = low;
		last[dim] = high;
	}

protected:
	void pack() {
		int s = 1;
		offset = 0;
		for(int i=0; i<dimension; i++) {
			stride[i] = s;
			offset += first[i] * stride[i];
			s *= size(i);
		}
		total = s;
	}

	inline int index(int i) const {
		return i*stride[0] - offset;
	}
	inline int index(int i, int j) const {
		return i*stride[0] + j*stride[1] - offset;
	}
	inline int index(int i, int j, int k) const {
		return i*stride[0] + j*stride[1] + k*stride[2] - offset;
	}
	inline int index(int i, int j, int k, int l) const {
		return i*stride[0] + j*stride[1] + k*stride[2] + l*stride[3] - offset;
	}
};

template<typename T>
class FortArrayRefBase : public FortDataIndexer
{
public:
	FortArrayRefBase(int dimension, T* pData) :
		FortDataIndexer(dimension),
		data(pData)
	{}

	inline T* getData() { return data; }
	inline const T* getData() const { return data; }
	inline void setData(T* pData) { data = pData; }
protected:
	T *data;

	inline T& getValue(int i) {
		return data[index(i)];
	}
	inline T& getValue(int i, int j) {
		return data[index(i,j)];
	}
	inline T& getValue(int i, int j, int k) {
		return data[index(i,j,k)];
	}
	inline T& getValue(int i, int j, int k, int l) {
		return data[index(i,j,k,l)];
	}

};

// forward
template<typename T, int Dimension>
class FortArrayWrap;

// 1-D
template<typename T>
class FortArrayWrap<T,1> : public FortArrayRefBase<T>
{
public:
	enum { ArrayDIM = 1, };
	typedef FortArrayRefBase<T> super_type;

	FortArrayWrap(T *data, int low, int high) :
		super_type(ArrayDIM, data)
	{
		super_type::setBound(0, low, high);

		super_type::pack();
	}

	inline T& operator() (int i) {
		return super_type::getValue(i);
	}
};

// 2-D
template<typename T>
class FortArrayWrap<T,2> : public FortArrayRefBase<T>
{
public:
	enum { ArrayDIM = 2, };
	typedef FortArrayRefBase<T> super_type;

	FortArrayWrap(T *data,
			int low1, int high1,
			int low2, int high2) :
		super_type(ArrayDIM, data)
	{
		super_type::setBound(0, low1, high1);
		super_type::setBound(1, low2, high2);

		super_type::pack();
	}

	inline T& operator() (int i, int j) {
		return super_type::getValue(i, j);
	}
};

// 3-D
template<typename T>
class FortArrayWrap<T,3> : public FortArrayRefBase<T>
{
public:
	enum { ArrayDIM = 3, };
	typedef FortArrayRefBase<T> super_type;

	FortArrayWrap(T *data,
			int low1, int high1,
			int low2, int high2,
			int low3, int high3) :
		super_type(ArrayDIM, data)
	{
		super_type::setBound(0, low1, high1);
		super_type::setBound(1, low2, high2);
		super_type::setBound(2, low3, high3);

		super_type::pack();
	}

	inline T& operator() (int i, int j, int k) {
		return super_type::getValue(i, j, k);
	}
};

// 4-D
template<typename T>
class FortArrayWrap<T,4> : public FortArrayRefBase<T>
{
public:
	enum { ArrayDIM = 4, };
	typedef FortArrayRefBase<T> super_type;

	FortArrayWrap(T *data,
			int low1, int high1,
			int low2, int high2,
			int low3, int high3,
			int low4, int high4) :
		super_type(ArrayDIM, data)
	{
		super_type::setBound(0, low1, high1);
		super_type::setBound(1, low2, high2);
		super_type::setBound(2, low3, high3);
		super_type::setBound(3, low4, high4);

		super_type::pack();
	}

	inline T& operator() (int i, int j, int k, int l) {
		return super_type::getValue(i, j, k, l);
	}
};


#define FORTRAN_ARRAY_2D(name, data, lo1,lo2, hi1,hi2) \
	FortArrayWrap<Real,2> name((data), (lo1),(hi1), (lo2),(hi2))
#define FORTRAN_ARRAY_3D(name, data, lo1,lo2,lo3, hi1,hi2,hi3) \
	FortArrayWrap<Real,3> name((data), (lo1),(hi1), (lo2),(hi2), (lo3),(hi3))



template<typename T, int Dimension>
class FortDataRefBase : public FortDataIndexer
{
public:
	FortDataRefBase() :
		FortDataIndexer(Dimension),
		data(NULL) { }

	inline T* getData() {
		return data;
	}
	inline const T* getData() const {
		return data;
	}

	inline void setData(T* pData) {
		data = pData;
	}

//	inline const int getDimension() const {
//		return Dimension;
//	}
//
//	inline int size(int dim) const {
//		return last[dim]+1 - first[dim];
//	}
//
//	inline const int* low() const {
//		return first;
//	}
//	inline int low(int dim) const {
//		return first[dim];
//	}
//	inline const int* high() const {
//		return last;
//	}
//	inline int high(int dim) const {
//		return last[dim];
//	}
//
//	inline void setBound(int dim, int low, int high) {
//		first[dim] = low;
//		last[dim] = high;
//	}
protected:
//	void pack() {
//		int s = 1;
//		offset = 0;
//		for(int i=0; i<Dimension; i++) {
//			stride[i] = s;
//			offset += first[i] * stride[i];
//			s *= size(i);
//		}
//		dataLength = s;
//	}

protected:
//	int first[Dimension];
//	int last[Dimension];
//	int stride[Dimension];
//
//	int offset;
//
//	int dataLength;

	T *data;
};

template<typename T, int Dimension>
class FortDataRef;



template<typename T, int Dimension>
class FabDataWrap;

template<typename T>
class FabDataWrap<T,2> : public FortDataRefBase<T,2+1>
{
public:
	typedef FortDataRefBase<T,2+1> super_type;

	FabDataWrap(T *dataPtr,
			int low1, int low2, int high1, int high2,
			int numComp=1, int iComp=0)
	{
		this->data = dataPtr;

		this->first[0] = low1;
		this->first[1] = low2;
		this->last[0] = high1;
		this->last[1] = high2;

		this->first[2] = iComp;
		this->last[2] = iComp+numComp-1;

		this->pack();
	}

//	inline T& operator()(int i, int j) {
//		return this->data[index(i,j)];
//	}
	inline T& operator()(int i, int j, int comp=0) {
		return this->data[index(i,j,comp)];
	}

protected:
//	inline int index(int i, int j) const {
//		return i*this->stride[0] + j*this->stride[1] - this->offset;
//	}
	inline int index(int i, int j, int comp) const {
//		return i*this->stride[0] + j*this->stride[1] +
//				comp*this->stride[2] - this->offset;
		return super_type::index(i, j, comp);
	}
};

//template<typename T>
//class FabDataWrap<T,3> : public FortDataRefBase<T, 3>
//{
//public:
//
//	FabDataWrap(T *dataPtr,
//			int low1, int low2, int low3,
//			int high1, int high2, int high3)
//	{
//		data = dataPtr;
//
//		first[0] = low1;
//		first[1] = low2;
//		first[2] = low3;
//		last[0] = high1;
//		last[1] = high2;
//		last[2] = high3;
//
//		pack();
//	}
//
//	inline T& operator()(int i, int j, int k) {
//		return data[index(i,j,k)];
//	}
//
//protected:
//	inline int index(int i, int j, int k) const {
//		return i*stride[0] + j*stride[1] + k*stride[2] - offset;
//	}
//};


typedef FabDataWrap<Real, BL_SPACEDIM> FabDataWrapper;


#endif /* FORTRANHELPER_HPP_ */
