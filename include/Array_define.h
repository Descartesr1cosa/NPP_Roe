#pragma once
#include <iostream>
#include <vector>
#include <cstdint>
/**
* @brief		Self-defined one-dimensional Array
* @param dim1	size of Array
* @param a1		vector, where data is saved
*/
template<typename T>
class Array1D
{
public:

	//Array1D() {};
	Array1D(){};

	Array1D(int32_t n1, int32_t ngg_in)
	{
		ngg = ngg_in;
		dim1 = n1;
		a1.resize(dim1);
	};

	~Array1D() {};

	T& operator()(int32_t n1)
	{
		return a1[n1 + ngg];
	};
	T& operator=(T a_in)
	{
		this->SetA1(a_in);
		return *this;
	};

	std::vector<T> GetA1() { return a1; };
	void SetA1(std::vector<T>& a1_in, int32_t n1)
	{
		a1 = a1_in;
		dim1 = n1;
		ngg = 0;
	};

	void SetA1(Array1D<T>& a1_in)
	{
		a1 = a1_in.GetA1();
		dim1 = a1_in.Getsize1();
		ngg = a1_in.Getngg();
	}

	int32_t Getsize1() { return dim1; };
	int32_t Getngg() { return ngg; };
	void SetSize(int32_t setdim1, int32_t nggin)
	{
		ngg = nggin;
		dim1 = setdim1;
		a1.resize(setdim1);
	};

private:
	int32_t dim1, ngg;
	std::vector<T> a1;
};

/**
* @brief		Self-defined two-dimensional Array

* @param dim1	size1 of Array
* @param dim2	size2 of Array
* @param a2		vector, where data is saved
*/
template<typename T>
class Array2D
{
public:

	//Array2D() {};
	Array2D() = default;

	Array2D(int32_t n1, int32_t n2, int32_t ngg_in)
	{
		dim1 = n1;
		dim2 = n2;
		ngg = ngg_in;
		a2.resize(dim1 * dim2);
	};

	~Array2D() {};

	T& operator()(int32_t n1, int32_t n2)
	{
		return a2[(n1+ngg) * dim2 + n2];
	};

	std::vector<T> GetA2() { return a2; };
	void SetA2(std::vector<T>& a2_in, int32_t n1, int32_t n2, int32_t ngg_in)
	{
		a2 = a2_in;
		dim1 = n1;
		dim2 = n2;
		ngg = ngg_in;
	};
	void SetA2(Array2D<T>& a2_in)
	{
		a2 = a2_in.GetA2();
		dim1 = a2_in.Getsize1();
		dim2 = a2_in.Getsize2();
		ngg = a2_in.Getngg();
	}

	int32_t Getsize1() { return dim1; };
	int32_t Getsize2() { return dim2; };
	int32_t Getngg() { return ngg; };

	void SetSize(int32_t setdim1, int32_t setdim2, int32_t ngg_in)
	{
		dim1 = setdim1;
		dim2 = setdim2;
		ngg = ngg_in;
		a2.resize(setdim1 * setdim2);
	};

private:
	int32_t dim1, dim2;
	int32_t ngg;
	std::vector<T> a2;
};

//typedef::
typedef Array1D<int> int1D;
typedef Array1D<double> double1D;

typedef Array2D<int> int2D;
typedef Array2D<double> double2D;

