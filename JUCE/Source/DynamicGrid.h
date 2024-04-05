/*
  ==============================================================================

	DynamicGrid.h
	Created: 4 Apr 2024 9:46:46am
	Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <concepts>
#include <JuceHeader.h>

template <typename T>
class DynamicGrid
{
public:

	void resize(T N_frac_new)
	{
		N_frac = N_frac_new;

		int N_new = static_cast<int>(N_frac_new);
		N = N_new;

		int M_u = std::ceil(0.5 * N);
		int M_w = std::floor(0.5 * N);

		u.resize(M_u + 1, 0);
		w.resize(M_w + 1, 0);

		size = M_u + 1 + M_w + 1;
	}

	int updateGrid(T N_frac_new)
	{
		int N_new = static_cast<int>(N_frac_new);

		if (N_new == N)
			return 0;
		if (N_new % N > 1)
			throw std::invalid_argument("Grid scaled too fast");

		N_frac = N_frac_new;

		T a = N_frac_new - N_new;

		T a1 = -((a * (a + 1)) / ((a + 2) * (a + 3)));
		T a2 = (2 * a) / (a + 2);
		T a3 = 2 / (a + 2);
		T a4 = -((2 * a) / ((a + 3) * (a + 2)));

		T I[4] = { a1, a2, a3, a4 };

		T* u_rit = u.rbegin();
		T* w_rit = w.rbegin();

		T v[4] = { u_rit[1], u_rit[0], w_rit[0], w_rit[1] };

		if (N_new > N)
		{
			if (N_new % 2 == 0) //Even
			{
				std::reverse(std::begin(I), std::end(I));

				T sum = multiplyVec(I, v);

				w.push_back(sum);
			}
			else //Odd
			{
				T sum = multiplyVec(I, v);

				u.push_back(sum);
			}
			size++;
			return 1;
		}
		else if (N_new < N)
		{
			if (N_new % 2 == 0) //Even
			{
				u.pop_back();
			}
			else //Odd
			{
				w.pop_back();
			}
			size--;
			return -1;
		}
		jassertfalse; //This point shouldn't be reachable
	}

	T operator[](int idx)
	{
		if (idx >= size)
			throw std::out_of_range("Index was out of range - index: " + std::to_string(idx) + " array size: " + std::to_string(size));

		if (idx < u.size() - 1) //Is element in u except last element
		{
			return u[idx];
		}
		else if (idx > u.size()) //Is element in w except first element
		{
			auto w_rit = w.rbegin();
			return w_rit[idx - u.size()];
		}
		else ()
		{
			//
			//If it's either last element of u or first element of w
		}



	}

	void set(T value)
	{

	}

private:

	//Vectors
	std::vector<T> u;
	std::vector<T> w;

	//Grid variables
	int N;
	T N_frac;

	int size;

	inline T multiplyVec(T I[4], T v[4])
	{
		T sum = 0.0;
		for (int i = 0; i < 4; i++)
		{
			sum += I[i] * v[i];
		}
		return sum;
	}
};