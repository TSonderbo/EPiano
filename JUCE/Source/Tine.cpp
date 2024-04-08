/*
  ==============================================================================

    Tine.cpp
    Created: 8 Mar 2024 10:33:32am
    Author:  Sonderbo

  ==============================================================================
*/

#include "Tine.h"

Tine::Tine()
{
	isActive = false;
	isStopped = true;
}

void Tine::prepareToPlay(double sampleRate, int tineNumber)
{
	L = TINE_LENGTHS[tineNumber];
	k = 1 / (sampleRate * config::oversampling);
	kSq = k * k;
	this->sampleRate = sampleRate * config::oversampling;

	h = sqrtf((4 * sigma_1 * k + sqrtf(pow((4 * sigma_1 * k), 2) + 16 * kappa * kappa * k * k)) / 2);

	N = static_cast<int>(L / h);
	N_frac = L / h;
	alpha = N_frac - N;

	M_u = ceil(0.5 * N);
	M_w = floor(0.5 * N);

	hSq = h * h;

	mu = (kappa * k) / hSq;
	muSq = mu * mu;

	freq = 1.426f * (juce::MathConstants<float>::pi * K) / (16 * L * L) * sqrtf(E / rho);

	uStates = std::vector<std::vector<float>>(3, std::vector<float>(M_u + 1, 0));
	wStates = std::vector<std::vector<float>>(3, std::vector<float>(M_w + 1, 0));

	u.resize(3);
	w.resize(3);

	for (int i = 0; i < 3; i++)
	{
		//TODO Change to allocate based on a calculated max requirement instead of hard-coded value
		uStates[i].reserve(size_t(100));
		wStates[i].reserve(size_t(100));

		u[i] = uStates[i].begin();
		w[i] = wStates[i].begin();
	}

	h_contact.resize(M_u + 1, 0.0f);

	h_contact[static_cast<int>(0.5 * M_u)] = 1;
	
	hammer.prepareToPlay(sampleRate, tineNumber, h_contact, M_u);
	h_ratio = hammer.getMass() / (rho * A * L);
	
	inactiveTimer = config::tine::inactiveTimer * sampleRate / config::oversampling;

	//TODO add Iterm and J calculation
	calculateInterpolation();
}

void Tine::startNote(float velocity)
{
	hammer.beginHammer(velocity);
	isActive = true;
	inactiveTimerIterator = inactiveTimer; //Inactive timer in samples
	isStopped = false;
}

void Tine::stopNote()
{
	isStopped = true;
}

float Tine::processSample()
{
	if (isActive == false)
		return 0.0f;
	if (isStopped == true)
	{
		inactiveTimerIterator = inactiveTimerIterator - 1;
		if (inactiveTimerIterator <= 0)
		{
			isActive = false;
			resetScheme();
			return 0.0f;
		}
	}

	float sample;

	for (size_t i = 0; i < config::oversampling; i++)
	{
		calculateScheme();
	}

	sample = w[1][0]; //Read output at tip

	float gain = 1.0f;

	if (isStopped == true)
		gain = (1.0f / inactiveTimer) * inactiveTimerIterator;

	return sample * gain;
}

void Tine::setLength(float scalar)
{
	float freq_new = freq * powf(2, scalar * 2 / 12);
	float L_new = 1.426f * sqrtf(((juce::MathConstants<float>::pi * (r / 2.0f) * K) / freq_new) / 16.0f);

	int Nprev = N;

	N = static_cast<int>(L_new / h);

	N_frac = L_new / h;
	alpha = N_frac - N;

	updateGridpoints(Nprev);
}

void Tine::calculateScheme()
{
	float force = hammer.calculateForce(u[1]);
	//TODO precalculate coefficients for update equations
	//Update loop - excluding outer and inner boundaries
	for (int l = 2; l <= M_u - 2; l++)
	{
		u[2][l] = ((2 - 6 * muSq - ((4 * sigma_1 * k) / hSq)) * u[1][l]
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (u[1][l + 1] + u[1][l - 1])
			- muSq * (u[1][l + 2] + u[1][l - 2]) 
			+ (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * u[0][l]
			- ((2 * sigma_1 * k) / hSq) * (u[0][l + 1] + u[0][l - 1]))
			+ (h_ratio * h_contact[l] * force * kSq)
			/ (1 + sigma_0 * k);

		w[2][l] = ((2 - 6 * muSq - ((4 * sigma_1 * k) / hSq)) * w[1][l]
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (w[1][l + 1] + w[1][l - 1])
			- muSq * (w[1][l + 2] + w[1][l - 2]) 
			+ (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * w[0][l]
			- ((2 * sigma_1 * k) / hSq) * (w[0][l + 1] + w[0][l - 1]))
			/ (1 + sigma_0 * k);
	}
	
	if (M_w > M_u) //In case of uneven grid lengths - W will always be greater
	{
		int l = M_w - 2;

		w[2][l] = ((2 - 6 * muSq - ((4 * sigma_1 * k) / hSq)) * w[1][l]
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (w[1][l + 1] + w[1][l - 1])
			- muSq * (w[1][l + 2] + w[1][l - 2])
			+ (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * w[0][l]
			- ((2 * sigma_1 * k) / hSq) * (w[0][l + 1] + w[0][l - 1]))
			/ (1 + sigma_0 * k);
	}

	//u inner boundary
	u[2][M_u - 1] = ((2 - 6 * muSq - ((4 * sigma_1 * k) / hSq)) * u[1][M_u - 1]
		+ (-J[1] * muSq + (2 * sigma_1 * k) / hSq) * u[1][M_u]
		+ (4 * muSq + (2 * sigma_1 * k) / hSq) * u[1][M_u - 2]
		- muSq * w[1][M_w]
		- J[3] * muSq * w[1][M_w - 1]
		+ (-1 + sigma_0 * k + (4 * sigma_1 * k) / hSq) * u[0][M_u - 1]
		- ((2 * sigma_1 * k) / hSq) * (u[0][M_u - 2] + u[0][M_u]))
		/ (1 + sigma_0 * k);

	u[2][M_u] = ((2 - J[0] * muSq + (Iterm - 2) * (2 * sigma_1 * k) / hSq) * u[1][M_u]
		+ (4 * muSq + (2 * sigma_1 * k) / hSq) * u[1][M_u - 1]
		- muSq * u[1][M_u - 2]
		- (J[1] * muSq - (2 * sigma_1 * k) / hSq) * w[1][M_w]
		- (J[2] * muSq + Iterm) * w[1][M_w - 1]
		- (J[3] * muSq) * w[1][M_w - 2]
		+ (-1 + sigma_0 * k - (Iterm - 2) * ((2 * sigma_1 * k) / hSq)) * u[0][M_u]
		- ((2 * sigma_1 * k) / hSq) * (u[0][M_u - 1] + w[0][M_w])
		- (-Iterm * ((2 * sigma_1 * k) / hSq)) * w[0][M_w - 1])
		/ (1 + sigma_0 * k);

	//w inner boundary
	w[2][M_w - 1] = ((2 - 6 * muSq - ((4 * sigma_1 * k) / hSq)) * w[1][M_w - 1]
		+ (-J[1] * muSq + (2 * sigma_1 * k) / hSq) * w[1][M_w]
		+ (4 * muSq + (2 * sigma_1 * k) / hSq) * w[1][M_w - 2]
		- muSq * u[1][M_u]
		- J[3] * muSq * u[1][M_u - 1]
		+ (-1 + sigma_0 * k + (4 * sigma_1 * k) / hSq) * w[0][M_w - 1]
		- ((2 * sigma_1 * k) / hSq) * (w[0][M_w - 2] + w[0][M_w]))
		/ (1 + sigma_0 * k);

	u[2][M_u] = ((2 - J[0] * muSq + (Iterm - 2) * (2 * sigma_1 * k) / hSq) * w[1][M_w]
		+ (4 * muSq + (2 * sigma_1 * k) / hSq) * w[1][M_w - 1]
		- muSq * w[1][M_w - 2]
		- (J[1] * muSq - (2 * sigma_1 * k) / hSq) * u[1][M_u]
		- (J[2] * muSq + Iterm) * u[1][M_u - 1]
		- (J[3] * muSq) * u[1][M_u - 2]
		+ (-1 + sigma_0 * k - (Iterm - 2) * ((2 * sigma_1 * k) / hSq)) * w[0][M_w]
		- ((2 * sigma_1 * k) / hSq) * (w[0][M_w - 1] + u[0][M_u])
		- (-Iterm * ((2 * sigma_1 * k) / hSq)) * u[0][M_u - 1])
		/ (1 + sigma_0 * k);

	//Free boundary update
	w[2][1] = ((2 - 5 * muSq - ((4 * sigma_1 * k) / hSq)) * w[1][1]
		+ ((2 * sigma_1 * k) / (hSq)+2 * muSq) * w[1][0]
		+ (4 * muSq + (2 * sigma_1 * k) / (hSq)) * w[1][2]
		- muSq * w[1][3]
		+ (-1 + sigma_0 * k + (4 * sigma_1 * k) / hSq) * w[0][1]
		- ((2 * sigma_1 * k) / (hSq)) * (w[0][2] + w[0][0]))
		/ (1 + sigma_0 * k);

	w[2][0] = ((2 - 2 * muSq - ((4 * sigma_1 * k) / hSq)) * w[1][0]
		+ (4 * muSq + ((4 * sigma_1 * k) / hSq)) * w[1][1]
		- muSq * 2 * w[1][2]
		+ (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * w[0][0]
		- ((2 * sigma_1 * k) / hSq) * (2 * w[0][1]))
		/ (1 + sigma_0 * k);

	updateStates();
}

void Tine::updateStates()
{
	auto uPrev = u[0];
	u[0] = u[1];
	u[1] = u[2];
	u[2] = uPrev;

	auto wPrev = w[0];
	w[0] = w[1];
	w[1] = w[2];
	w[2] = wPrev;
}

bool Tine::getIsActive()
{
	return isActive;
}

void Tine::resetScheme()
{
	for (int i = 0; i < 3; i++)
	{
		std::fill(uStates[i].begin(), uStates[i].end(), 0);
		std::fill(wStates[i].begin(), wStates[i].end(), 0);
		u[i] = uStates[i].begin();
		w[i] = wStates[i].begin();
	}
}

void Tine::updateGridpoints(int Nprev)
{
	if (abs(Nprev - N) > 1)
		throw std::invalid_argument("Addition of multiple new grid points is invalid");

	if (Nprev < N) //Addition
	{
		if (N % 2 == 0) //Even - Add to w
		{
			addGridpoint(wStates, uStates);
		}
		else //Odd - Add to u
		{
			addGridpoint(uStates, wStates);
		}
	}
	else //Removal
	{
		if (N % 2 == 0) //Even  - Remove from u
		{
			removeGridpoint(uStates);
		}
		else //Odd - Remove from w
		{
			removeGridpoint(wStates);
		}
	}
	M_u = uStates[0].size() - 1;
	M_w = wStates[0].size() - 1;
}

void Tine::removeGridpoint(std::vector<std::vector<float>>& grid)
{
	for (int i = 0; i < grid.size(); i++)
	{
		grid[i].pop_back();
	}
	//TODO Displacement correction
}

void Tine::addGridpoint(std::vector<std::vector<float>>& main, std::vector<std::vector<float>>& sec)
{
	interp[0] = -alpha * (alpha + 1.0) / ((alpha + 2.0) * (alpha + 3.0));
	interp[1] = 2.0 * alpha / (alpha + 2.0);
	interp[2] = 2.0 / (alpha + 2.0);
	interp[3] = -2.0 * alpha / ((alpha + 3.0) * (alpha + 2.0));

	for (int i = 0; i < main.size(); i++)
	{
		int M = main[i].size() - 1;
		int S = sec[i].size() - 1;
		float v[4] = { main[i][M - 1], main[i][M], sec[i][S], sec[i][S - 1]}; //Build v vector

		float point = 0.0f;
		for (int j = 0; j < 4; j++)
			point += interp[j] * v[j];

		main[i].push_back(point);
	}
}

void Tine::calculateInterpolation()
{
	Iterm = (alpha - 1) / (alpha + 1);
	float itermSq = Iterm * Iterm;
	J[0] = itermSq - 4 * Iterm + 6;
	J[1] = Iterm - 4;
	J[2] = -itermSq + 4 * Iterm + 1;
	J[3] = -Iterm;
}