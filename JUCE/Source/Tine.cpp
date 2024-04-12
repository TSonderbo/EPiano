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
	isStopped = true;
	isPrepared = false;
}

void Tine::prepareToPlay(double sampleRate)
{
	k = 1 / (sampleRate * config::oversampling);
	kSq = k * k;
	this->sampleRate = sampleRate * config::oversampling;

	h = sqrtf((4 * sigma_1 * k + sqrtf(pow((4 * sigma_1 * k), 2) + 16 * kappa * kappa * k * k)) / 2);

	hSq = h * h;

	mu = (kappa * k) / hSq;
	muSq = mu * mu;

	isPrepared = true;
}

void Tine::prepareGrid(float freq)
{
	L = calculateLength(freq);

	N = static_cast<int>(L / h);
	N_frac = L / h;
	alpha = N_frac - N;

	calculateInterpolation();

	M_u = ceil(0.5 * N);
	M_w = floor(0.5 * N);

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

	h_contact[static_cast<int>(0.33 * M_u)] = 1;

	hammer.prepareToPlay(sampleRate, h_contact, M_u);
	h_ratio = hammer.getMass() / (rho * A * L);
}

void Tine::noteStarted()
{
	jassert(currentlyPlayingNote.isValid());
	jassert(currentlyPlayingNote.keyState == juce::MPENote::keyDown
		|| currentlyPlayingNote.keyState == juce::MPENote::keyDownAndSustained);
	jassert(isPrepared);

	float freq = currentlyPlayingNote.getFrequencyInHertz();
	float velocity = currentlyPlayingNote.noteOnVelocity.asUnsignedFloat() * config::mpe::maxInputVelocity;

	prepareGrid(freq);

	hammer.beginHammer(velocity);

	isStopped = false;
}

void Tine::noteStopped(bool allowTailOff)
{
	isStopped = true;
}

void Tine::notePressureChanged()
{

}
void Tine::noteTimbreChanged()
{

}
void Tine::noteKeyStateChanged()
{

}

void Tine::notePitchbendChanged()
{
	//float scalar = currentlyPlayingNote.pitchbend.asSignedFloat();
	float baseFreq = juce::MidiMessage::getMidiNoteInHertz(currentlyPlayingNote.initialNote);
	float freq = currentlyPlayingNote.getFrequencyInHertz();


	//float freq_new = freq * powf(2, scalar * 2 / 12);
	float L_new = calculateLength(freq);

	int Nprev = N;

	N = static_cast<int>(L_new / h);

	N_frac = L_new / h;
	alpha = N_frac - N;
	calculateInterpolation();

	updateGridpoints(Nprev);
}

float Tine::processSample()
{
	float sample;

	for (size_t i = 0; i < config::oversampling; i++)
	{
		calculateScheme();
	}

	sample = w[1][0]; //Read output at tip

	//TODO pickup filtering goes here...

	return limit(sample * 100.0f);
}

void Tine::renderNextBlock(juce::AudioBuffer<float>& outputBuffer, int startSample, int numSamples)
{
	while (--numSamples >= 0)
	{
		auto currentSample = processSample();

		for (auto i = outputBuffer.getNumChannels(); --i >= 0;)
			outputBuffer.addSample(i, startSample, currentSample);

		++startSample;
	}
}

void Tine::calculateScheme()
{
	float force = hammer.calculateForce(u[1]);
	//TODO precalculate coefficients for update equations
	//Update loop - excluding outer and inner boundaries
	for (int l = 2; l <= M_w - 2; l++)
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

	if (M_u > M_w) //In case of uneven grid lengths - U will always be greater
	{
		int l = M_u - 2;

		u[2][l] = ((2 - 6 * muSq - ((4 * sigma_1 * k) / hSq)) * u[1][l]
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (u[1][l + 1] + u[1][l - 1])
			- muSq * (u[1][l + 2] + u[1][l - 2])
			+ (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * u[0][l]
			- ((2 * sigma_1 * k) / hSq) * (u[0][l + 1] + u[0][l - 1]))
			/ (1 + sigma_0 * k);
	}

	//Inner boundaries
	u[2][M_u - 1] = ((2 - 6 * muSq - (4 * sigma_1 * k) / hSq) * u[1][M_u - 1]
		+ (-J[1] * muSq + (2 * sigma_1 * k) / hSq) * u[1][M_u]
		+ (4 * muSq + (2 * sigma_1 * k) / hSq) * u[1][M_u - 2]
		- muSq * u[1][M_u - 3]
		- muSq * w[1][M_w]
		- J[3] * muSq * w[1][M_w - 1]
		+ (-1 + sigma_0 * k + (4 * sigma_1 * k) / hSq) * u[0][M_u - 1]
		- ((2 * sigma_1 * k) / hSq) * (u[0][M_u - 2] + u[0][M_u])
		) / (1 + sigma_0 * k);

	u[2][M_u] = ((2 - J[0] * muSq + (Iterm - 2) * ((2 * sigma_1 * k) / hSq)) * u[1][M_u]
		+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * u[1][M_u - 1]
		- muSq * u[1][M_u - 2]
		+ (((2 * sigma_1 * k) / hSq) - J[1] * muSq) * w[1][M_w]
		+ (J[3] * ((2 * sigma_1 * k) / hSq) - J[2] * muSq) * w[1][M_w - 1]
		- (J[3] * muSq) * w[1][M_w - 2]
		+ (-1 + sigma_0 * k - (Iterm - 2) * ((2 * sigma_1 * k) / hSq)) * u[0][M_u]
		- ((2 * sigma_1 * k) / hSq) * (u[0][M_u - 1] + w[0][M_w])
		- (J[3] * ((2 * sigma_1 * k) / hSq)) * w[0][M_w - 1]
		) / (1 + sigma_0 * k);

	w[2][M_w - 1] = ((2 - 6 * muSq - (4 * sigma_1 * k) / hSq) * w[1][M_w - 1]
		+ (-J[1] * muSq + (2 * sigma_1 * k) / hSq) * w[1][M_w]
		+ (4 * muSq + (2 * sigma_1 * k) / hSq) * w[1][M_w - 2]
		- muSq * w[1][M_w - 3]
		- muSq * u[1][M_u]
		- J[3] * muSq * u[1][M_u - 1]
		+ (-1 + sigma_0 * k + (4 * sigma_1 * k) / hSq) * w[0][M_w - 1]
		- ((2 * sigma_1 * k) / hSq) * (w[0][M_w - 2] + w[0][M_w])
		) / (1 + sigma_0 * k);

	w[2][M_w] = ((2 - J[0] * muSq + (Iterm - 2) * ((2 * sigma_1 * k) / hSq)) * w[1][M_w]
		+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * w[1][M_w - 1]
		- muSq * w[1][M_w - 2]
		+ (((2 * sigma_1 * k) / hSq) - J[1] * muSq) * u[1][M_u]
		+ (J[3] * ((2 * sigma_1 * k) / hSq) - J[2] * muSq) * u[1][M_u - 1]
		- (J[3] * muSq) * u[1][M_u - 2]
		+ (-1 + sigma_0 * k - (Iterm - 2) * ((2 * sigma_1 * k) / hSq)) * w[0][M_w]
		- ((2 * sigma_1 * k) / hSq) * (w[0][M_w - 1] + u[0][M_u])
		- (J[3] * ((2 * sigma_1 * k) / hSq)) * u[0][M_u - 1]
		) / (1 + sigma_0 * k);

	//Free Boundary
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
	jassert(abs(N - Nprev) <= 1);

	if (Nprev == N) //Same amount of grid points
		return;

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
	//This will run for the 3 state vectors of the given grid
	for (int i = 0; i < grid.size(); i++)
	{
		grid[i].pop_back();
	}
	//TODO Displacement correction
}

void Tine::addGridpoint(std::vector<std::vector<float>>& main, std::vector<std::vector<float>>& sec)
{
	interp[0] = -alpha * (alpha + 1.0f) / ((alpha + 2.0f) * (alpha + 3.0f));
	interp[1] = 2.0f * alpha / (alpha + 2.0f);
	interp[2] = 2.0f / (alpha + 2.0f);
	interp[3] = -2.0f * alpha / ((alpha + 3.0f) * (alpha + 2.0f));

	for (int i = 0; i < main.size(); i++)
	{
		int M = main[i].size() - 1;
		int S = sec[i].size() - 1;
		float v[4] = { main[i][M - 1], main[i][M], sec[i][S], sec[i][S - 1] }; //Build v vector

		float point = 0.0f;
		for (int j = 0; j < 4; j++)
			point += interp[j] * v[j];

		main[i].push_back(point);
	}
}

void Tine::calculateInterpolation()
{
	Iterm = (alpha - 1.0f) / (alpha + 1.0f);
	float itermSq = Iterm * Iterm;
	J[0] = itermSq - 4.0f * Iterm + 6.0f;
	J[1] = Iterm - 4.0f;
	J[2] = -itermSq + 4.0f * Iterm + 1.0f;
	J[3] = -Iterm;
}

float Tine::limit(float sample)
{
	if (sample > 1.0f)
	{
		DBG("Sample Exceeded 1.0");
		return 1.0f;
	}
	else if (sample < -1.0f)
	{
		DBG("Sample Exceeded -1.0");
		return -1.0f;
	}
	else
		return sample;
}

float Tine::calculateLength(float freq)
{
	return sqrt(((1.426f * juce::MathConstants<float>::pi * K * kappa_1) / freq) / 8);
}