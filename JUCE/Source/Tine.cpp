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
	pickup.prepareToPlay(sampleRate);
	this->sampleRate = sampleRate;
	tine_gain.reset(sampleRate, config::tine::smoothDuration);
	tine_gain.setCurrentAndTargetValue(config::tine::tineGain);
}

void Tine::prepareGrid(float freq)
{
	oversampling = calculateOversampling();

	k = 1.0f / (sampleRate * oversampling);
	kSq = k * k;

	calculateGridSpacing();

	hSq = h * h;

	mu = (kappa * k) / hSq;
	muSq = mu * mu;

	//Calculate Coefficients
	S = (2.0f * sigma_1 * k) / hSq;
	S_2 = 2.0f * S;

	C_0 = (2.0f - 6.0f * muSq - S_2);
	C_1 = 4.0f * muSq + S;
	B_0 = (-1.0f + sigma_0 * k + S_2);
	D = 1.0f / (1.0f + sigma_0 * k);

	F_0 = (2.0f - 5.0f * muSq - S_2);
	F_1 = (S + 2.0f * muSq);
	F_2 = (2.0f - 2.0f * muSq - S_2);

	L = calculateLength(freq);

	N = static_cast<int>(L / h);
	N_frac = L / h;
	alpha = N_frac - N;

	//Ensure minimum grid points
	jassert(N >= 10);

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
	
	int contact_point = static_cast<int>(0.33f * M_u);

	if (contact_point < 4)
		contact_point = 4;

	h_contact[contact_point] = 1;

	hammer.prepareToPlay(k, h_contact, M_u);
	h_ratio = hammer.getMass() / (rho * A * L);

	//Pre-calculate hammer contact
	h_contact[contact_point] = 1 * kSq * h_ratio; 

	d_loc = static_cast<int>(M_w * 0.5f);
}

void Tine::noteStarted()
{
	jassert(currentlyPlayingNote.isValid());
	jassert(currentlyPlayingNote.keyState == juce::MPENote::keyDown
		|| currentlyPlayingNote.keyState == juce::MPENote::keyDownAndSustained);

	if (isNoteValid() == false)
		return;

	float freq = currentlyPlayingNote.getFrequencyInHertz();
	float velocity = currentlyPlayingNote.noteOnVelocity.asUnsignedFloat() * config::mpe::maxInputVelocity;

	prepareGrid(freq);

	pickup.reset();
	pickup.setFreq(freq);

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
	if (isNoteValid() == false)
		return;

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

	pickup.setFreq(freq);
}

float Tine::processSample()
{
	if (isNoteValid() == false)
		return 0;

	float sample;

	for (size_t i = 0; i < oversampling; i++)
	{
		calculateScheme();
	}

	sample = w[1][0]; //Read output at tip

	//TODO pickup filtering goes here...

	sample = sample * tine_gain.getNextValue();

	if (pickup_bypass == false)
	{
		sample = pickup.processSample(sample);
	}

	sample = limit(sample);

	return sample;
}

void Tine::setParameters(juce::NamedValueSet paramValueSet)
{
	//Set pickup values
	pickup.setParameters(paramValueSet);

	using namespace config::parameter;

	tine_gain.setTargetValue(paramValueSet[id_tine_gain]);

	pickup_bypass = static_cast<bool>(paramValueSet[id_pickup_bypass]);
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
		u[2][l] = (C_0 * u[1][l]
			+ C_1 * (u[1][l + 1] + u[1][l - 1])
			- muSq * (u[1][l + 2] + u[1][l - 2])
			+ B_0 * u[0][l]
			- S * (u[0][l + 1] + u[0][l - 1])
			+ (h_contact[l] * force))
			* D;

		w[2][l] = (C_0 * w[1][l]
			+ C_1 * (w[1][l + 1] + w[1][l - 1])
			- muSq * (w[1][l + 2] + w[1][l - 2])
			+ B_0 * w[0][l]
			- S * (w[0][l + 1] + w[0][l - 1]))
			* D;
	}

	if (M_u > M_w) //In case of uneven grid lengths - U will always be greater
	{
		int l = M_u - 2;

		u[2][l] = (C_0 * u[1][l]
			+ C_1 * (u[1][l + 1] + u[1][l - 1])
			- muSq * (u[1][l + 2] + u[1][l - 2])
			+ B_0 * u[0][l]
			- S * (u[0][l + 1] + u[0][l - 1])
			+ (h_ratio * h_contact[l] * force * kSq))
			* D;
	}

	//Inner boundaries
	u[2][M_u - 1] = (C_0 * u[1][M_u - 1]
		+ E_1 * u[1][M_u]
		+ C_1 * u[1][M_u - 2]
		- muSq * u[1][M_u - 3]
		- muSq * w[1][M_w]
		- A_3 * w[1][M_w - 1]
		+ B_0 * u[0][M_u - 1]
		- S * (u[0][M_u - 2] + u[0][M_u]))
		* D;

	u[2][M_u] = (A_0 * u[1][M_u]
		+ C_1 * u[1][M_u - 1]
		- muSq * u[1][M_u - 2]
		+ A_1 * w[1][M_w]
		+ A_2 * w[1][M_w - 1]
		- A_3 * w[1][M_w - 2]
		+ A_4 * u[0][M_u]
		- S * (u[0][M_u - 1] + w[0][M_w])
		- A_5 * w[0][M_w - 1])
		* D;

	w[2][M_w - 1] = (C_0 * w[1][M_w - 1]
		+ E_1 * w[1][M_w]
		+ C_1 * w[1][M_w - 2]
		- muSq * w[1][M_w - 3]
		- muSq * u[1][M_u]
		- A_3 * u[1][M_u - 1]
		+ B_0 * w[0][M_w - 1]
		- S * (w[0][M_w - 2] + w[0][M_w]))
		* D;

	w[2][M_w] = (A_0 * w[1][M_w]
		+ C_1 * w[1][M_w - 1]
		- muSq * w[1][M_w - 2]
		+ A_1 * u[1][M_u]
		+ A_2 * u[1][M_u - 1]
		- A_3 * u[1][M_u - 2]
		+ A_4 * w[0][M_w]
		- S * (w[0][M_w - 1] + u[0][M_u])
		- A_5 * u[0][M_u - 1])
		* D;

	//Free Boundary
	w[2][1] = (F_0 * w[1][1]
		+ F_1 * w[1][0]
		+ C_1 * w[1][2]
		- muSq * w[1][3]
		+ B_0 * w[0][1]
		- S * (w[0][2] + w[0][0]))
		* D;

	w[2][0] = (F_2 * w[1][0]
		+ C_1 * w[1][1]
		- muSq * 2 * w[1][2]
		+ B_0 * w[0][0]
		- S * (2 * w[0][1]))
		* D;

	if (isStopped)
		applyDamping();

	updateStates();
}

void Tine::applyDamping()
{
	float d_eta = w[1][d_loc];
	float d_etaPrev = w[0][d_loc];
	float wStar = w[2][d_loc];


	float rPlus = K1 / 4 + K3 * (d_eta * d_eta) / 2 + R / (2 * k);
	float rMinus = K1 / 4 + K3 * (d_eta * d_eta) / 2 - R / (2 * k);

	float d_force = ((wStar + (K1 / (2 * rPlus) * d_eta) + (rMinus / (rPlus)*d_etaPrev)) 
		/ (1 / rPlus + (1 / h * kSq) / (rho * A * (1 + sigma_0 * k)))) * 0.5;

	w[2][d_loc] = wStar - (1 / h) * ((kSq * d_force) / (rho * A * (1 + sigma_0 * k)));
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

	if (abs(N - Nprev) > 1)
	{
		if (Nprev < N) //Addition
		{
			updateGridpoints(Nprev + 1);
		}
		else
		{
			updateGridpoints(Nprev - 1);
		}
	}

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

	//Recalculate Coefficients with interpolation values
	E_1 = (-J[1] * muSq + (2 * sigma_1 * k) / hSq);
	A_0 = (2 - J[0] * muSq + (Iterm - 2) * S);
	A_1 = (S - J[1] * muSq);
	A_2 = (J[3] * S - J[2] * muSq);
	A_3 = (J[3] * muSq);
	A_4 = (-1 + sigma_0 * k - (Iterm - 2) * S);
	A_5 = (J[3] * S);
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
	return sqrt((1.426f * juce::MathConstants<float>::pi * K * kappa_1 / freq) / 8);
}

bool Tine::isNoteValid()
{
	auto note = getCurrentlyPlayingNote().initialNote;
	return  note >= config::mpe::minNote && note <= config::mpe::maxNote;
}

void inline Tine::calculateGridSpacing()
{
	h = sqrtf((4.0f * sigma_1 * k + sqrtf(pow((4.0f * sigma_1 * k), 2.0f) + 16.0f * kappa * kappa * kSq)) / 2.0f);
}

int Tine::calculateOversampling()
{
	auto note = getCurrentlyPlayingNote().initialNote;

	if (note < 40)
	{
		return 4;
	}
	else if (note < 52)
	{
		return 8;
	}
	else if (note < 64)
	{
		return 16;
	}
	else
	{
		return 32;
	}
}