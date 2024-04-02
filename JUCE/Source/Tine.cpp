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
	k = 1 / (sampleRate * OVERSAMPLING);
	kSq = k * k;
	this->sampleRate = sampleRate * OVERSAMPLING;

	h = sqrtf((4 * sigma_1 * k + sqrtf(pow((4 * sigma_1 * k), 2) + 16 * kappa * kappa * k * k)) / 2);

	N = static_cast<int>(L / h);

	h = L / N;
	hSq = h * h;

	mu = (kappa * k) / hSq;
	muSq = mu * mu;

	uStates = std::vector<std::vector<float>>(3, std::vector<float>(N + 1, 0));
	u.resize(3, nullptr);

	for (int i = 0; i < 3; i++)
	{
		u[i] = &uStates[i][0];
	}

	h_contact.resize(N + 1, 0.0f);

	//TODO Calculate h_contact - static for now
	h_contact[static_cast<int>(0.33 * N)] = 1;
	
	hammer.prepareToPlay(sampleRate, tineNumber, h_contact, N);
	h_ratio = hammer.getMass() / (rho * A * L);
	//Read ouput at tip
	outLoc = N;

	inactiveTimer = TINE_INACTIVE_TIMER * sampleRate / OVERSAMPLING;
}

void Tine::startNote(float velocity)
{
	hammer.beginHammer(velocity);
	isActive = true;
	inactiveTimerItterator = inactiveTimer; //Inactive timer in samples
	isStopped = false;
}

void Tine::stopNote()
{
	isStopped = true;
	// apply damper
}

float Tine::processSample()
{
	if (isActive == false)
		return 0.0f;
	if (isStopped == true)
	{
		inactiveTimerItterator = inactiveTimerItterator - 1;
		if (inactiveTimerItterator <= 0)
		{
			isActive = false;
			resetScheme();
		}
			
	}

	float sample;

	for (size_t i = 0; i < OVERSAMPLING; i++)
	{
		calculateScheme();
	}

	sample = u[1][outLoc];

	float gain = 1.0f;

	if (isStopped == true)
		gain = (1 / inactiveTimer) * inactiveTimerItterator;

	return sample * gain;
}

void Tine::calculateScheme()
{

	float force = hammer.calculateForce(u[1]);

	//Update loop
	for (int l = 2; l < N - 1; l++)
	{
		u[2][l] = ((2 - 6 * muSq - ((4 * sigma_1 * k) / hSq)) * u[1][l]
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (u[1][l + 1] + u[1][l - 1])
			- muSq * (u[1][l + 2] + u[1][l - 2]) + (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * u[0][l]
			- ((2 * sigma_1 * k) / hSq) * (u[0][l + 1] + u[0][l - 1]))
			+ (h_ratio * h_contact[l] * force * kSq)
			/ (1 + sigma_0 * k);
	}

	//Free boundary update
	u[2][N - 1] = ((2 - 5 * muSq - ((4 * sigma_1 * k) / hSq)) * u[1][N - 1]
		+ ((2 * sigma_1 * k) / (hSq)+2 * muSq) * u[1][N]
		+ (4 * muSq + (2 * sigma_1 * k) / (hSq)) * u[1][N - 2]
		- muSq * u[1][N - 3]
		+ (-1 + sigma_0 * k + (4 * sigma_1 * k) / hSq) * u[0][N - 1]
		- ((2 * sigma_1 * k) / (hSq)) * (u[0][N - 2] + u[0][N]));

	u[2][N] = ((2 - 2 * muSq - ((4 * sigma_1 * k) / hSq)) * u[1][N]
		+ (4 * muSq + ((4 * sigma_1 * k) / hSq)) * u[1][N - 1]
		- muSq * 2 * u[1][N - 2]
		+ (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * u[0][N]
		- ((2 * sigma_1 * k) / hSq) * (2 * u[0][N - 1]))
		/ (1 + sigma_0 * k);

	updateStates();
}

void Tine::updateStates()
{
	auto uPrev = u[0];
	u[0] = u[1];
	u[1] = u[2];
	u[2] = uPrev;
}

bool Tine::getIsActive()
{
	return isActive;
}

void Tine::resetScheme()
{
	uStates = std::vector<std::vector<float>>(3, std::vector<float>(N + 1, 0));

	for (int i = 0; i < 3; i++)
	{
		u[i] = &uStates[i][0];
	}
}