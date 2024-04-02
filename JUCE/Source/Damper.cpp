/*
  ==============================================================================

    Damper.cpp
    Created: 19 Mar 2024 3:00:26pm
    Author:  Sonderbo

  ==============================================================================
*/

#include "Damper.h"

void Damper::prepareToPlay(double sampleRate, int N)
{

	psiStates = std::vector<std::vector<float>>(3, std::vector<float>(N, 0));
	psi.resize(3, nullptr);

	for (int i = 0; i < 3; i++)
	{
		psi[i] = &psiStates[i][0];
	}
}

void Damper::calculateDamping(const float* tine_u)
{

}

