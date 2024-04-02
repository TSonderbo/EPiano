/*
  ==============================================================================

	Tine.h
	Created: 8 Mar 2024 10:33:32am
	Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Hammer.h"
#include "TineLengths.h"
#include "Definitions.h"

class Tine
{
public:

	Tine();
	//==============================================================================
	void prepareToPlay(double sampleRate, int tineNumber);
	void startNote(float velocity);
	void stopNote();
	float processSample();
	bool getIsActive();

private:

	//Util
	double sampleRate; //Sample rate
	float k; //Sampling period
	float kSq; //Sampling period squared
	int outLoc;
	int inactiveTimer;
	int inactiveTimerItterator;

	//Material properties
	float L; //Length - From Tine Length Config
	float rho = 7850.0f; //Material density
	float r = 1.524 * pow(10, -3); //Radius
	float A = juce::MathConstants<float>::pi * pow(r, 2); //Cross-sectional area
	float E = 2 * pow(10, 11); //Young's modulus
	float I = juce::MathConstants<float>::pi * pow(r, 4) / 4; //Inertia

	//Coefficients
	float kappa = sqrtf((E*I)/(rho*A)); //Stiffness coefficient
	float mu; //Coefficient
	float muSq; //Coefficient

	//Damping
	float sigma_0 = 0.005; //Freq. dependent damping
	float sigma_1 = 0.0001; //Freq. independent damping

	//Grid
	int N; //Number of grid intervals
	float h; //Grid spacing
	float hSq; //Grid spacing squared

	//States
	std::vector<std::vector<float>> uStates; //3 by N+1 state matrix
	std::vector<float*> u; //uStates pointer
	bool isActive;
	bool isStopped;

	//Excitation
	Hammer hammer; //Hammer
	std::vector<float> h_contact; //Hammer contact distribution
	float h_ratio;

	//==============================================================================
	void calculateScheme();
	void updateStates();
	void resetScheme();
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Tine)
};