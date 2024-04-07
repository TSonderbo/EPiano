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
#include "Configuration.h"

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
	int inactiveTimerIterator;

	//Material properties
	float L; //Length - From Tine Length Config
	float rho = 7850.0f; //Material density
	float r = 1.524f * pow(10, -3); //Radius
	float A = juce::MathConstants<float>::pi * pow(r, 2); //Cross-sectional area
	float E = 2.0f * pow(10, 11); //Young's modulus
	float I = juce::MathConstants<float>::pi * pow(r, 4) / 4.0f; //Inertia

	//Coefficients
	float kappa = sqrtf((E*I)/(rho*A)); //Stiffness coefficient
	float mu; //Coefficient
	float muSq; //Coefficient

	//Damping
	float sigma_0 = 0.05f; //Freq. dependent damping
	float sigma_1 = 0.0001f; //Freq. independent damping

	//Grid
	int N; //Number of grid intervals
	float h; //Grid spacing
	float hSq; //Grid spacing squared

	//States
	std::vector<std::vector<float>> uStates; //3 by Mu+1 state matrix
	std::vector<std::vector<float>> wStates; //3 by Mw+1 state matrix

	std::vector<std::vector<float>::iterator> u; //uStates iterators
	std::vector<std::vector<float>::iterator> w; //wStates iterators

	//Dynamic Grid
	float N_frac; //Fractional grid intervals
	float alpha;
	float interp[4]; //Cubic inpterpolation coefficient matrix
	float interpFlip[4]; //Flipped cubic inpterpolation coefficient matrix
	int M_u;
	int M_w;


	//Activity states
	bool isActive;
	bool isStopped;

	//Excitation
	Hammer hammer;
	std::vector<float> h_contact; //Hammer contact distribution
	float h_ratio; //Hammer-Tine mass ratio

	//==============================================================================
	void calculateScheme();
	void updateStates();
	void resetScheme();
	void updateGridpoints(int Nnew);
	void addGridpoint(std::vector<std::vector<float>>& main, std::vector<std::vector<float>>& sec);
	void removeGridpoint(std::vector<std::vector<float>>& grid);
	void calculateInterpolation();

	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Tine)
};