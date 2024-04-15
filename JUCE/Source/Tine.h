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

class Tine : public juce::MPESynthesiserVoice
{
public:

	Tine();
	//Overrides ====================================================================
	void noteStarted() override;
	void noteStopped(bool allowTailOff) override;
	void notePitchbendChanged() override;
	void notePressureChanged() override;
	void noteTimbreChanged() override;
	void noteKeyStateChanged() override;
	void renderNextBlock(juce::AudioBuffer<float>& outputBuffer,
		int startSample,
		int numSamples) override;
	//===============================================================================
	void prepareToPlay(double sampleRate);
	float processSample();

private:

	//Util
	double sampleRate; //Sample rate
	float k; //Sampling period
	float kSq; //Sampling period squared
	int outLoc;
	bool isPrepared;

	//Material properties
	float L; //Length - From Tine Length Config
	float rho = 7850.0f; //Material density
	float r = 1.524f * powf(10, -3); //Radius
	float A = juce::MathConstants<float>::pi * powf(r, 2); //Cross-sectional area
	float E = 2.0f * powf(10, 11); //Young's modulus
	float I = juce::MathConstants<float>::pi * powf(r, 4) / 4.0f; //Inertia
	float K = r/2.0f;
	float kappa_1 = sqrtf(E / rho);

	//Coefficients
	float kappa = sqrtf((E*I)/(rho*A)); //Stiffness coefficient
	float mu; //Coefficient
	float muSq; //Coefficient

	//Damping
	float sigma_0 = 0.0005f; //Freq. independent damping
	float sigma_1 = 0.001f; //Freq. dependent damping

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
	float Iterm; //Interpolation term
	float J[4]; //
	float interp[4]; //Cubic inpterpolation coefficient matrix
	
	int M_u;
	int M_w;

	//Activity states
	bool isStopped;

	//Excitation
	Hammer hammer;
	std::vector<float> h_contact; //Hammer contact distribution
	float h_ratio; //Hammer-Tine mass ratio

	//Material coefficients
	float C_0;
	float C_1;
	float B_0;
	float S;
	float D;

	float A_0;
	float A_1;
	float A_2;
	float A_3;
	float A_4;
	float A_5;

	float E_1;

	float F_0;
	float F_1;
	float F_2;


	//==============================================================================
	void calculateScheme();
	void updateStates();
	void resetScheme();
	void updateGridpoints(int Nnew);
	void addGridpoint(std::vector<std::vector<float>>& main, std::vector<std::vector<float>>& sec);
	void removeGridpoint(std::vector<std::vector<float>>& grid);
	void calculateInterpolation();
	void prepareGrid(float freq);
	float calculateLength(float freq);
	float limit(float sample);
	bool isNoteValid();
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Tine)
};