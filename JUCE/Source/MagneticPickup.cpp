/*
  ==============================================================================

	MagneticPickup.cpp
	Created: 8 Mar 2024 11:44:28am
	Author:  Sonderbo

  ==============================================================================
*/

#include "MagneticPickup.h"

void Pickup::prepareToPlay(double sampleRate)
{
	this->sampleRate = sampleRate;

	juce::dsp::ProcessSpec spec;
	spec.sampleRate = sampleRate;

	gain.prepare(spec);
	gain.setGainDecibels(config::pickup::gain);
	gain.setRampDurationSeconds(config::pickup::smoothDuration);

	symmetry_gain.prepare(spec);
	symmetry_gain.setGainDecibels(config::pickup::symmetry);
	symmetry_gain.setRampDurationSeconds(config::pickup::smoothDuration);

	lowpass.prepare(spec);
	highpass.prepare(spec);

	highpass.setType(juce::dsp::StateVariableTPTFilterType::highpass);
	highpass.setResonance(config::pickup::resonance);

	lowpass.setType(juce::dsp::StateVariableTPTFilterType::lowpass);
	lowpass.setCutoffFrequency(config::pickup::lpCutoff);
	lowpass.setResonance(config::pickup::resonance);

	smooLpCutoff.reset(sampleRate, config::pickup::smoothDuration);
	smooLpCutoff.setCurrentAndTargetValue(config::pickup::lpCutoff);

	smooLpResonance.reset(sampleRate, config::pickup::smoothDuration);
	smooLpResonance.setCurrentAndTargetValue(config::pickup::resonance);

	smooHpResonance.reset(sampleRate, config::pickup::smoothDuration);
	smooHpResonance.setCurrentAndTargetValue(config::pickup::resonance);

	bypass = false;
}

float Pickup::processSample(float sample)
{
	if (bypass)
		return sample;

	smoothParametersChanges();

	//Resonant Lowpass
	sample = lowpass.processSample(0, sample);
	lowpass.snapToZero();

	//Apply gain
	sample = gain.processSample(sample);

	//Apply soft-clipper
	sample = juce::dsp::FastMathApproximations::tanh(sample);

	float path1 = sample;
	float path2;

	//Path 1
	path1 = symmetry_gain.processSample(path1);

	path1 = powf(2, path1);
	path1 -= 1; //Scale signal

	//Path 2
	path2 = powf(2, symmetry_gain.getGainLinear());

	path2 = 1 / path2;

	//Combine
	sample = path1 * path2;

	//Buzz
	sample = sample - (powf(sample, 3.0f) / 3.0f);

	//Resonant Highpass
	//TODO - Configure SVF filter
	sample = highpass.processSample(0, sample);
	highpass.snapToZero();

	return sample;
}

void Pickup::reset()
{
	gain.reset();
	symmetry_gain.reset();
	lowpass.reset();
	highpass.reset();
}

void Pickup::setFreq(float freq)
{
	highpass.setCutoffFrequency(freq);
}

void Pickup::setParameters(juce::NamedValueSet paramValueSet)
{
	using namespace config::parameter;

	float lp_cutoff = paramValueSet[id_pickup_lowpass_cutoff];
	float lp_res = paramValueSet[id_pickup_lowpass_resonance];

	float hp_res = paramValueSet[id_pickup_highpass_resonance];

	float g = paramValueSet[id_pickup_gain];
	float sym = paramValueSet[id_pickup_symmetry];

	smooLpCutoff.setTargetValue(lp_cutoff);
	smooLpResonance.setTargetValue(lp_res);
	smooHpResonance.setTargetValue(hp_res);

	gain.setGainDecibels(g);
	symmetry_gain.setGainDecibels(sym);

	bypass = static_cast<bool>(paramValueSet[id_pickup_bypass]);
}

void Pickup::smoothParametersChanges()
{
	if ((smooLpCutoff.isSmoothing() || smooLpResonance.isSmoothing() || smooHpResonance.isSmoothing()) == false) //If no parameters is currently changing
		return;

	float lp_cutoff = smooLpCutoff.getNextValue();
	float lp_res = smooLpResonance.getNextValue();
	float hp_res = smooHpResonance.getNextValue();

	lowpass.setCutoffFrequency(lp_cutoff);
	lowpass.setResonance(lp_res);
	highpass.setResonance(hp_res);
}