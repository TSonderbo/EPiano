/*
  ==============================================================================

    MagneticPickup.h
    Created: 8 Mar 2024 11:44:28am
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

class Pickup
{
public:

    void prepareToPlay(double sampleRate);
    float processSample(float sample);
    void reset();

private:

    double sampleRate;

    juce::dsp::Gain<float> gain;
    juce::dsp::Gain<float> symmetry_gain;

    juce::dsp::StateVariableFilter::Filter<float> lowpass;
    juce::dsp::StateVariableFilter::Filter<float> highpass;
};