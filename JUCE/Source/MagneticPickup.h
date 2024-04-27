/*
  ==============================================================================

    MagneticPickup.h
    Created: 8 Mar 2024 11:44:28am
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Configuration.h"

class Pickup
{
public:

    void prepareToPlay(double sampleRate);
    float processSample(float sample);
    void reset();
    void setFreq(float freq);
    void setParameters(juce::NamedValueSet paramValueSet);

private:

    double sampleRate;

    bool bypass;

    juce::dsp::Gain<float> gain;
    juce::dsp::Gain<float> symmetry_gain;

    juce::dsp::StateVariableTPTFilter<float> lowpass;
    juce::dsp::StateVariableTPTFilter<float> highpass;

    //Smoothed values
    juce::SmoothedValue<float, juce::ValueSmoothingTypes::Multiplicative> smooLpCutoff;
    juce::SmoothedValue<float> smooLpResonance;
    juce::SmoothedValue<float> smooHpResonance;

    void smoothParametersChanges();
};