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
    gain.setGainDecibels(12.0f);
    gain.setRampDurationSeconds(20.0); //20 seconds?

    symmetry_gain.prepare(spec);
    symmetry_gain.setGainDecibels(15.0f);
    symmetry_gain.setRampDurationSeconds(20.0);

    lowpass.prepare(spec);
    highpass.prepare(spec);
}

float Pickup::processSample(float sample)
{
    //Resonant Lowpass
    //TODO - Configure SVF filter
    //sample = lowpass.processSample(sample);

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
    sample = sample - (powf(sample, 3) / 3);

    //Resonant Highpass
    //TODO - Configure SVF filter
    //sample = highpass.processSample(sample);

    return sample;
}

void Pickup::reset()
{
    gain.reset();
    symmetry_gain.reset();
    lowpass.reset();
    highpass.reset();
}
