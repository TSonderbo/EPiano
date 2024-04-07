/*
  ==============================================================================

    EPianoManager.h
    Created: 8 Mar 2024 11:46:19am
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Configuration.h"
#include "Tine.h"

class EPianoManager : public juce::MPESynthesiserBase
{
public:
    EPianoManager();

    void prepareToPlay(double sampleRate);

    void noteAdded(juce::MPENote newNote) override;
    //void notePressureChanged(juce::MPENote changedNote) override;
    void notePitchbendChanged(juce::MPENote changedNote) override;
    //void noteTimbreChanged(juce::MPENote changedNote) override;
    //void noteKeyStateChanged(juce::MPENote changedNote) override;
    void noteReleased(juce::MPENote finishedNote) override;
    //void zoneLayoutChanged() override;

    void renderNextSubBlock(juce::AudioBuffer<float>& outputAudio,
        int startSample,
        int numSamples) override;

    void setMaximumInputVelocity(float newVelocityMax);

private:

    std::vector<Tine> tines;

    float v_max = 10.0f; //Maximum input velocity

    float renderNextSample();
    float limit(float sample);

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(EPianoManager)
};