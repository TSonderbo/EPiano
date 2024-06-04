/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"
#include "ScopeComponent.h"
#include "ToneComponent.h"
#include "VolumeComponent.h"
//==============================================================================
/**
*/
class EPianoAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    EPianoAudioProcessorEditor (EPianoAudioProcessor&);
    ~EPianoAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    EPianoAudioProcessor& audioProcessor;

    ScopeComponent oscilloscope;
    ToneComponent toneComponent;
    VolumeComponent volumeComponent;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (EPianoAudioProcessorEditor)
};
