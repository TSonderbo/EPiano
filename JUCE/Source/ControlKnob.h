/*
  ==============================================================================

    ControlKnob.h
    Created: 26 Apr 2024 2:38:16pm
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/*
*/
class ControlKnob  : public juce::Component
{
public:
    ControlKnob();
    ~ControlKnob() override;

    void paint (juce::Graphics&) override;
    void resized() override;

    void configure(EPianoAudioProcessor& p, juce::String identifier, juce::String name);

    juce::Slider knob;
    juce::Label knobLabel;

private:

    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> attachment;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ControlKnob)
};
