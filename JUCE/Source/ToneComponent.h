/*
  ==============================================================================

    ToneComponent.h
    Created: 26 Apr 2024 9:18:01am
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"
#include "ControlKnob.h"
#include "Configuration.h"

//==============================================================================
/*
*/
class ToneComponent  : public juce::Component
{
public:
    ToneComponent(EPianoAudioProcessor& p);
    ~ToneComponent() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:

    ControlKnob lowpassCutoffKnob;
    ControlKnob lowpassResonanceKnob;
    ControlKnob highpassResonanceKnob;

    ControlKnob gainKnob;
    ControlKnob symmetryKnob;

    juce::ToggleButton pickupToggle;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> attachment;

    EPianoAudioProcessor& processor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ToneComponent)
};
