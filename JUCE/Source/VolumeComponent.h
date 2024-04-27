/*
  ==============================================================================

    VolumeComponent.h
    Created: 26 Apr 2024 9:18:50am
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class VolumeComponent  : public juce::Component
{
public:
    VolumeComponent();
    ~VolumeComponent() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeComponent)
};
