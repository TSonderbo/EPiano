/*
  ==============================================================================

    ExcitationComponent.h
    Created: 26 Apr 2024 9:18:24am
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class ExcitationComponent  : public juce::Component
{
public:
    ExcitationComponent();
    ~ExcitationComponent() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ExcitationComponent)
};
