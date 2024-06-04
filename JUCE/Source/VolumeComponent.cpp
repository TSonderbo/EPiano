/*
  ==============================================================================

    VolumeComponent.cpp
    Created: 26 Apr 2024 9:18:50am
    Author:  Sonderbo

  ==============================================================================
*/

#include <JuceHeader.h>
#include "VolumeComponent.h"

//==============================================================================
VolumeComponent::VolumeComponent(EPianoAudioProcessor& p) : processor(p)
{
    // In your constructor, you should add any child components, and
    // initialise any special settings that your component needs.
    using namespace config::parameter;

    masterKnob.configure(p, id_master_volume, name_master_volume);
    tineGainKnob.configure(p, id_tine_gain, name_tine_gain);

    masterKnob.knob.setSliderStyle(juce::Slider::SliderStyle::LinearBarVertical);
    tineGainKnob.knob.setSliderStyle(juce::Slider::SliderStyle::LinearBarVertical);

    

    addAndMakeVisible(masterKnob);
    addAndMakeVisible(tineGainKnob);
}

VolumeComponent::~VolumeComponent()
{
}

void VolumeComponent::paint (juce::Graphics& g)
{
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));   // clear the background

    g.setColour(juce::Colours::grey);

    auto bounds = getLocalBounds();
    bounds.reduce(4.0f, 4.0f);

    g.drawRect(bounds, 2); // draw an outline around the component

    g.setColour(juce::Colours::white);
    g.setFont(15.0f);
    g.drawFittedText("Tone Controls", 8.0f, 8.0f, 100.0f, 15.0f, juce::Justification::centred, 1);
}

void VolumeComponent::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..

    juce::FlexBox fb;

    fb.items.add(juce::FlexItem(masterKnob)
        .withMinWidth(50.0f)
        .withMinHeight(100.0f)
        .withFlex(1));

    fb.items.add(juce::FlexItem(tineGainKnob)
        .withMinWidth(50.0f)
        .withMinHeight(100.0f)
        .withFlex(1));

    auto rect = getLocalBounds();
    rect.reduce(8, 8);
    fb.performLayout(rect);
}
