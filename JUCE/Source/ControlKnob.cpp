/*
  ==============================================================================

    ControlKnob.cpp
    Created: 26 Apr 2024 2:38:16pm
    Author:  Sonderbo

  ==============================================================================
*/

#include <JuceHeader.h>
#include "ControlKnob.h"

//==============================================================================
ControlKnob::ControlKnob()
{
    setSize(60, 110);
}

void ControlKnob::configure(EPianoAudioProcessor& p, juce::String identifier, juce::String name)
{
    knob.setSliderStyle(juce::Slider::SliderStyle::Rotary);
    knob.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 50, 25);

    knobLabel.setColour(juce::Label::ColourIds::textColourId, juce::Colours::white);
    knobLabel.setFont(15.0f);
    knobLabel.setJustificationType(juce::Justification::centred);
    knobLabel.setText(name, juce::dontSendNotification);

    attachment.reset(new juce::AudioProcessorValueTreeState::SliderAttachment(p.apvts, identifier, knob));

    addAndMakeVisible(knobLabel);
    addAndMakeVisible(knob);
}

ControlKnob::~ControlKnob()
{
}

void ControlKnob::paint (juce::Graphics& g)
{
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));
}

void ControlKnob::resized()
{
    juce::FlexBox fb;
    fb.flexDirection = juce::FlexBox::Direction::column;

    fb.items.add(juce::FlexItem(knobLabel)
        .withMinWidth(50.0f)
        .withMinHeight(10.0f)
        .withFlex(1));

    fb.items.add(juce::FlexItem(knob)
        .withMinWidth(50.0f)
        .withMinHeight(50.0f)
        .withFlex(1));

    fb.performLayout(getLocalBounds().toFloat());
}
