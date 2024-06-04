/*
  ==============================================================================

	ToneComponent.cpp
	Created: 26 Apr 2024 9:18:01am
	Author:  Sonderbo

  ==============================================================================
*/

#include <JuceHeader.h>
#include "ToneComponent.h"

using namespace config::parameter;

//==============================================================================
ToneComponent::ToneComponent(EPianoAudioProcessor& p) : processor(p)
{
	lowpassCutoffKnob.configure(p, id_pickup_lowpass_cutoff, name_pickup_lowpass_cutoff);
	lowpassResonanceKnob.configure(p, id_pickup_lowpass_resonance, name_pickup_lowpass_resonance);
	highpassResonanceKnob.configure(p, id_pickup_highpass_resonance, name_pickup_highpass_resonance);
	gainKnob.configure(p, id_pickup_gain, name_pickup_gain);
	symmetryKnob.configure(p, id_pickup_symmetry, name_pickup_symmetry);

	pickupToggle.setButtonText("Bypass");
	attachment.reset(new juce::AudioProcessorValueTreeState::ButtonAttachment(p.apvts, id_pickup_bypass, pickupToggle));

	addAndMakeVisible(lowpassCutoffKnob);
	addAndMakeVisible(lowpassResonanceKnob);
	addAndMakeVisible(highpassResonanceKnob);
	addAndMakeVisible(gainKnob);
	addAndMakeVisible(symmetryKnob);
	addAndMakeVisible(pickupToggle);
}

ToneComponent::~ToneComponent()
{
}

void ToneComponent::paint(juce::Graphics& g)
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

void ToneComponent::resized()
{
	pickupToggle.setBounds(12, 30, 100, 20);

	juce::FlexBox fb;
	fb.flexDirection = juce::FlexBox::Direction::column;

	juce::FlexBox filtersFB;
	filtersFB.flexWrap = juce::FlexBox::Wrap::noWrap;
	filtersFB.justifyContent = juce::FlexBox::JustifyContent::spaceBetween;
	filtersFB.alignContent = juce::FlexBox::AlignContent::center;
	filtersFB.alignItems = juce::FlexBox::AlignItems::center;
	filtersFB.flexDirection = juce::FlexBox::Direction::row;

	filtersFB.items.add(juce::FlexItem(lowpassCutoffKnob)
		.withMinWidth(50.0f)
		.withMinHeight(100.0f)
		.withFlex(1));

	filtersFB.items.add(juce::FlexItem(lowpassResonanceKnob)
		.withMinWidth(50.0f)
		.withMinHeight(100.0f)
		.withFlex(1));

	filtersFB.items.add(juce::FlexItem(highpassResonanceKnob)
		.withMinWidth(50.0f)
		.withMinHeight(100.0f)
		.withFlex(1));

	juce::FlexBox gainsFB;
	gainsFB.flexWrap = juce::FlexBox::Wrap::noWrap;
	gainsFB.justifyContent = juce::FlexBox::JustifyContent::spaceBetween;
	gainsFB.alignContent = juce::FlexBox::AlignContent::center;
	gainsFB.alignItems = juce::FlexBox::AlignItems::center;
	gainsFB.flexDirection = juce::FlexBox::Direction::row;

	gainsFB.items.add(juce::FlexItem(gainKnob)
		.withMinWidth(50.0f)
		.withMinHeight(100.0f)
		.withFlex(1));

	gainsFB.items.add(juce::FlexItem(symmetryKnob)
		.withMinWidth(50.0f)
		.withMinHeight(100.0f)
		.withFlex(1));

	fb.items.add(juce::FlexItem(gainsFB)
		.withMinWidth(50.0f)
		.withMinHeight(100.0f)
		.withFlex(1));

	fb.items.add(juce::FlexItem(filtersFB)
		.withMinWidth(50.0f)
		.withMinHeight(100.0f)
		.withFlex(1));
	
	auto rect = getLocalBounds();
	rect.reduce(8, 8);
	fb.performLayout(rect);
}
