/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
EPianoAudioProcessorEditor::EPianoAudioProcessorEditor (EPianoAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p), oscilloscope(audioProcessor.getAudioBufferQueue())
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.

    addAndMakeVisible(oscilloscope);

    setSize(1280, 720);
}

EPianoAudioProcessorEditor::~EPianoAudioProcessorEditor()
{
    
}

//==============================================================================
void EPianoAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
}

void EPianoAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..

    auto visualiserCompWidth = 2800;
    auto visualiserCompHeight = 300;

    juce::FlexBox main;
    main.flexDirection = juce::FlexBox::Direction::column;

    main.items.add(juce::FlexItem(oscilloscope)
        .withMinWidth(200.0f)
        .withMinHeight(200.0f)
        .withFlex(1.0f));

    main.performLayout(getLocalBounds().toFloat());
}