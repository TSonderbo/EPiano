/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "EPianoManager.h"
#include "ScopeDataCollector.h"
#include "AudioBufferQueue.h"
#include "Definitions.h"

//==============================================================================
/**
*/
class EPianoAudioProcessor  : public juce::AudioProcessor
                            #if JucePlugin_Enable_ARA
                             , public juce::AudioProcessorARAExtension
                            #endif
{
public:
    //==============================================================================
    EPianoAudioProcessor();
    ~EPianoAudioProcessor() override;

    juce::AudioProcessorValueTreeState apvts;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    AudioBufferQueue& getAudioBufferQueue();

private:

    EPianoManager ePianoManager;
    AudioBufferQueue audioBufferQueue;
    ScopeDataCollector scopeDataCollector{ audioBufferQueue };

    juce::MidiBuffer uiMidiBuffer;


    //==============================================================================

    juce::AudioProcessorValueTreeState::ParameterLayout createParams();

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (EPianoAudioProcessor)
};
