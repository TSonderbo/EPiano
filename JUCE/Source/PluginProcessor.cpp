/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
EPianoAudioProcessor::EPianoAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
	: AudioProcessor(BusesProperties()
#if ! JucePlugin_IsMidiEffect
#if ! JucePlugin_IsSynth
		.withInput("Input", juce::AudioChannelSet::stereo(), true)
#endif
		.withOutput("Output", juce::AudioChannelSet::stereo(), true)
#endif
	), apvts(*this, nullptr, "Parameters", createParams())
#endif
{
	for (int i = 0; i < config::mpe::numVoices; i++)
	{
		synth.addVoice(new Tine());
	}

	synth.setVoiceStealingEnabled(true);

	synth.enableLegacyMode(24);

	updateParams = false;
}

EPianoAudioProcessor::~EPianoAudioProcessor()
{

}

//==============================================================================
const juce::String EPianoAudioProcessor::getName() const
{
	return JucePlugin_Name;
}

bool EPianoAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
	return true;
#else
	return false;
#endif
}

bool EPianoAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
	return true;
#else
	return false;
#endif
}

bool EPianoAudioProcessor::isMidiEffect() const
{
#if JucePlugin_IsMidiEffect
	return true;
#else
	return false;
#endif
}

double EPianoAudioProcessor::getTailLengthSeconds() const
{
	return 0.0;
}

int EPianoAudioProcessor::getNumPrograms()
{
	return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
	// so this should be at least 1, even if you're not really implementing programs.
}

int EPianoAudioProcessor::getCurrentProgram()
{
	return 0;
}

void EPianoAudioProcessor::setCurrentProgram(int index)
{
}

const juce::String EPianoAudioProcessor::getProgramName(int index)
{
	return {};
}

void EPianoAudioProcessor::changeProgramName(int index, const juce::String& newName)
{
}

//==============================================================================
void EPianoAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
	synth.setCurrentPlaybackSampleRate(sampleRate);

	for (int i = 0; i < synth.getNumVoices(); i++)
	{
		if (auto voice = dynamic_cast<Tine*>(synth.getVoice(i)))
		{
			voice->prepareToPlay(sampleRate);
		}
	}

	
	startTimerHz(60);
}

void EPianoAudioProcessor::releaseResources()
{
	// When playback stops, you can use this as an opportunity to free up any
	// spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool EPianoAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
#if JucePlugin_IsMidiEffect
	juce::ignoreUnused(layouts);
	return true;
#else
	// This is the place where you check if the layout is supported.
	// In this template code we only support mono or stereo.
	// Some plugin hosts, such as certain GarageBand versions, will only
	// load plugins that support stereo bus layouts.
	if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
		&& layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
		return false;

	// This checks if the input layout matches the output layout
#if ! JucePlugin_IsSynth
	if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
		return false;
#endif

	return true;
#endif
}
#endif

void EPianoAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
	juce::ScopedNoDenormals noDenormals;
	auto totalNumInputChannels = getTotalNumInputChannels();
	auto totalNumOutputChannels = getTotalNumOutputChannels();

	auto numSamples = buffer.getNumSamples();

	// In case we have more outputs than inputs, this code clears any output
	// channels that didn't contain input data, (because these aren't
	// guaranteed to be empty - they may contain garbage).
	// This is here to avoid people getting screaming feedback
	// when they first compile a plugin, but obviously you don't need to keep
	// this code if your algorithm always overwrites all the output channels.
	for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
		buffer.clear(i, 0, buffer.getNumSamples());

	if (updateParams)
	{
		for (int i = 0; i < synth.getNumVoices(); i++)
		{
			if (auto voice = dynamic_cast<Tine*>(synth.getVoice(i)))
			{
				voice->setParameters(paramValueSet);
			}
		}
		updateParams = false;
	}

	// This is the place where you'd normally do the guts of your plugin's
	// audio processing...
	// Make sure to reset the state if your inner loop is processing
	// the samples and the outer loop is handling the channels.
	// Alternatively, you can process the samples with the channels
	// interleaved by keeping the same state.

	synth.renderNextBlock(buffer, midiMessages, 0, numSamples);

	scopeDataCollector.process(buffer.getReadPointer(0), (size_t)numSamples);
}

//==============================================================================
bool EPianoAudioProcessor::hasEditor() const
{
	return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* EPianoAudioProcessor::createEditor()
{
	return new EPianoAudioProcessorEditor(*this);
}

//==============================================================================
void EPianoAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
	// You should use this method to store your parameters in the memory block.
	// You could do that either as raw data, or use the XML or ValueTree classes
	// as intermediaries to make it easy to save and load complex data.

	auto state = apvts.copyState();
	std::unique_ptr<juce::XmlElement> xml(state.createXml());
	copyXmlToBinary(*xml, destData);
}

void EPianoAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
	// You should use this method to restore your parameters from this memory block,
	// whose contents will have been created by the getStateInformation() call.

	std::unique_ptr<juce::XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));
	if (xmlState != nullptr) {
		if (xmlState->hasTagName(apvts.state.getType())) {
			apvts.replaceState(juce::ValueTree::fromXml(*xmlState));
		}
	}
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
	return new EPianoAudioProcessor();
}

void EPianoAudioProcessor::timerCallback()
{
	updateParams |= checkParameterValues();
}

juce::AudioProcessorValueTreeState::ParameterLayout EPianoAudioProcessor::createParams()
{
	std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;
	using namespace config::parameter;

	//Inserts audio params here

	//Tone control
	params.push_back(std::make_unique<juce::AudioParameterFloat>(id_pickup_lowpass_cutoff, name_pickup_lowpass_cutoff, juce::NormalisableRange<float>(200.0f, 2000.0f), config::pickup::lpCutoff));
	params.push_back(std::make_unique<juce::AudioParameterFloat>(id_pickup_lowpass_resonance, name_pickup_lowpass_resonance, juce::NormalisableRange<float>(0.1f, 1.0f), config::pickup::resonance));
	params.push_back(std::make_unique<juce::AudioParameterFloat>(id_pickup_highpass_resonance, name_pickup_highpass_resonance, juce::NormalisableRange<float>(0.1f, 1.0f), config::pickup::resonance));
	params.push_back(std::make_unique<juce::AudioParameterFloat>(id_pickup_gain, name_pickup_gain, juce::NormalisableRange<float>(0.1f, 30.0f), config::pickup::gain));
	params.push_back(std::make_unique<juce::AudioParameterFloat>(id_pickup_symmetry, name_pickup_symmetry, juce::NormalisableRange<float>(0.1f, 30.0f), config::pickup::symmetry));
	params.push_back(std::make_unique<juce::AudioParameterBool>(id_pickup_bypass, name_pickup_bypass, false));
	

	//Volume Controls
	params.push_back(std::make_unique<juce::AudioParameterFloat>(id_master_volume, name_master_volume, juce::NormalisableRange<float>(0.1f, 2.0f), config::master_volume));
	params.push_back(std::make_unique<juce::AudioParameterFloat>(id_tine_gain, name_tine_gain, juce::NormalisableRange<float>(1.0f, 20000.0f), config::tine::tineGain));


	return { params.begin(), params.end() };
}

AudioBufferQueue& EPianoAudioProcessor::getAudioBufferQueue()
{
	return audioBufferQueue;
}

bool EPianoAudioProcessor::checkParameterValues()
{
	bool paramChanged = false;

	using namespace config::parameter;
	//Tone control
	paramChanged |= paramValueSet.set(id_pickup_lowpass_cutoff, apvts.getRawParameterValue(id_pickup_lowpass_cutoff)->load());
	paramChanged |= paramValueSet.set(id_pickup_lowpass_resonance, apvts.getRawParameterValue(id_pickup_lowpass_resonance)->load());
	paramChanged |= paramValueSet.set(id_pickup_highpass_resonance, apvts.getRawParameterValue(id_pickup_highpass_resonance)->load());
	paramChanged |= paramValueSet.set(id_pickup_gain, apvts.getRawParameterValue(id_pickup_gain)->load());
	paramChanged |= paramValueSet.set(id_pickup_symmetry, apvts.getRawParameterValue(id_pickup_symmetry)->load());
	paramChanged |= paramValueSet.set(id_pickup_bypass, apvts.getRawParameterValue(id_pickup_bypass)->load());
	//Volume Controls
	paramChanged |= paramValueSet.set(id_master_volume, apvts.getRawParameterValue(id_master_volume)->load());
	paramChanged |= paramValueSet.set(id_tine_gain, apvts.getRawParameterValue(id_tine_gain)->load());

	return paramChanged;
}