/*
  ==============================================================================

	Configuration.h
	Created: 8 Mar 2024 10:36:28am
	Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

namespace config
{
	//Universal
	unsigned const int oversampling = 32;

	//Tines
	namespace tine
	{
		unsigned const int count = 73;
		const float inactiveTimer = 0.2f;
	}
	
	//MPE-MIDI
	namespace mpe
	{
		unsigned const int numVoices = 3;
		unsigned const int minNote = 28;
		unsigned const int maxNote = 100;
		const float maxInputVelocity = 4.0f;
	}

	namespace pickup
	{
		const float lpCutoff = 1000.0f;
		const float resonance = 0.404061f;
		const float smoothDuration = 0.2f;
		const float gain = 15.0f;
		const float symmetry = 15.0f;
	}

	//Configuration
	namespace parameter
	{
		//Input velocity
		const juce::String id_velocity("velocity");
		const juce::String name_velocity("Velocity");
		
		//Amplitude
		const juce::String id_amplitude("amplitude");
		const juce::String name_amplitude("Amplitude");

		//Magnetic pickup lowpass cutoff
		const juce::String id_pickup_lowpass_cutoff("tone_cutoff");
		const juce::String name_pickup_lowpass_cutoff("Tone");

		//Magnetic pickup lowpass resonance
		const juce::String id_pickup_lowpass_resonance("tone_resonance_1");
		const juce::String name_pickup_lowpass_resonance("Resonance 1");

		//Magnetic pickup highpass resonance
		const juce::String id_pickup_highpass_resonance("tone_resonance_2");
		const juce::String name_pickup_highpass_resonance("Resonance 2");

		//Magnetic pickup highpass resonance
		const juce::String id_pickup_gain("tone_gain");
		const juce::String name_pickup_gain("Gain");

		//Magnetic pickup highpass resonance
		const juce::String id_pickup_symmetry("tone_symmetry");
		const juce::String name_pickup_symmetry("Symmetry");

		//Magnetic pickup bypass
		const juce::String id_pickup_bypass("tone_bypass");
		const juce::String name_pickup_bypass("Tone bypass");
	}
}

