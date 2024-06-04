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
	unsigned const int master_volume = 1.0f;


	//Tines
	namespace tine
	{
		unsigned const int count = 73;
		const float inactiveTimer = 0.2f;
		unsigned const int tineGain = 10000.0f;
		const float smoothDuration = 0.2f;
	}
	
	//MPE-MIDI
	namespace mpe
	{
		unsigned const int numVoices = 1;
		unsigned const int minNote = 28;
		unsigned const int maxNote = 100;
		const float maxInputVelocity = 6.0f;
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

#pragma region pickup
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
#pragma endregion 

		//Master Volume
		const juce::String id_master_volume("master_volume");
		const juce::String name_master_volume("Master");

		//Pre Pickup Volume
		const juce::String id_tine_gain("tine_gain");
		const juce::String name_tine_gain("Tine");
	}
}

