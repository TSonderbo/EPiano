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
		unsigned const int minNote = 28;
		unsigned const int maxNote = 99;
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
	}
}

