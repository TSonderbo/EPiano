/*
  ==============================================================================

    EPianoManager.cpp
    Created: 8 Mar 2024 11:46:19am
    Author:  Sonderbo

  ==============================================================================
*/

#include "EPianoManager.h"



EPianoManager::EPianoManager() : tines(TINE_COUNT)
{

}

void EPianoManager::prepareToPlay(double sampleRate)
{
    setCurrentPlaybackSampleRate(sampleRate);

    for (int i = 0; i < TINE_COUNT; i++)
    {
        tines[i].prepareToPlay(sampleRate, i);
    }
}

void EPianoManager::noteAdded(juce::MPENote newNote)
{
    if (newNote.initialNote < MPE_MIN_NOTE || newNote.initialNote > MPE_MAX_NOTE) // Below lowest tine or above highest
        return;

    int noteNumber = newNote.initialNote - 28;

    float velocity = newNote.noteOnVelocity.asUnsignedFloat() * v_max;

    tines[noteNumber].startNote(velocity);
}

void EPianoManager::notePitchbendChanged(juce::MPENote changedNote)
{
    //TODO Dynamic Grid for this
}

void EPianoManager::noteReleased(juce::MPENote finishedNote)
{
    if (finishedNote.initialNote < MPE_MIN_NOTE || finishedNote.initialNote > MPE_MAX_NOTE) // Below lowest tine or above highest
        return;

    int noteNumber = finishedNote.initialNote - MPE_MIN_NOTE;

    tines[noteNumber].stopNote();
}

void EPianoManager::renderNextSubBlock(juce::AudioBuffer<float>& outputAudio, int startSample, int numSamples)
{
    float sample = 0.0f;
    for (int i = 0; i < numSamples; i++)
    {
        sample = renderNextSample();

        for (int channel = 0; channel < outputAudio.getNumChannels(); ++channel)
        {
            outputAudio.addSample(channel, i + startSample, sample);
        }
    }
}

void EPianoManager::setMaximumInputVelocity(float newVelocityMax)
{
    v_max = newVelocityMax;
}

float EPianoManager::renderNextSample()
{
    int activeCount = 0;
    float sample = 0.0f;
    for (auto& tine : tines) {
        if (tine.getIsActive())
        {
            activeCount++;
            sample += tine.processSample();
        }
    }

    if (activeCount == 0)
        return sample;
    else
        return limit(500.0f * sample);
}

float EPianoManager::limit(float sample)
{
    if (sample > 1.0f)
    {
        DBG("Sample Exceeded 1.0");
        return 1.0f;
    }
    else if (sample < -1.0f)
    {
        DBG("Sample Exceeded -1.0");
        return -1.0f;
    }
    else
        return sample;
}