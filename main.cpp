#include <stdio.h>
#include <math.h>

#include "SoundChain.hpp"
#include "PGSSoundChain.hpp"
#include "TBLSoundChain.hpp"
#include "SoundChainWaveFile.hpp"
#include "DAESoundChain.hpp"
#include "MiniAudioSCP.hpp"

#define SAMPLERATE 44100
#define CHANNELCOUNT 1

int main (int argc, char** argv) {

	TBLAdditiveOscSoundChain testOsc;
	testOsc.SetFrequency(220.0f);

	WavReaderSoundChain wavin;
	// wavin.SetPrevious(&testOsc);

	TrimSoundChain trim;
	auto trimParams = trim.GetParameters();
	// trimParams.amplitude = 0.5f;
	trim.SetParameters(trimParams);
	trim.SetPrevious(&testOsc);

	DEA::DSPSoundChain<DEA::DSP_Bypass> deaEffect1;
	auto deaEffect1Parameters = deaEffect1.GetParameters();
	deaEffect1.SetParameters(deaEffect1Parameters);
	deaEffect1.SetPrevious(&trim);

	DEA::DSPSoundChain<DEA::DSP_Bypass> deaEffect2;
	auto deaEffect2Parameters = deaEffect2.GetParameters();
	deaEffect2.SetParameters(deaEffect2Parameters);
	deaEffect2.SetPrevious(&deaEffect1);

	DEA::DSPSoundChain<DEA::DSP_Bypass> deaEffect3;
	auto deaEffect3Parameters = deaEffect3.GetParameters();
	deaEffect3.SetParameters(deaEffect3Parameters);
	deaEffect3.SetPrevious(&deaEffect2);

	DEA::DSPSoundChain<DEA::DSP_Bypass> deaEffect4;
	auto deaEffect4Parameters = deaEffect4.GetParameters();
	deaEffect4.SetParameters(deaEffect4Parameters);
	deaEffect4.SetPrevious(&deaEffect3);

	WavWriterSoundChain wavout;
	wavout.SetPrevious(&deaEffect4);

	SoundChainPlatformSettings soundChainSettings;
	soundChainSettings.SampleRate = SAMPLERATE;
	soundChainSettings.Channels = CHANNELCOUNT;

	MiniAudioSCP platform;
	platform.SetPrevious(&wavout);
	platform.Initialize(soundChainSettings);

	wavout.StartRecording();

	float nextUpdateTime = 0.0f;
	while (platform.GetTime() <= 2.0) {
		float delta = platform.GetTime() / 2.0;

		// deaEffect1Parameters.frequency = delta*delta*delta * 9850.0f + 150;
		// deaEffect1.SetParameters(deaEffect1Parameters);

		if (platform.GetTime() > nextUpdateTime) {
			testOsc.SetFrequency(10000.0f * delta*delta + 10.0f);
			nextUpdateTime += 0.01f;
		}

		// char input;
		// std::cin >> input;
		// if (input == 'q') {
		// 	break;
		// } else if (input == 'p') {
		// 	testOsc.SetFrequency(testOsc.GetFrequency() * 1.0594f);
		// } else if (input == 'o') {
		// 	testOsc.SetFrequency(testOsc.GetFrequency() * 0.9438f);
		// }
	};

	wavout.StopRecording();
	platform.Terminate();
}

// g++ -I./lib/miniaudio main.cpp ./lib/miniaudio/miniaudio.c -o SCWin