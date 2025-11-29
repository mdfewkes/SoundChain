#include <stdio.h>
#include <math.h>

#include "SoundChain.hpp"
#include "SoundChainWaveFile.hpp"
#include "DAESoundChain.hpp"
#include "PortAudioSCP.hpp"
// #include "WaveOutSCP.hpp"

#define SAMPLERATE 44100
#define CHANNELCOUNT 1

int main (int argc, char** argv) {

	WhiteNoiseSoundChain noise;

	WavReaderSoundChain wavin;
	// wavin.SetPrevious(&noise);

	TrimSoundChain trim;
	auto trimParams = trim.GetParameters();
	// trimParams.amplitude = 0.25f;
	trim.SetParameters(trimParams);
	trim.SetPrevious(&wavin);

	DEA::DSPSoundChain<DEA::DSP_Biquad> deaEffect1;
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

	PortAudioSCP platform;
	// WaveOutSCP platform;
	platform.SetPrevious(&wavout);
	platform.Initialize(soundChainSettings);

	wavout.StartRecording();

	while (platform.GetTime() <= 6.0) {
		// printf("%f\n", platform.time());
		float delta = platform.GetTime() / 6.0;

		deaEffect1Parameters.frequency = delta*delta*delta * 9850.0f + 150;
		deaEffect1.SetParameters(deaEffect1Parameters);

		// platform.RecordForSeconds(0.01);
	};

	wavout.StopRecording();
	platform.Terminate();
}