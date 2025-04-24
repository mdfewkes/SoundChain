#include <stdio.h>
#include <math.h>

#include "SoundChain.hpp"
#include "SoundChainWaveFile.hpp"
#include "DAESoundChain.hpp"
// #include "PortAudioSCP.hpp"
#include "WaveOutSCP.hpp"

#define SAMPLERATE 44100
#define CHANNELCOUNT 1

int main (int argc, char** argv) {

	SoundChainSettings soundChainSettings;
	soundChainSettings.SampleRate = SAMPLERATE;
	soundChainSettings.Channels = CHANNELCOUNT;

	WhiteNoiseSoundChain noise;

	WavReaderSoundChain wavin;
	// wavin.SetPrevious(&noise);

	TrimSoundChain trim;
	auto trimParams = trim.GetParameters();
	// trimParams.amplitude = 0.25f;
	trim.SetParameters(trimParams);
	trim.SetPrevious(&wavin);

	DEA::DSPSoundChain<DEA::DSP_AnalogFIRFilter> deaEffect1;
	auto deaEffect1Parameters = deaEffect1.GetParameters();
	deaEffect1.SetParameters(deaEffect1Parameters);
	deaEffect1.SetPrevious(&trim);

	DEA::DSPSoundChain<DEA::DSP_3BandCompressor> deaEffect2;
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

	// WavWriterSoundChain wavout;
	// wavout.SetPrevious(&deaEffect4);

	// PortAudioSCP platform;
	WaveOutSCP platform;
	platform.SetPrevious(&deaEffect4);
	platform.Initialize(soundChainSettings);

	// wavout.StartRecording();

	while (platform.time() <= 6) {
		float delta = platform.time() / 6.0;

		deaEffect1Parameters.frequency = delta * 10000.0f + 1;
		deaEffect1.SetParameters(deaEffect1Parameters);

		platform.RecordForSeconds(0.01);
	};

	// wavout.StopRecording();
	platform.Terminate();
}