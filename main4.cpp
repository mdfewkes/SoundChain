#include <stdio.h>
#include <math.h>

#include "SoundChain.hpp"
#include "SoundChainWaveFile.hpp"
#include "DAESoundChain.hpp"
#include "PortAudioSCP.hpp"
#include "WaveOutSCP.hpp"

#define SAMPLERATE 44100
#define CHANNELCOUNT 2

int main (int argc, char** argv) {

	SoundChainSettings soundChainSettings;
	soundChainSettings.SampleRate = SAMPLERATE;
	soundChainSettings.Channels = CHANNELCOUNT;

	// WhiteNoiseSoundChain noise;

	WavReaderSoundChain wavin;
	// wavin.SetPrevious(&noise);

	TrimSoundChain trim;
	auto trimParams = trim.GetParameters();
	trim.SetPrevious(&wavin);

	DelaySoundChain delay;
	delay.SetPrevious(&trim);

	WavWriterSoundChain wavout;
	wavout.SetPrevious(&delay);

	PortAudioSCP platform;
	// WaveOutSCP platform;
	platform.SetPrevious(&wavout);
	platform.Initialize(soundChainSettings);

	wavout.StartRecording();

	trimParams.amplitude = 1.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.05) {
		printf("%s: %ld\n", "Samples: ",  platform.sampleElapsed());
		printf("%s: %f\n", "Seconds: ",  platform.time());
	};
	trimParams.amplitude = 0.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.25) {};

	trimParams.amplitude = 1.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.3) {};
	trimParams.amplitude = 0.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.5) {};

	trimParams.amplitude = 1.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.55) {};
	trimParams.amplitude = 0.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.75) {};

	trimParams.amplitude = 1.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.8) {};
	trimParams.amplitude = 0.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 0.95) {};

	trimParams.amplitude = 1.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 1.0) {};
	trimParams.amplitude = 0.0f;
	trim.SetParameters(trimParams);
	while (platform.time() < 1.5) {
		printf("%s: %ld\n", "Samples: ",  platform.sampleElapsed());
		printf("%s: %f\n", "Seconds: ",  platform.time());
	};

	wavout.StopRecording();
	platform.Terminate();
}