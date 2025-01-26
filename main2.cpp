#include <stdio.h>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include <portaudio.h>

#include "SoundChain.hpp"
#include "SoundChainWaveFile.hpp"

#define SAMPLERATE 44100.
#define CHANNELCOUNT 2

int audio_callback(const void *input, void *output, 
	unsigned long frameCount, 
	const PaStreamCallbackTimeInfo *timeInfo, 
	PaStreamCallbackFlags statuseFlags, 
	void *userData) {

	// SineOsc *osc = (SineOsc *)userData;
	SoundChainBase *data = (SoundChainBase *)userData;
	float *outp = (float *)output;

	data->ReadSamples(outp, frameCount);
	
	return paContinue;
}

static void checkErr(PaError err) {
	if (err != paNoError) {
		printf("PortAudio error: %s", Pa_GetErrorText(err));
		exit(EXIT_FAILURE);
	}
}

int main (int argc, char** argv) {

	SoundChainSettings soundChainSettings;
	soundChainSettings.SampleRate = SAMPLERATE;
	soundChainSettings.Channels = CHANNELCOUNT;

	SinewaveSoundChain* osc1 = new SinewaveSoundChain(330.0);
	SinewaveSoundChain* osc2 = new SinewaveSoundChain(220.0);
	osc2->SetPrevious(osc1);

	WavReaderSoundChain* wavin = new WavReaderSoundChain();
	// wavin->SetPrevious(osc2);
	WavWriterSoundChain* wavout = new WavWriterSoundChain();
	wavout->SetPrevious(wavin);

	TrimSoundChain* out = new TrimSoundChain();
	out->SetAmplitude(1);
	out->SetPrevious(wavout);
	
	out->Initialize(soundChainSettings);

	PaError err;
	PaDeviceIndex outputDevice;
	PaStreamParameters inputParameters, outputParameters;
	PaStream *paStreamHandle;

	err = Pa_Initialize();
	checkErr(err);

	int numDevices = Pa_GetDeviceCount();
	// printf("Number of devices: %d\n", numDevices);
	if (numDevices < 0) {
		printf("Error getting device count.\n");
		exit(EXIT_FAILURE);
	} else if (numDevices == 0) {
		printf("There are no available audio devices on this machine.");
		exit(EXIT_SUCCESS);
	}
	// for (int i = 0; i < numDevices; i++) {
	// 	const PaDeviceInfo *deviceInfo;
	// 	deviceInfo = Pa_GetDeviceInfo(i);
	// 	printf("Device %d:\n", i);
	// 	printf("  name: %s\n", deviceInfo->name);
	// 	printf("  maxInputChannels: %d\n", deviceInfo->maxInputChannels);
	// 	printf("  maxOutputChannels: %d\n", deviceInfo->maxOutputChannels);
	// 	printf("  defaultSampleRate: %f\n", deviceInfo->defaultSampleRate);
	// }

	// printf("choose device for output: ");
	// scanf("%d", &outputDevice);
	outputDevice = Pa_GetDefaultOutputDevice();

	memset(&outputParameters, 0, sizeof(PaStreamParameters));
	outputParameters.device = outputDevice;
	outputParameters.channelCount = CHANNELCOUNT;
	outputParameters.sampleFormat = paFloat32;

	err = Pa_OpenStream(&paStreamHandle, NULL, &outputParameters, 
						SAMPLERATE, 1024, paNoFlag, audio_callback, out);
	checkErr(err);

	wavout->StartRecording();

	err = Pa_StartStream(paStreamHandle);
	checkErr(err);

	Pa_Sleep(2 * 1000);

	wavout->StopRecording();

	err = Pa_StopStream(paStreamHandle);
	checkErr(err);

	err = Pa_CloseStream(paStreamHandle);
	checkErr(err);

	err = Pa_Terminate();
	checkErr(err);
}