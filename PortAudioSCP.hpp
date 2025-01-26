#pragma once

#include <portaudio.h>
#include <stdlib.h>
#include <cstring>

#include "SoundChainPlatform.hpp"

class PortAudioSCP : public SoundChainPlatform {
public:
	PortAudioSCP() {}
	~PortAudioSCP() {}

	static int audio_callback(const void *input, void *output, 
		unsigned long frameCount, 
		const PaStreamCallbackTimeInfo *timeInfo, 
		PaStreamCallbackFlags statuseFlags, 
		void *userData) {

		float *outp = (float *)output;
		PortAudioSCP *scPlatform = (PortAudioSCP *)userData;
		SoundChainBase *chain = scPlatform->GetPrevious();

		chain->ReadSamples(outp, frameCount);
		
		scPlatform->_samplesElapsed += frameCount;

		return paContinue;
	};

private:
	PaStream *_paStreamHandle;

	void Setup() override {
		PaError err;
		PaDeviceIndex outputDevice;
		PaStreamParameters outputParameters;

		err = Pa_Initialize();
		checkErr(err);

		int numDevices = Pa_GetDeviceCount();
		if (numDevices < 0) {
			printf("Error getting device count.\n");
			exit(EXIT_FAILURE);
		} else if (numDevices == 0) {
			printf("There are no available audio devices on this machine.");
			exit(EXIT_SUCCESS);
		}

		outputDevice = Pa_GetDefaultOutputDevice();

		memset(&outputParameters, 0, sizeof(PaStreamParameters));
		outputParameters.device = outputDevice;
		outputParameters.channelCount = _settings.Channels;
		outputParameters.sampleFormat = paFloat32;

		err = Pa_OpenStream(&_paStreamHandle, NULL, &outputParameters, 
							_settings.SampleRate, 0, paNoFlag, PortAudioSCP::audio_callback, this);
		checkErr(err);
	}

	void Start() override {
		PaError err;

		err = Pa_StartStream(_paStreamHandle);
		checkErr(err);
	}

	void End() override {
		PaError err;

		err = Pa_StopStream(_paStreamHandle);
		checkErr(err);

		err = Pa_CloseStream(_paStreamHandle);
		checkErr(err);

		err = Pa_Terminate();
		checkErr(err);
	}

	void checkErr(PaError err) {
		if (err != paNoError) {
			printf("PortAudio error: %s", Pa_GetErrorText(err));
			exit(EXIT_FAILURE);
		}
	}
};