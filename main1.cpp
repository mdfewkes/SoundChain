#include <stdio.h>
#include <time.h>
#include <cstring>
#include <stdlib.h>
#include <portaudio.h>

#define SAMPLERATE 48000.
#define BUFFERSIZE 1024

typedef struct _myUserData {
	float delay[(int)SAMPLERATE/2];
	int rp;
} myUserData;

int audio_callback(const void *input, void *output, 
	unsigned long frameCount, 
	const PaStreamCallbackTimeInfo *timeInfo, 
	PaStreamCallbackFlags statuseFlags, 
	void *userData) {

	myUserData *p = (myUserData *)userData;
	int rp = p->rp;
	float out, *delay = p->delay;
	float *inp = (float *)input, *outp = (float *)output;
	for (int i = 0; i < frameCount; i++) {
		out = delay[rp];
		delay[rp++] = inp[i] + out*0.5;
		if (rp >= SAMPLERATE/2) rp = 0;
		outp[i] = out + inp[i];
	}
	p->rp = rp;
	return paContinue;
}

static void checkErr(PaError err) {
	if (err != paNoError) {
		printf("PortAudio error: %s", Pa_GetErrorText(err));
		exit(EXIT_FAILURE);
	}
}

int main (int argc, char** argv) {

	PaError err;
	PaDeviceIndex inputDevice, outputDevice;
	PaStreamParameters inputParameters, outputParameters;
	PaStream *paStreamHandle;
	myUserData *data = (myUserData *) calloc(sizeof(myUserData), 1);

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

	// printf("choose device for input: ");
	// scanf("%d", &inputDevice);
	inputDevice = 0; //Pa_GetDefaultInputDevice();
	// printf("choose device for output: ");
	// scanf("%d", &outputDevice);
	outputDevice = 0; //Pa_GetDefaultOutputDevice();

	memset(&inputParameters, 0, sizeof(PaStreamParameters));
	inputParameters.device = inputDevice;
	inputParameters.channelCount = 1;
	inputParameters.sampleFormat = paFloat32;
	memset(&outputParameters, 0, sizeof(PaStreamParameters));
	outputParameters.device = outputDevice;
	outputParameters.channelCount = 1;
	outputParameters.sampleFormat = paFloat32;

	err = Pa_OpenStream(&paStreamHandle, &inputParameters, &outputParameters, 
						SAMPLERATE, BUFFERSIZE, paNoFlag, audio_callback, data);
	checkErr(err);

	err = Pa_StartStream(paStreamHandle);
	checkErr(err);

	Pa_Sleep(10 * 1000);

	err = Pa_StopStream(paStreamHandle);
	checkErr(err);

	err = Pa_CloseStream(paStreamHandle);
	checkErr(err);

	err = Pa_Terminate();
	checkErr(err);

	free(data);

	return EXIT_SUCCESS;
}