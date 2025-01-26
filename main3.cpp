#include "SoundChain.hpp"
#include "SoundChainWaveFile.hpp"
#include "SoundChainBuffers.hpp"

#define SAMPLERATE 44100
#define CHANNELCOUNT 2
#define BUFFERSIZE 1024


int main(int argc, char** argv) {
	SoundChainSettings soundChainSettings;
	soundChainSettings.SampleRate = SAMPLERATE;
	soundChainSettings.Channels = CHANNELCOUNT;

	auto* osc = new SinewaveSoundChain(440.0);
	//auto* reader = new WavReaderSoundChain();
	auto* doubleBuffer = new MultiBufferSoundChain();
	auto* writer = new WavWriterSoundChain();
	doubleBuffer->SetBufferInput(osc);
	writer->SetPrevious(doubleBuffer);
	
	writer->Initialize(soundChainSettings);

	float* buffer = new float[BUFFERSIZE * CHANNELCOUNT]; 

	writer->StartRecording();
	int numberOfFrames = BUFFERSIZE;
	for (int samplesToWrite = SAMPLERATE * 2; samplesToWrite > 0; samplesToWrite -= BUFFERSIZE) {
		if (samplesToWrite < BUFFERSIZE) {
			numberOfFrames = samplesToWrite;
		}

		writer->ReadSamples(buffer, numberOfFrames);
		// osc->SetFrequency(osc->GetFrequency() * 1.02);
	}
	writer->StopRecording();

	return 0;
}