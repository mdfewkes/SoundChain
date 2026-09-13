#pragma once

#include "SoundChainPlatform.hpp"
#include "SoundChainWaveFile.hpp"

class WaveWriteSCP : public SoundChainPlatform {
public:
	WaveWriteSCP() {}
	~WaveWriteSCP() {}

	void RecordForSeconds(double secondsToRecord) {
		int framesToWrite = (int)(secondsToRecord * GetSoundChainSettings().SampleRate);

		while (framesToWrite > _sizeInFrames) {
			FillBuffer(_buffer, _sizeInFrames);
			framesToWrite -= _sizeInFrames;
		}
		FillBuffer(_buffer, framesToWrite);
	}

private:
	WavWriterSoundChain* _wavWriter;
	float* _buffer;
	int _sizeInFrames;

	void Setup() override {
		_sizeInFrames = GetSoundChainSettings().SampleRate * GetSoundChainSettings().Channels;
		_buffer = new float[_sizeInFrames];

		_wavWriter = new WavWriterSoundChain();
		_wavWriter->SetPrevious(GetPrevious());
		SetPrevious(_wavWriter);
	}

	void Start() override {
		printf("Start Recording\n");
		_wavWriter->StartRecording();
	}

	void End() override {
		_wavWriter->StopRecording();
		printf("Stop Recording\n");
		delete _wavWriter;
		delete[] _buffer;
	}
};
