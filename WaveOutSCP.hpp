#pragma once

#include "SoundChainPlatform.hpp"
#include "SoundChainWaveFile.hpp"

class WaveOutSCP : public SoundChainPlatform {
public:
	WaveOutSCP() {}
	~WaveOutSCP() {}

	void RecordForSeconds(double secondsToRecord) {
		int framesToWrite = (int)(secondsToRecord * _settings.SampleRate);

		while (framesToWrite > _sizeInFrames) {
			_wavWriter->ReadSamples(_buffer, _sizeInFrames);
			framesToWrite -= _sizeInFrames;
		}
		_wavWriter->ReadSamples(_buffer, framesToWrite);

		_samplesElapsed += framesToWrite;
	}

private:
	WavWriterSoundChain* _wavWriter;
	float*  _buffer;
	int _sizeInFrames;

	void Setup() override {
		_sizeInFrames = _settings.SampleRate;
		_buffer = new float[_sizeInFrames * _settings.Channels];

		_wavWriter = new WavWriterSoundChain();
		_wavWriter->Initialize(_settings);
		_wavWriter->SetPrevious(GetPrevious());
	}

	void Start() override {
		_wavWriter->StartRecording();
	}

	void End() override {
		_wavWriter->StopRecording();
		delete[] _buffer;
	}
};